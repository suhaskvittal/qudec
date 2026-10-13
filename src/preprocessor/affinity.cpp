/*
 *  author: Suhas Vittal
 *  contributor: OpenAI Codex
 *  date:   24 September 2026
 * */

#include "preprocessor/affinity.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace pp
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

using dem_type = stim::DetectorErrorModel;

AffinityTable::value_type
_quantize(int distance)
{
    // One bucket per hop; reserve the largest value for saturated/disconnected paths.
    return static_cast<AffinityTable::value_type>(std::min(distance, int(AffinityTable::MASK)));
}

size_t
_data_width(size_t num_detectors)
{
    constexpr size_t CELLS_PER_WORD = 8 / AffinityTable::BIT_WIDTH;
    if (num_detectors == std::numeric_limits<size_t>::max())
        throw std::length_error("affinity table is too large");
    const size_t n = num_detectors + 1; // Include the shared boundary.
    if (n > std::numeric_limits<size_t>::max() / n)
        throw std::length_error("affinity table is too large");
    const size_t cells = n * n;
    return cells / CELLS_PER_WORD + (cells % CELLS_PER_WORD != 0);
}

void
_read(gzFile& file, void* data, size_t size)
{
    auto* bytes = static_cast<unsigned char*>(data);
    while (size)
    {
        const unsigned chunk = static_cast<unsigned>(std::min(size, size_t(1 << 20)));
        if (gzread(file, bytes, chunk) != static_cast<int>(chunk))
            throw std::runtime_error("unable to read affinity table");
        bytes += chunk;
        size -= chunk;
    }
}

void
_write(gzFile& file, const void* data, size_t size)
{
    const auto* bytes = static_cast<const unsigned char*>(data);
    while (size)
    {
        const unsigned chunk = static_cast<unsigned>(std::min(size, size_t(1 << 20)));
        if (gzwrite(file, bytes, chunk) != static_cast<int>(chunk))
            throw std::runtime_error("unable to write affinity table");
        bytes += chunk;
        size -= chunk;
    }
}

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

AffinityTable::AffinityTable(const dem_type& dem)
    :num_detectors_(dem.count_detectors()),
    data_(_data_width(num_detectors_))
{
    aff::hg_type gr = aff::build_hypergraph_from_dem(dem);

    // Run floyd-warshall to get all pairs distances
    const size_t n = gr.N();
    std::vector<int> dist(n*n, int(AffinityTable::MASK));

    // define `dist` to be `1` wherever two nodes share a hyperedge
    for (hg::id_type v = 0; v < n; v++)
        dist[v*n+v] = 0;

    for (hg::id_type e = 0; e < gr.M(); e++)
    {
        auto supp = gr.support(e);
        for (auto v : supp)
            for (auto w : supp)
                if (v != w && !gr.v(v).is_boundary)  // do not set dist[B][*] = 1 
                    dist[v*n+w] = 1;
    }

    for (size_t k = 0; k < n; k++)
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                dist[i*n+j] = std::min(dist[i*n+j], dist[i*n+k] + dist[k*n+j]);
    

    // now we have the distances, we can populate `data_`
    for (size_t i = 0; i < n; i++)
        for (size_t j = i+1; j < n; j++)
            set( i, j, _quantize(dist[i*n+j]) );
}

AffinityTable::AffinityTable(gzFile& f)
    :num_detectors_{0}
{
    uint32_t count;
    _read(f, &count, sizeof(count));
    num_detectors_ = count;
    data_.resize(_data_width(num_detectors_));
    _read(f, data_.data(), data_.size() * sizeof(value_type));
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

AffinityTable::value_type
AffinityTable::operator()(size_t i, size_t j) const
{
    if (i > num_detectors_ || j > num_detectors_)
        throw std::out_of_range("affinity table detector ID is out of range");
    if (i == j)
        return 0;
    if (i > j)
        std::swap(i, j);
    const auto [word_idx, bit_idx] = index(i,j);
    return (static_cast<uint8_t>(data_[word_idx]) >> bit_idx) & MASK;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
AffinityTable::write_to_file(gzFile& f)
{
    if (num_detectors_ > std::numeric_limits<uint32_t>::max())
        throw std::length_error("affinity table detector count exceeds file format");
    const uint32_t count = static_cast<uint32_t>(num_detectors_);
    _write(f, &count, sizeof(count));
    _write(f, data_.data(), data_.size() * sizeof(value_type));
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
AffinityTable::set(size_t i, size_t j, value_type x)
{
    auto [word_idx, bit_idx] = index(i,j);
    const auto word = static_cast<uint8_t>(data_[word_idx]);
    data_[word_idx] = static_cast<value_type>(
            (word & ~(MASK << bit_idx)) | ((x & MASK) << bit_idx));
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::tuple<size_t, size_t>
AffinityTable::index(size_t i, size_t j) const
{
    constexpr size_t BITS_PER_WORD{sizeof(value_type) * 8};
    constexpr size_t CELLS_PER_WORD{BITS_PER_WORD / BIT_WIDTH};
    const size_t n = num_detectors() + 1;
    const size_t cell = i * n + j;                 // linear cell index
    return { cell / CELLS_PER_WORD,                // word index in data_
             (cell % CELLS_PER_WORD) * BIT_WIDTH }; // bit offset within word
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace aff
{

// OpenAI GPT-6: Convert complete DEM faults into parity-reduced hyperedges,
// attaching one shared boundary to odd supports and merging independent faults
// with identical detector support by their odd-parity probability.
hg_type
build_hypergraph_from_dem(const dem_type& dem)
{
    const uint64_t num_detectors = dem.count_detectors();
    if (num_detectors > static_cast<uint64_t>(std::numeric_limits<hg::id_type>::max()))
        throw std::runtime_error("too many detectors for Hypergraph IDs");

    hg_type gr;
    for (uint64_t i = 0; i < num_detectors; i++)
        gr.add_vertex(v_data_type{});
    const auto boundary_id = gr.add_vertex(v_data_type{true});

    std::map<std::vector<hg::id_type>, double> merged;
    dem.iter_flatten_error_instructions(
            [&] (const stim::DemInstruction& instruction)
            {
                const double pr = instruction.arg_data[0];
                if (!std::isfinite(pr) || pr < 0.0 || pr > 1.0)
                    throw std::runtime_error("invalid DEM error probability");
                if (pr == 0.0)
                    return;

                std::vector<hg::id_type> support;
                // Logical targets do not affect detector affinity. Separators
                // describe components of the same fault, not independent errors.
                for (const auto& target : instruction.target_data)
                    if (target.is_relative_detector_id())
                    {
                        if (target.val() >= num_detectors)
                            throw std::runtime_error("DEM detector target is out of range");
                        support.push_back(static_cast<hg::id_type>(target.val()));
                    }

                std::sort(support.begin(), support.end());
                std::vector<hg::id_type> odd_support;
                for (size_t i = 0; i < support.size();)
                {
                    size_t j = i + 1;
                    while (j < support.size() && support[j] == support[i])
                        j++;
                    if ((j - i) % 2 != 0)
                        odd_support.push_back(support[i]);
                    i = j;
                }
                if (odd_support.empty())
                    return;
                if (odd_support.size() % 2 != 0)
                    odd_support.push_back(boundary_id);
                if (odd_support.size() > HG_ORDER)
                    throw std::runtime_error("DEM error support exceeds HG_ORDER after boundary insertion");

                auto [it, inserted] = merged.try_emplace(std::move(odd_support), pr);
                if (!inserted)
                    it->second = it->second * (1.0 - pr) + pr * (1.0 - it->second);
            });

    for (const auto& [support, pr] : merged)
        if (pr > 0.0)
            gr.add_edge(support, e_data_type{pr});
    return gr;
}

} // namespace aff

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace pp
