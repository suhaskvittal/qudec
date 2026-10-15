/*
 *  author: Suhas Vittal
 *  date:   15 October 2026
 * */

#define PHANTOM_DEBUG

#include "decoder/phantom.h"

#include <algorithm>
#include <iomanip>
#include <numeric>
#include <stdexcept>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

constexpr size_t ASTREA_MAX_HW{10};

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Phantom::Phantom(const stim::DetectorErrorModel& dem, size_t code_distance)
    : code_distance(code_distance), d(code_distance), dem(dem),
      num_detectors(dem.count_detectors()), num_observables(dem.count_observables()),
      internal_matcher_(std::make_unique<PyMatching>(dem)), aff_(dem)
{
    if (code_distance == 0)
        throw std::invalid_argument("Phantom requires a positive code distance");
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Result
Phantom::decode(SyndromeRef s, ObsRef obs)
{
#if defined(PHANTOM_DEBUG)
    std::cout << "---------------------------------------------------------\n";
#endif

    detectors_.clear();
    is_matched_.clear();
    pairwise_affinity_.clear();

    detectors_.reserve(s.popcnt());
    for (det_id_type i = 0; i < num_detectors; i++)
        if (s[i])
            detectors_.push_back(i);
    const size_t hw = detectors_.size();
    is_matched_.resize(num_detectors, false);

    // build `pairwise_affinity_`
    pairwise_affinity_.resize(hw*hw);
    for (size_t i = 0; i < hw; i++)
    {
        for (size_t j = 0; j < hw; j++)
        {
            size_t idx = i*hw + j;
            if (i == j) pairwise_affinity_[idx] = aff_type{0};
            else        pairwise_affinity_[idx] = aff_(detectors_[i], detectors_[j]);
        }
    }

#if defined(PHANTOM_DEBUG)
    std::cout << "dets = ";
    for (size_t d : detectors_)
        std::cout << " " << d;
    std::cout << "\naffinity matrix:\n";
    for (size_t i = 0; i < hw; i++)
    {
        std::cout << detectors_[i] << "\t" << std::hex;
        for (size_t j = 0; j < hw; j++)
            std::cout << " " << (pairwise_affinity_[i*hw+j]+0);
        std::cout << std::dec << "\n";
    }
#endif

    // do decode:
    Result out;
    out.flipped_obs = ObsType(num_observables);

    for (size_t w = 1; w <= (d-1)/2 && remaining() > ASTREA_MAX_HW; w++)
    {
        size_t rem_b = remaining();
        out.flipped_obs ^= filter(w).flipped_obs;
        if (remaining() == rem_b)
            break;
    }

    const size_t max_aff_start = std::max(size_t{1}, (d-1)/2);
    for (size_t a = max_aff_start; a > 0 && remaining() > ASTREA_MAX_HW; a--)
    {
        auto clusters = build_clusters(a);
        if (clusters.size() > 0)
        {
            // run pymatching on each cluster and accumulate results.
            for (const auto& cluster : clusters)
            {
                SyndromeType cs(num_detectors);
                for (auto i : cluster) 
                    cs[i] = 1;
                auto sub_result = internal_matcher_->decode(cs, obs);
                out.flipped_obs ^= sub_result.flipped_obs;

                for (auto i : cluster)
                    is_matched_[i] = true;

#if defined(PHANTOM_DEUBG)
                std::cout << "(a = " << a < << ") cluster: ";
                for (size_t i : cluster)
                    std::cout << " " << i;
                std::cout << ", obs = ";
                for (size_t i = 0; i < num_observables; i++)
                    if (sub_result.flipped_obs[i])
                        std::cout << " " << i;
                std::cout << ", remaining = ";
                for (size_t d : detectors_)
                    if (!is_matched_[d])
                        std::cout << " " << d;
                std::cout << "\n";
#endif
            }
        }
    }

    // complete any remaining detectors
    if (remaining() > 0)
    {
        SyndromeType rem_s(num_detectors);
        for (size_t d : detectors_)
            if (!is_matched_[d])
                rem_s[d] ^= 1;
        auto rem_out = internal_matcher_->decode(rem_s, obs);
        out.flipped_obs ^= rem_out.flipped_obs;
    }

#if defined(PHANTOM_DEBUG)
    auto ref_out = internal_matcher_->decode(s, obs);
    bool error{false};
    for (size_t i = 0; i < num_observables; i++)
        error = error || (out.flipped_obs[i] != ref_out.flipped_obs[i]);
    if (error)
        std::cout << "logical error\n";
#endif

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
Phantom::print_stats(std::ostream& out) const
{
    (void)out;
}

void
Phantom::mpi_accumulate()
{
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Result
Phantom::filter(size_t w)
{
    const size_t hw = detectors_.size();

    std::vector<size_t> unique_neighbor(hw),
                        count(hw, 0);
    for (size_t i = 0; i < hw; i++)
    {
        if (is_matched_[detectors_[i]])
            continue;
        for (size_t j = 0; j < hw; j++)
        {
            if (i == j || is_matched_[detectors_[j]])
                continue;
            size_t idx = i*hw + j;
            if (pairwise_affinity_[idx] <= w)
            {
                unique_neighbor[i] = j;
                count[i]++;
            }
        }
    }

    SyndromeType s(num_detectors);
    [[ maybe_unused ]] const size_t b = num_detectors;
    for (size_t i = 0; i < hw; i++)
    {
        if (is_matched_[detectors_[i]])
            continue;
        if (count[i] == 1 && count[unique_neighbor[i]] == 1)
        {
            s[detectors_[i]] ^= 1;
            is_matched_[detectors_[i]] = true; 

#if defined(PHANTOM_DEBUG)
            std::cout << "filter: pairing " << detectors_[i] << " and " << detectors_[unique_neighbor[i]] << "\n";
#endif
        }
    }

    if (s.popcnt() == 0)
        return Result{};

    ObsType dummy(num_observables);
    auto out = internal_matcher_->decode(s, dummy);
#if defined(PHANTOM_DEBUG)
    std::cout << "filter flipped obs: ";
    for (size_t i = 0; i < num_observables; i++)
        if (out.flipped_obs[i])
            std::cout << i;
    std::cout << ", remaining dets =";
    for (size_t i = 0; i < hw; i++)
        if (!is_matched_[detectors_[i]])
            std::cout << " " << detectors_[i];
    std::cout << "\n";
#endif
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::vector<std::vector<Phantom::det_id_type>>
Phantom::build_clusters(size_t max_aff)
{
    const size_t hw = detectors_.size();
    std::vector<size_t> parent(hw);
    std::iota(parent.begin(), parent.end(), size_t{0});

    auto root = [&](size_t i)
    {
        while (parent[i] != i)
        {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };

    // Examine every unordered pair, including opposite pairs when hw is even.
    // Join roots so an edge never overwrites an earlier connection.
    for (size_t i = 0; i < hw; i++)
    {
        if (is_matched_[detectors_[i]])
            continue;
        for (size_t j = i + 1; j < hw; j++)
        {
            if (is_matched_[detectors_[j]])
                continue;
            if (pairwise_affinity_[i*hw + j] <= max_aff)
            {
                const size_t ri = root(i), rj = root(j);
                parent[std::max(ri, rj)] = std::min(ri, rj);
            }
        }
    }

    std::vector<std::vector<det_id_type>> out(hw);
    for (size_t i = 0; i < hw; i++)
    {
        if (is_matched_[detectors_[i]])
            continue;
        out[root(i)].push_back(detectors_[i]);
    }
    auto it = std::remove_if(out.begin(), out.end(), 
                    [] (const auto& c) { return c.empty() || c.size() > ASTREA_MAX_HW; });
    out.erase(it, out.end());
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec
