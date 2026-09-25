/*
 * author: OpenAI GPT-6
 * date: 24 September 2026
 * purpose: Record SI1000 surface-code syndromes, PyMatching outcomes, and
 *          effective detector-partner counts in a compact binary format.
 */

#include "circuit_generator.h"
#include "decoder/surface_code.h"
#include "hypergraph.h"
#include "preprocessor/affinity.h"

#include <stim.h>
#include <stim/util_top/circuit_to_dem.h>

#include <argparse/argparse.h>
#include <lzma.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace
{

struct VertexData { };
struct EdgeData { double pr; };

// Reject raw DEM terms exceeding the graph's supported hyperedge order.
constexpr size_t MAX_ERROR_ORDER{16};
using AffinityGraph = Hypergraph<VertexData, EdgeData, MAX_ERROR_ORDER>;
constexpr size_t BATCH_SIZE{64};
constexpr std::string_view FILE_MAGIC{"QDAFF001"};

// OpenAI GPT-6: Encode the versioned file in little-endian fixed-width fields.
// The flags byte uses bit 0 for the odd-weight boundary rule and bit 1 for
// opposite-basis detector inclusion.
void append_u8(std::string& out, uint8_t value) { out.push_back(static_cast<char>(value)); }

void
append_u16(std::string& out, uint16_t value)
{
    for (size_t i = 0; i < 2; ++i)
        append_u8(out, static_cast<uint8_t>(value >> (8 * i)));
}

void
append_u64(std::string& out, uint64_t value)
{
    for (size_t i = 0; i < 8; ++i)
        append_u8(out, static_cast<uint8_t>(value >> (8 * i)));
}

void
append_f64(std::string& out, double value)
{
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
    append_u64(out, std::bit_cast<uint64_t>(value));
}

// OpenAI GPT-6: Compress binary batches as they arrive when the path ends in
// .xz, without keeping earlier batches in memory.
class BinaryWriter
{
public:
    explicit BinaryWriter(const std::string& path)
        : output_(path, std::ios::binary | std::ios::trunc),
          compressed_(path.ends_with(".xz"))
    {
        if (!output_)
            throw std::runtime_error("failed to open output file: " + path);
        if (compressed_)
        {
            const auto status = lzma_easy_encoder(&stream_, 1, LZMA_CHECK_CRC64);
            if (status != LZMA_OK)
                throw std::runtime_error("failed to initialize LZMA encoder: " + std::to_string(status));
        }
    }

    BinaryWriter(const BinaryWriter&) = delete;
    BinaryWriter& operator=(const BinaryWriter&) = delete;

    ~BinaryWriter()
    {
        if (compressed_)
            lzma_end(&stream_);
    }

    void write(std::string_view data)
    {
        if (finished_)
            throw std::runtime_error("cannot write after binary output is finished");
        if (!compressed_)
        {
            output_.write(data.data(), static_cast<std::streamsize>(data.size()));
        }
        else
        {
            stream_.next_in = reinterpret_cast<const uint8_t*>(data.data());
            stream_.avail_in = data.size();
            while (stream_.avail_in > 0)
                encode(LZMA_RUN);
        }
        if (!output_)
            throw std::runtime_error("failed while writing syndrome records");
    }

    void finish()
    {
        if (finished_)
            return;
        if (compressed_)
        {
            lzma_ret status;
            do
            {
                status = encode(LZMA_FINISH);
            } while (status != LZMA_STREAM_END);
        }
        output_.flush();
        if (!output_)
            throw std::runtime_error("failed while finishing syndrome output");
        finished_ = true;
    }

private:
    lzma_ret encode(lzma_action action)
    {
        stream_.next_out = buffer_.data();
        stream_.avail_out = buffer_.size();
        const auto status = lzma_code(&stream_, action);
        if (status != LZMA_OK && status != LZMA_STREAM_END)
            throw std::runtime_error("LZMA compression failed: " + std::to_string(status));
        const auto produced = buffer_.size() - stream_.avail_out;
        output_.write(reinterpret_cast<const char*>(buffer_.data()),
                      static_cast<std::streamsize>(produced));
        if (!output_)
            throw std::runtime_error("failed while writing compressed syndrome records");
        return status;
    }

    std::ofstream output_;
    bool compressed_;
    bool finished_{false};
    lzma_stream stream_ = LZMA_STREAM_INIT;
    std::array<uint8_t, 1 << 16> buffer_{};
};

// Merge independent mechanisms with identical detector support by odd parity.
double
xor_probability(double a, double b)
{
    return a * (1.0 - b) + b * (1.0 - a);
}

// OpenAI GPT-6: Attach a shared boundary vertex to every odd-order raw error
// term, including single-detector terms, before building the affinity graph.
AffinityGraph
build_affinity_graph(const stim::DetectorErrorModel& dem)
{
    if (dem.count_detectors() > static_cast<uint64_t>(std::numeric_limits<hg::id_type>::max()))
        throw std::runtime_error("too many detectors for Hypergraph IDs");

    AffinityGraph graph;
    for (uint64_t i = 0; i < dem.count_detectors(); ++i)
        graph.add_vertex(VertexData{});
    const auto boundary_id = static_cast<hg::id_type>(dem.count_detectors());
    if (graph.add_vertex(VertexData{}) != boundary_id)
        throw std::runtime_error("failed to assign boundary detector ID");

    std::map<std::vector<hg::id_type>, double> merged;
    dem.iter_flatten_error_instructions(
            [&](const stim::DemInstruction& instruction)
            {
                const double pr = instruction.arg_data[0];
                if (!(pr >= 0.0 && pr <= 1.0))
                    throw std::runtime_error("invalid DEM error probability");
                if (pr == 0.0)
                    return;

                std::vector<hg::id_type> support;
                for (const auto& target : instruction.target_data)
                    if (target.is_relative_detector_id())
                        support.push_back(static_cast<hg::id_type>(target.val()));

                // Repeated detector targets cancel by parity within one fault.
                std::sort(support.begin(), support.end());
                std::vector<hg::id_type> odd_support;
                for (size_t i = 0; i < support.size();)
                {
                    size_t j = i + 1;
                    while (j < support.size() && support[j] == support[i])
                        ++j;
                    if ((j - i) & 1)
                        odd_support.push_back(support[i]);
                    i = j;
                }
                if (odd_support.empty())
                    return;
                if (odd_support.size() % 2 != 0)
                    odd_support.push_back(boundary_id);
                if (odd_support.size() > MAX_ERROR_ORDER)
                    throw std::runtime_error("undecomposed DEM term exceeds MAX_ERROR_ORDER");

                auto [it, inserted] = merged.try_emplace(std::move(odd_support), pr);
                if (!inserted)
                    it->second = xor_probability(it->second, pr);
            });

    for (const auto& [support, pr] : merged)
        if (pr > 0.0)
            graph.add_edge(support, EdgeData{pr});
    return graph;
}

// OpenAI GPT-6: Reduce each detector's affinity row to its effective partner
// count, excluding self-affinity even when a decoder returns a nonzero diagonal.
// The boundary detector is included as the final row for odd physical weight.
void
write_shot(std::string& out, bool truth, bool prediction,
           size_t physical_weight, const std::vector<hg::id_type>& detector_ids,
           const pp::AffinityResult& affinity)
{
    const size_t matrix_size = detector_ids.size();
    if (physical_weight > std::numeric_limits<uint16_t>::max() ||
        matrix_size != physical_weight + physical_weight % 2 ||
        affinity.size() != matrix_size * matrix_size)
        throw std::runtime_error("unexpected affinity matrix size");

    append_u16(out, static_cast<uint16_t>(physical_weight));
    append_u8(out, static_cast<uint8_t>(physical_weight % 2));
    append_u8(out, static_cast<uint8_t>(truth != prediction));
    for (const auto id : detector_ids)
    {
        if (id > std::numeric_limits<uint16_t>::max())
            throw std::runtime_error("detector ID exceeds binary format range");
        append_u16(out, static_cast<uint16_t>(id));
    }

    for (size_t i = 0; i < matrix_size; ++i)
    {
        double row_sum{0.0};
        double row_square_sum{0.0};
        for (size_t j = 0; j < matrix_size; ++j)
        {
            if (i == j)
                continue;
            const double value = affinity[i * matrix_size + j];
            if (!(value >= 0.0 && value <= 1.0) || !std::isfinite(value))
                throw std::runtime_error("measure_affinity returned a value outside [0, 1]");
            row_sum += value;
            row_square_sum += value * value;
        }
        const double effective_partners = row_square_sum > 0.0
                ? row_sum * row_sum / row_square_sum : 0.0;
        if (!std::isfinite(effective_partners))
            throw std::runtime_error("effective partner count is not finite");
        append_f64(out, effective_partners);
    }
}

// Rank 0 owns the output file. Gather one bounded batch of binary records at a
// time so high-shot MPI runs do not keep all syndromes in memory.
void
write_batch(const std::unique_ptr<BinaryWriter>& output, const std::string& local_records,
            int rank, int world_size)
{
#if defined(ENABLE_MPI)
    if (local_records.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("one rank's binary batch exceeds MPI count range");
    const int send_count = static_cast<int>(local_records.size());
    std::vector<int> counts(rank == 0 ? world_size : 0);
    MPI_Gather(&send_count, 1, MPI_INT, rank == 0 ? counts.data() : nullptr,
               1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> offsets;
    std::vector<char> gathered;
    if (rank == 0)
    {
        offsets.resize(world_size);
        int64_t total{0};
        for (int i = 0; i < world_size; ++i)
        {
            if (total > std::numeric_limits<int>::max())
                throw std::runtime_error("combined binary batch exceeds MPI count range");
            offsets[i] = static_cast<int>(total);
            total += counts[i];
        }
        if (total > std::numeric_limits<int>::max())
            throw std::runtime_error("combined binary batch exceeds MPI count range");
        gathered.resize(static_cast<size_t>(total));
    }
    MPI_Gatherv(local_records.data(), send_count, MPI_CHAR,
                rank == 0 ? gathered.data() : nullptr,
                rank == 0 ? counts.data() : nullptr,
                rank == 0 ? offsets.data() : nullptr,
                MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank == 0)
        output->write(std::string_view(gathered.data(), gathered.size()));
#else
    (void)rank;
    (void)world_size;
    output->write(local_records);
#endif
}

} // namespace

int
run_analysis(int argc, char* argv[], int rank, int world_size)
{
    std::string output_path;
    int64_t distance{0};
    int64_t rounds{-1};
    int64_t shot_count{10000};
    int64_t seed{0};
    double p{1e-3};
    bool opposite_basis{false};

    ARGPARSE()
        .required("output-file", "Binary output path (.bin.xz or .bin)", output_path)
        .required("code-distance", "Distance of the rotated surface code", distance)
        .optional("-r", "--rounds", "Syndrome rounds (-1 = distance)", rounds, -1)
        .optional("-p", "--physical-error-rate", "SI1000 physical error rate", p, 1e-3)
        .optional("-s", "--shots", "Monte Carlo shots", shot_count, 10000)
        .optional("", "--seed", "Monte Carlo seed", seed, 0)
        .optional("", "--include-opposite-basis-detectors", "Include opposite-basis checks", opposite_basis, false)
        .parse(argc, argv);

    if (rounds == -1)
        rounds = distance;
    if (distance < 3 || distance % 2 == 0 ||
        distance > std::numeric_limits<uint32_t>::max() ||
        rounds < 1 || rounds > std::numeric_limits<uint32_t>::max() ||
        shot_count < 1 || seed < 0 || !(p > 0.0 && p < 1.0))
        throw std::invalid_argument("invalid distance, rounds, shots, seed, or physical error rate");
    if (!(output_path.ends_with(".bin.xz") || output_path.ends_with(".bin")))
        throw std::invalid_argument("output file must end in .bin.xz or .bin");
    if (distance > std::numeric_limits<uint8_t>::max() ||
        rounds > std::numeric_limits<uint8_t>::max() ||
        world_size > std::numeric_limits<uint16_t>::max())
        throw std::invalid_argument("distance, rounds, or MPI ranks exceed binary format range");

    const auto circuit = sc_si1000(static_cast<uint32_t>(distance),
                                   static_cast<uint32_t>(rounds), p, false, opposite_basis);
    // The sampled DEM keeps raw hyperedges; PyMatching receives a separate
    // decomposed DEM with the same detector and observable numbering.
    const auto affinity_dem = stim::circuit_to_dem(circuit, {.decompose_errors = false});
    const auto matching_dem = stim::circuit_to_dem(circuit, {.decompose_errors = true});
    // OpenAI GPT-6: Stim scans the DEM to count detectors, so cache the count
    // before the per-shot detector scan.
    const auto detector_count = affinity_dem.count_detectors();
    if (detector_count > std::numeric_limits<uint16_t>::max())
        throw std::invalid_argument("detector count exceeds binary format range");
    if (detector_count != matching_dem.count_detectors() ||
        affinity_dem.count_observables() != matching_dem.count_observables() ||
        affinity_dem.count_observables() != 1)
        throw std::runtime_error("affinity and matching DEM layouts differ");

    auto graph = build_affinity_graph(affinity_dem);
    dec::PyMatching decoder(matching_dem, false);
    stim::DemSampler<64> sampler(affinity_dem,
                                std::mt19937_64(static_cast<uint64_t>(seed) + rank), BATCH_SIZE);

    const uint64_t requested = static_cast<uint64_t>(shot_count);
    const uint64_t base = requested / world_size;
    const uint64_t remainder = requested % world_size;
    const uint64_t local_shots = base + (static_cast<uint64_t>(rank) < remainder);
    const uint64_t max_local_shots = base + (remainder > 0);
    const uint64_t batches = (max_local_shots + BATCH_SIZE - 1) / BATCH_SIZE;

    std::unique_ptr<BinaryWriter> output;
    if (rank == 0)
    {
        output = std::make_unique<BinaryWriter>(output_path);
        std::string header(FILE_MAGIC);
        append_u8(header, static_cast<uint8_t>(distance));
        append_u8(header, static_cast<uint8_t>(rounds));
        append_f64(header, p);
        append_u16(header, static_cast<uint16_t>(world_size));
        append_u16(header, static_cast<uint16_t>(detector_count));
        append_u8(header, static_cast<uint8_t>(1 | (opposite_basis ? 2 : 0)));
        output->write(header);
    }

    uint64_t completed{0};
    for (uint64_t batch = 0; batch < batches; ++batch)
    {
        std::string records;
        const size_t rows = static_cast<size_t>(
            std::min<uint64_t>(BATCH_SIZE, local_shots - completed));
        if (rows > 0)
        {
            sampler.resample(false);
            auto detectors = sampler.det_buffer.transposed();
            auto observables = sampler.obs_buffer.transposed();
            for (size_t row = 0; row < rows; ++row, ++completed)
            {
                std::vector<hg::id_type> detector_ids;
                for (size_t det = 0; det < detector_count; ++det)
                    if (detectors[row][det])
                        detector_ids.push_back(static_cast<hg::id_type>(det));

                // OpenAI GPT-6: The physical Hamming weight excludes the
                // boundary, while the affinity matrix includes it when odd.
                const size_t physical_weight = detector_ids.size();
                if (physical_weight % 2 != 0)
                    detector_ids.push_back(static_cast<hg::id_type>(detector_count));
                const auto decoded = decoder.decode(detectors[row], observables[row]);
                const auto affinity = pp::measure_affinity(graph, detector_ids);
                write_shot(records, observables[row][0], decoded.flipped_obs[0],
                           physical_weight, detector_ids, affinity);
            }
        }
        write_batch(output, records, rank, world_size);
    }

    uint64_t total_completed = completed;
#if defined(ENABLE_MPI)
    MPI_Reduce(&completed, &total_completed, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    if (rank == 0)
    {
        if (total_completed != requested)
            throw std::runtime_error("MPI shot counts do not sum to requested total");
        output->finish();
    }
    return 0;
}

int
main(int argc, char* argv[])
{
    int rank{0};
    int world_size{1};
#if defined(ENABLE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
    try
    {
        const int result = run_analysis(argc, argv, rank, world_size);
#if defined(ENABLE_MPI)
        MPI_Finalize();
#endif
        return result;
    }
    catch (const std::exception& error)
    {
        std::cerr << "analyze_syndrome rank " << rank << ": " << error.what() << '\n';
#if defined(ENABLE_MPI)
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
        return 1;
    }
}
