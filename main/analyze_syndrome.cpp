/*
 * author: OpenAI GPT-6
 * date: 24 September 2026
 * purpose: Record SI1000 surface-code syndromes, PyMatching outcomes, and
 *          unnormalized detector-affinity matrices as JSON Lines.
 */

#include "circuit_generator.h"
#include "decoder/surface_code.h"
#include "hypergraph.h"
#include "preprocessor/affinity.h"

#include <stim.h>
#include <stim/util_top/circuit_to_dem.h>

#include <argparse/argparse.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
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

// Preserve the raw row-major array returned by measure_affinity(). Detector
// IDs identify its row and column order; no normalization is applied here.
void
write_shot(std::ostream& out, uint64_t shot_index, bool truth, bool prediction,
           size_t physical_weight, const std::vector<hg::id_type>& detector_ids,
           const pp::AffinityResult& affinity)
{
    const size_t matrix_size = detector_ids.size();
    if (matrix_size != physical_weight + physical_weight % 2 ||
        affinity.size() != matrix_size * matrix_size)
        throw std::runtime_error("unexpected affinity matrix size");

    out << "{\"type\":\"shot\",\"shot_index\":" << shot_index
        << ",\"hamming_weight\":" << physical_weight
        << ",\"boundary_added\":" << (physical_weight % 2 ? "true" : "false")
        << ",\"logical_error\":" << (truth != prediction ? "true" : "false")
        << ",\"true_observable\":" << int(truth)
        << ",\"pymatching_prediction\":" << int(prediction)
        << ",\"detector_ids\":[";
    for (size_t i = 0; i < matrix_size; ++i)
    {
        if (i > 0) out << ',';
        out << detector_ids[i];
    }
    out << "],\"affinity\":[";
    for (size_t i = 0; i < affinity.size(); ++i)
    {
        const double value = affinity[i];
        if (!(value >= 0.0 && value <= 1.0) || !std::isfinite(value))
            throw std::runtime_error("measure_affinity returned a value outside [0, 1]");
        if (i > 0) out << ',';
        out << value;
    }
    out << "]}\n";
}

// Rank 0 owns the output file. Gather one bounded batch of JSONL records at a
// time so high-shot MPI runs do not keep all syndrome matrices in memory.
void
write_batch(std::ofstream& output, const std::string& local_records,
            int rank, int world_size)
{
#if defined(ENABLE_MPI)
    if (local_records.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("one rank's JSONL batch exceeds MPI count range");
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
                throw std::runtime_error("combined JSONL batch exceeds MPI count range");
            offsets[i] = static_cast<int>(total);
            total += counts[i];
        }
        if (total > std::numeric_limits<int>::max())
            throw std::runtime_error("combined JSONL batch exceeds MPI count range");
        gathered.resize(static_cast<size_t>(total));
    }
    MPI_Gatherv(local_records.data(), send_count, MPI_CHAR,
                rank == 0 ? gathered.data() : nullptr,
                rank == 0 ? counts.data() : nullptr,
                rank == 0 ? offsets.data() : nullptr,
                MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank == 0)
        output.write(gathered.data(), static_cast<std::streamsize>(gathered.size()));
#else
    (void)rank;
    (void)world_size;
    output << local_records;
#endif
    if (rank == 0)
    {
        output.flush();
        if (!output)
            throw std::runtime_error("failed while writing syndrome records");
    }
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
        .required("output-file", "JSON Lines output path", output_path)
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

    const auto circuit = sc_si1000(static_cast<uint32_t>(distance),
                                   static_cast<uint32_t>(rounds), p, false, opposite_basis);
    // The sampled DEM keeps raw hyperedges; PyMatching receives a separate
    // decomposed DEM with the same detector and observable numbering.
    const auto affinity_dem = stim::circuit_to_dem(circuit, {.decompose_errors = false});
    const auto matching_dem = stim::circuit_to_dem(circuit, {.decompose_errors = true});
    if (affinity_dem.count_detectors() != matching_dem.count_detectors() ||
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
    const uint64_t rank_start = static_cast<uint64_t>(rank) * base +
                                std::min<uint64_t>(static_cast<uint64_t>(rank), remainder);
    const uint64_t max_local_shots = base + (remainder > 0);
    const uint64_t batches = (max_local_shots + BATCH_SIZE - 1) / BATCH_SIZE;

    std::ofstream output;
    if (rank == 0)
    {
        output.open(output_path, std::ios::out | std::ios::trunc);
        if (!output)
            throw std::runtime_error("failed to open output file: " + output_path);
        output << std::setprecision(17)
               << "{\"type\":\"metadata\",\"schema\":2"
               << ",\"distance\":" << distance
               << ",\"rounds\":" << rounds
               << ",\"physical_error_rate\":" << p
               << ",\"memory_basis\":\"Z\""
               << ",\"include_opposite_basis_detectors\":" << (opposite_basis ? "true" : "false")
               << ",\"requested_shots\":" << requested
               << ",\"seed\":" << seed
               << ",\"mpi_ranks\":" << world_size
               << ",\"rank_seed\":\"seed_plus_rank\""
               << ",\"detector_count\":" << affinity_dem.count_detectors()
               << ",\"boundary_detector_id\":" << affinity_dem.count_detectors()
               << ",\"affinity_graph_edges\":" << graph.edge_count()
               << ",\"affinity_dem_decompose_errors\":false"
               << ",\"pymatching_dem_decompose_errors\":true"
               << ",\"affinity_layout\":\"row-major square matrix in detector_ids order\""
               << ",\"affinity_normalized\":false"
               << ",\"boundary_rule\":\"append boundary to odd-order errors and odd-weight syndromes\"}\n";
        if (!output)
            throw std::runtime_error("failed while writing metadata");
    }

    uint64_t completed{0};
    for (uint64_t batch = 0; batch < batches; ++batch)
    {
        std::ostringstream records;
        records << std::setprecision(17);
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
                for (size_t det = 0; det < affinity_dem.count_detectors(); ++det)
                    if (detectors[row][det])
                        detector_ids.push_back(static_cast<hg::id_type>(det));

                // OpenAI GPT-6: The physical Hamming weight excludes the
                // boundary, while the affinity matrix includes it when odd.
                const size_t physical_weight = detector_ids.size();
                if (physical_weight % 2 != 0)
                    detector_ids.push_back(static_cast<hg::id_type>(affinity_dem.count_detectors()));
                const auto decoded = decoder.decode(detectors[row], observables[row]);
                const auto affinity = pp::measure_affinity(graph, detector_ids);
                write_shot(records, rank_start + completed,
                           observables[row][0], decoded.flipped_obs[0],
                           physical_weight, detector_ids, affinity);
            }
        }
        write_batch(output, records.str(), rank, world_size);
    }

    uint64_t total_completed = completed;
#if defined(ENABLE_MPI)
    MPI_Reduce(&completed, &total_completed, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    if (rank == 0)
    {
        if (total_completed != requested)
            throw std::runtime_error("MPI shot counts do not sum to requested total");
        output << "{\"type\":\"summary\",\"completed_shots\":" << total_completed << "}\n";
        output.flush();
        if (!output)
            throw std::runtime_error("failed while writing syndrome summary");
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
