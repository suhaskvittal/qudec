/*
 *  contributor: OpenAI Codex
 *  purpose: Evaluate two-level first-decoder/Libra routing by complementary or cluster gap.
 * */

#include "circuit_generator.h"
#include "decoder/libra.h"
#include "decoder/pymatching.h"
#include "decoder/uf.h"
#include "sampler.h"
#include "stats.h"

#include <argparse/argparse.h>
#include <stim/util_top/circuit_to_dem.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace
{

enum Metric : size_t
{
    SHOTS, LOW_SHOTS, PM_ERRORS, HIGH_PM_ERRORS, LOW_PM_ERRORS,
    LIBRA_ERRORS, FRAME_CHANGES, FIXED, INTRODUCED,
    HIGH_PM_NS, LOW_PM_NS, LOW_LIBRA_NS, METRIC_COUNT
};

uint64_t
sum_across_ranks(uint64_t value)
{
#if defined(ENABLE_MPI)
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
    return value;
}

struct SavedShot
{
    dec::SyndromeType syndrome;
    bool truth;
    bool pymatching_frame;
    double gap_db;
    uint64_t pymatching_ns;
    size_t level;
};

// Stage one is compatible with the existing Monte Carlo sampler. The rare
// event driver below additionally sets current_level for probability weighting.
class GapStudy
{
public:
    dec::PyMatching pymatching;
    std::unique_ptr<dec::UnionFind> union_find;
    StatsHistogram<double> pm_high, pm_low, libra_low;
    size_t current_level{0};

    GapStudy(const stim::DetectorErrorModel& dem, double threshold, size_t levels,
             int64_t max_us, int64_t bin_us, const std::filesystem::path& output, int rank,
             const std::string& first_decoder, const std::string& route_metric, int64_t uf_weight_scale)
        :pymatching(dem, {.enable_gap_estimation=route_metric == "complementary" && first_decoder == "pymatching",
                          .enable_cluster_gap=route_metric == "cluster" && first_decoder == "pymatching"}),
        pm_high("PYMATCHING_HIGH_GAP_LATENCY_US", 0, double(max_us), max_us / bin_us),
        pm_low("PYMATCHING_LOW_GAP_LATENCY_US", 0, double(max_us), max_us / bin_us),
        libra_low("LIBRA_ON_PYMATCHING_LOW_GAP_LATENCY_US", 0, double(max_us), max_us / bin_us),
        threshold_(threshold), counters_(levels * METRIC_COUNT, 0), output_(output), rank_(rank),
        bits_(output / ("low_gap_syndromes_rank_" + std::to_string(rank) + ".b8"), std::ios::binary)
    {
        if (first_decoder == "uf")
            union_find = std::make_unique<dec::UnionFind>(dem, dec::UnionFind::config_type{
                .weight_scale=static_cast<uint64_t>(uf_weight_scale), .enable_cluster_gap=true});
        first_decoder_ = first_decoder;
        route_metric_ = route_metric;
        pm_high.name = "FIRST_DECODER_HIGH_GAP_LATENCY_US";
        pm_low.name = "FIRST_DECODER_LOW_GAP_LATENCY_US";
        libra_low.name = "LIBRA_ON_FIRST_DECODER_LOW_GAP_LATENCY_US";
        if (!bits_)
            throw std::runtime_error("Cannot create low-gap syndrome file");
    }

    dec::Result decode(dec::SyndromeRef syndrome, dec::ObsRef truth)
    {
        auto result = union_find ? union_find->decode(syndrome, truth) : pymatching.decode(syndrome, truth);
        // The sign of the complementary gap reveals the true frame, so route on
        // its magnitude. Cluster gap is already nonnegative and truth independent.
        if (route_metric_ == "cluster" && (!result.has_cluster_gap || !std::isfinite(result.cluster_gap_db) ||
                                            result.cluster_gap_db < 0))
            throw std::runtime_error("First decoder did not return a finite nonnegative cluster gap");
        const double gap_db = route_metric_ == "cluster" ? result.cluster_gap_db
                                                          : std::abs(result.matching_data.gap);
        const bool low = gap_db < threshold_;
        const bool failed = result.flipped_obs[0] != truth[0];
        add(SHOTS, 1);
        add(PM_ERRORS, failed);
        add(low ? LOW_PM_ERRORS : HIGH_PM_ERRORS, failed);
        add(low ? LOW_PM_NS : HIGH_PM_NS, result.decode_latency_ns);
        (low ? pm_low : pm_high).add(double(result.decode_latency_ns) / 1000.0);
        if (low)
        {
            add(LOW_SHOTS, 1);
            SavedShot saved{dec::SyndromeType(pymatching.num_detectors), bool(truth[0]),
                            bool(result.flipped_obs[0]), gap_db, result.decode_latency_ns, current_level};
            saved.syndrome = syndrome;
            saved_.push_back(std::move(saved));
            // Store only detector bits, in Stim b8 order, with zero padding.
            for (size_t i = 0; i < (pymatching.num_detectors + 7) / 8; i++)
            {
                uint8_t byte = uint8_t(syndrome.u64[i / 8] >> (8 * (i % 8)));
                if (i == pymatching.num_detectors / 8 && pymatching.num_detectors % 8)
                    byte &= uint8_t((1u << (pymatching.num_detectors % 8)) - 1);
                bits_.put(char(byte));
            }
            if (!bits_)
                throw std::runtime_error("Cannot record low-gap syndrome");
        }
        return result;
    }

    void decode_saved(const stim::DetectorErrorModel& dem, dec::Libra::Config config)
    {
        bits_.close();
        if (!bits_)
            throw std::runtime_error("Cannot finish low-gap syndrome file");
        std::ofstream records(output_ / ("low_gap_records_rank_" + std::to_string(rank_) + ".csv"));
        records << "index,error_count_level,gap_db,first_decoder_latency_ns,libra_latency_ns,"
                   "first_decoder_frame,libra_frame,true_frame\n";
        records << std::setprecision(17);
        if (!saved_.empty())
        {
            dec::Libra libra(dem, config);
            libra.prepare_ensemble();
            dec::ObsType truth(1);
            for (size_t i = 0; i < saved_.size(); i++)
            {
                auto& saved = saved_[i];
                current_level = saved.level;
                truth[0] = saved.truth;
                auto result = libra.decode(saved.syndrome, truth);
                const bool frame = result.flipped_obs[0];
                const bool pm_failed = saved.pymatching_frame != saved.truth;
                const bool libra_failed = frame != saved.truth;
                add(LIBRA_ERRORS, libra_failed);
                add(FRAME_CHANGES, frame != saved.pymatching_frame);
                add(FIXED, pm_failed && !libra_failed);
                add(INTRODUCED, !pm_failed && libra_failed);
                add(LOW_LIBRA_NS, result.decode_latency_ns);
                libra_low.add(double(result.decode_latency_ns) / 1000.0);
                records << i << ',' << saved.level << ',' << saved.gap_db << ',' << saved.pymatching_ns
                        << ',' << result.decode_latency_ns << ',' << saved.pymatching_frame
                        << ',' << frame << ',' << saved.truth << '\n';
            }
        }
        records.close();
        if (!records)
            throw std::runtime_error("Cannot write low-gap latency records");
    }

    std::vector<uint64_t> global_counters() const
    {
        auto counters = counters_;
#if defined(ENABLE_MPI)
        MPI_Allreduce(MPI_IN_PLACE, counters.data(), int(counters.size()), MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
        return counters;
    }

private:
    double threshold_;
    std::vector<uint64_t> counters_;
    std::vector<SavedShot> saved_;
    std::filesystem::path output_;
    int rank_;
    std::string first_decoder_;
    std::string route_metric_;
    std::ofstream bits_;

    void add(Metric metric, uint64_t value) { counters_[current_level * METRIC_COUNT + metric] += value; }
};

// Use the corrected heterogeneous conditional sampler. Keep count levels with
// each retained shot so conditional latency averages can be weighted by P(K).
double
sample_rare(GapStudy& study, const ConditionalErrorSampler& sampler, const ExperimentConfig& conf,
            int rank, int workers, size_t& first_omitted, double& omitted_bound)
{
    const auto& probabilities = sampler.count_probabilities();
    const size_t limit = probabilities.size();
    omitted_bound = sampler.tail_probability_bound(limit);
    if (omitted_bound > conf.rare_event.tail_tolerance)
        throw std::invalid_argument("Rare-event omitted-count bound exceeds tolerance; increase --rare-max-level");
    std::vector<double> remaining_bounds(limit + 1, omitted_bound);
    long double remaining = omitted_bound;
    for (size_t k = limit; k-- > 0; )
    {
        remaining += probabilities[k];
        remaining_bounds[k] = std::nextafter(double(std::min(remaining, 1.0L)), 1.0);
    }
    first_omitted = limit;
    Rng rng(uint64_t(conf.seed) + rank);
    double ler = 0;
    for (size_t k = 0; k < limit; k++)
    {
        if (conf.rare_event.relative_tail_tolerance > 0 && ler > 0 &&
            remaining_bounds[k] < conf.rare_event.relative_tail_tolerance * ler)
        {
            first_omitted = k;
            omitted_bound = remaining_bounds[k];
            break;
        }
        if (!(probabilities[k] > 0))
            continue;
        study.current_level = k;
        uint64_t shots = 0, errors = 0;
        while (shots < uint64_t(conf.rare_event.samples_per_level) &&
               errors < uint64_t(conf.rare_event.max_errors_per_level))
        {
            const uint64_t batch = std::min<uint64_t>(conf.rare_event.samples_per_level - shots, 256ULL * workers);
            const size_t local = batch / workers + (uint64_t(rank) < batch % workers);
            auto [syndromes, truths] = sampler.sample_syndromes(k, local, rng);
            uint64_t new_errors = 0;
            for (size_t i = 0; i < local; i++)
            {
                auto result = study.decode(syndromes[i], truths[i]);
                new_errors += result.flipped_obs[0] != truths[i][0];
            }
            errors += sum_across_ranks(new_errors);
            shots += batch;
        }
        ler += probabilities[k] * double(errors) / shots;
        if (rank == 0 && conf.print_progress)
            std::cout << "[K=" << k << "] samples=" << shots << " first_decoder_errors=" << errors
                      << " probability=" << probabilities[k] << '\n';
    }
    return ler;
}

std::filesystem::path
prepare_output(std::string directory, int64_t distance, int rank)
{
    if (rank == 0)
    {
        if (directory.empty())
        {
            const auto timestamp = std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::system_clock::now().time_since_epoch()).count();
            directory = "out/qudec_eval_gap_d" + std::to_string(distance) + "_" + std::to_string(timestamp);
        }
        // Reject existing run directories so recordings cannot be overwritten.
        if (!std::filesystem::create_directories(directory))
            throw std::runtime_error("Output run directory already exists: " + directory);
    }
#if defined(ENABLE_MPI)
    uint64_t length = directory.size();
    MPI_Bcast(&length, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (length > std::numeric_limits<int>::max())
        throw std::invalid_argument("Output path is too long");
    directory.resize(length);
    MPI_Bcast(directory.data(), int(length), MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
    return directory;
}

void
print_histogram(std::ostream& out, const StatsHistogram<double>& histogram)
{
    out << histogram.to_string_full();
}

int
run(int argc, char** argv, int rank, int workers)
{
    int64_t distance, rounds, ensemble_size = 60, ensemble_seed = 0x4c49425241ULL;
    int64_t latency_max_us = 20000, latency_bin_us = 1;
    double p, threshold_db;
    bool memory_x = false;
    bool print_histograms = false;
    std::string output_directory;
    std::string first_decoder = "pymatching", route_metric = "complementary";
    int64_t uf_weight_scale = 1024;
    ExperimentConfig conf;
    dec::Libra::Config libra_config;
    ARGPARSE()
        .required("distance", "Surface-code distance", distance)
        .optional("-r", "--rounds", "Rounds (-1 = distance)", rounds, -1)
        .optional("-p", "--physical-error-rate", "SI1000 physical error rate", p, 0.001)
        .optional("", "--memory-x", "Use X memory (default: Z)", memory_x, false)
        .optional("", "--first-decoder", "First-stage decoder: pymatching or uf", first_decoder, "pymatching")
        .optional("", "--gap-metric", "Routing metric: complementary or cluster", route_metric, "complementary")
        .optional("", "--uf-weight-scale", "UF quantization units per natural-log weight", uf_weight_scale, 1024)
        .optional("-m", "--method", "monte_carlo or rare_event", conf.method, "monte_carlo")
        .optional("", "--gap-threshold-db", "Forward to Libra when selected first-decoder gap is below this threshold", threshold_db, 20.0)
        .optional("", "--seed", "Sampling seed; ranks use distinct streams", conf.seed, 0)
        .optional("", "--mc-max-samples", "Global Monte Carlo shot budget", conf.monte_carlo.max_samples, 1000000)
        .optional("", "--mc-stop-at-errors", "Stop at this many first-stage PyMatching errors", conf.monte_carlo.stop_at_error_count, std::numeric_limits<int64_t>::max())
        .optional("", "--rare-samples-per-level", "Global samples per active-error count", conf.rare_event.samples_per_level, 10000)
        .optional("", "--rare-max-errors-per-level", "Stop a count level at this many PyMatching errors", conf.rare_event.max_errors_per_level, 25)
        .optional("", "--rare-max-level", "Exclusive active-error count limit", conf.rare_event.max_level, 128)
        .optional("", "--rare-tail-tolerance", "Maximum omitted-count bound beyond configured limit", conf.rare_event.tail_tolerance, 1e-12)
        .optional("", "--rare-relative-tail-tolerance", "Remaining mass / running PyMatching LER stop threshold; 0 disables", conf.rare_event.relative_tail_tolerance, 0.1)
        .optional("", "--libra-ensemble-size", "Libra ensemble size", ensemble_size, 60)
        .optional("", "--libra-seed", "Shared Libra ensemble seed", ensemble_seed, ensemble_seed)
        .optional("", "--libra-gap-threshold-db", "Libra's own ensemble trigger gap", libra_config.gap_threshold_db, 20.0)
        .optional("", "--libra-validate", "Validate Libra's corrections", libra_config.validate, false)
        .optional("", "--latency-bin-us", "Latency histogram bin width (integer microseconds)", latency_bin_us, 1)
        .optional("", "--latency-max-us", "Latency histogram exclusive upper limit", latency_max_us, 20000)
        .optional("", "--print-histograms", "Print latency histogram bins to console (always saved in summary.txt)", print_histograms, false)
        .optional("", "--output-dir", "New run directory (default: unique directory under out/)", output_directory, "")
        .optional("", "--print-progress", "Print first-stage sampling progress", conf.print_progress, false)
        .parse(argc, argv);

    if (rounds == -1) rounds = distance;
    if (distance < 2 || distance > std::numeric_limits<uint32_t>::max() || rounds < 1 ||
        rounds > std::numeric_limits<uint32_t>::max() || !std::isfinite(p) || p <= 0 || p >= 0.5 ||
        !std::isfinite(threshold_db) || threshold_db < 0 || conf.seed < 0 ||
        conf.seed > std::numeric_limits<int64_t>::max() - (workers - 1) || ensemble_seed < 0 || ensemble_size < 0 ||
        !std::isfinite(libra_config.gap_threshold_db) || libra_config.gap_threshold_db < 0 ||
        conf.monte_carlo.max_samples <= 0 || conf.monte_carlo.stop_at_error_count <= 0 ||
        conf.rare_event.samples_per_level <= 0 || conf.rare_event.max_errors_per_level <= 0 ||
        conf.rare_event.max_level <= 0 || conf.rare_event.max_level > std::numeric_limits<int>::max() / METRIC_COUNT ||
        !std::isfinite(conf.rare_event.tail_tolerance) || conf.rare_event.tail_tolerance < 0 || conf.rare_event.tail_tolerance >= 1 ||
        !std::isfinite(conf.rare_event.relative_tail_tolerance) || conf.rare_event.relative_tail_tolerance < 0 || conf.rare_event.relative_tail_tolerance > 1 ||
        latency_bin_us <= 0 || latency_max_us <= 0 || latency_max_us % latency_bin_us != 0 ||
        latency_max_us / latency_bin_us > std::numeric_limits<int>::max() - 2)
        throw std::invalid_argument("Invalid code, gap, sampling, or latency histogram parameters");
    if (conf.method != "monte_carlo" && conf.method != "rare_event")
        throw std::invalid_argument("Unknown sampling method: " + conf.method);
    if (first_decoder != "pymatching" && first_decoder != "uf")
        throw std::invalid_argument("First decoder must be pymatching or uf");
    if (route_metric != "complementary" && route_metric != "cluster")
        throw std::invalid_argument("Gap metric must be complementary or cluster");
    if (route_metric == "complementary" && first_decoder != "pymatching")
        throw std::invalid_argument("Complementary-gap routing requires --first-decoder pymatching");
    if (uf_weight_scale <= 0)
        throw std::invalid_argument("UF weight scale must be positive");
    libra_config.ensemble_size = size_t(ensemble_size);
    libra_config.seed = uint64_t(ensemble_seed);
    auto circuit = sc_si1000(distance, rounds, p, memory_x, true);
    auto dem = stim::circuit_to_dem(circuit, {.decompose_errors=true});
    const auto output = prepare_output(output_directory, distance, rank);
    std::unique_ptr<ConditionalErrorSampler> conditional;
    std::vector<double> probabilities{1.0};
    if (conf.method == "rare_event")
    {
        conditional = std::make_unique<ConditionalErrorSampler>(dem, std::min<size_t>(conf.rare_event.max_level, dem.count_errors() + 1));
        probabilities = conditional->count_probabilities();
    }
    GapStudy study(dem, threshold_db, probabilities.size(), latency_max_us, latency_bin_us, output, rank,
                   first_decoder, route_metric, uf_weight_scale);
    if (rank == 0)
        std::cout << "Two-level gap study: distance=" << distance << " rounds=" << rounds << " p=" << p
                  << " opposite_basis_detectors=true first_decoder=" << first_decoder << " gap_metric=" << route_metric
                  << " correlated_pymatching=false gap_threshold_db=" << threshold_db
                  << " sampler=" << conf.method << " ranks=" << workers << " output=" << output << '\n';

    size_t first_omitted = probabilities.size();
    double omitted_bound = 0;
    const double pm_ler = conditional ? sample_rare(study, *conditional, conf, rank, workers, first_omitted, omitted_bound)
                                     : estimate_logical_error_rate(dem, study, conf);
    const auto stage_one = study.global_counters();
    study.pm_high.mpi_accumulate();
    study.pm_low.mpi_accumulate();
    if (rank == 0)
    {
        uint64_t shots = 0, low = 0;
        for (size_t k = 0; k < probabilities.size(); k++)
        {
            shots += stage_one[k * METRIC_COUNT + SHOTS];
            low += stage_one[k * METRIC_COUNT + LOW_SHOTS];
        }
        std::cout << "Stage 1 complete: shots=" << shots << " low_gap_shots=" << low
                  << " first_decoder=" << first_decoder << " logical_error_rate=" << pm_ler << '\n';
        if (print_histograms)
        {
            print_histogram(std::cout, study.pm_high);
            print_histogram(std::cout, study.pm_low);
        }
        std::cout << "Stage 2: decoding retained low-gap syndromes with Libra\n";
    }
    study.decode_saved(dem, libra_config);
    const auto counts = study.global_counters();
    study.libra_low.mpi_accumulate();
    if (rank != 0)
        return 0;

    std::vector<uint64_t> totals(METRIC_COUNT, 0);
    double retained_probability = 0, low_probability = 0, weighted_ler = 0;
    double weighted_latency = 0, sequential_latency = 0;
    std::ofstream levels(output / "error_count_levels.csv");
    levels << "level,probability,shots,low_shots,first_decoder_errors,two_level_errors,high_first_decoder_ns,low_first_decoder_ns,low_libra_ns\n";
    levels << std::setprecision(17);
    for (size_t k = 0; k < probabilities.size(); k++)
    {
        const auto* row = counts.data() + k * METRIC_COUNT;
        for (size_t m = 0; m < METRIC_COUNT; m++) totals[m] += row[m];
        if (!row[SHOTS]) continue;
        const double weight = probabilities[k] / double(row[SHOTS]);
        retained_probability += probabilities[k];
        low_probability += weight * row[LOW_SHOTS];
        weighted_ler += weight * (row[HIGH_PM_ERRORS] + row[LIBRA_ERRORS]);
        weighted_latency += weight * (double(row[HIGH_PM_NS]) + double(row[LOW_LIBRA_NS])) / 1000.0;
        sequential_latency += weight * (double(row[HIGH_PM_NS]) + double(row[LOW_PM_NS]) + double(row[LOW_LIBRA_NS])) / 1000.0;
        levels << k << ',' << probabilities[k] << ',' << row[SHOTS] << ',' << row[LOW_SHOTS] << ',' << row[PM_ERRORS]
               << ',' << row[HIGH_PM_ERRORS] + row[LIBRA_ERRORS] << ',' << row[HIGH_PM_NS] << ',' << row[LOW_PM_NS]
               << ',' << row[LOW_LIBRA_NS] << '\n';
    }
    levels.close();
    if (!levels) throw std::runtime_error("Cannot write count-level statistics");
    if (!totals[SHOTS]) throw std::runtime_error("No sampled shots");
    std::ostringstream report;
    report << std::setprecision(12)
           << "shots=" << totals[SHOTS] << " high_gap_shots=" << totals[SHOTS] - totals[LOW_SHOTS]
           << " low_gap_shots=" << totals[LOW_SHOTS] << '\n'
           << "first_decoder_logical_error_rate=" << pm_ler << '\n'
           << "first_decoder_errors=" << totals[PM_ERRORS] << " high_gap_errors=" << totals[HIGH_PM_ERRORS]
           << " low_gap_errors=" << totals[LOW_PM_ERRORS] << " libra_low_gap_errors=" << totals[LIBRA_ERRORS] << '\n'
           << "libra_frame_changes=" << totals[FRAME_CHANGES] << " errors_fixed=" << totals[FIXED]
           << " errors_introduced=" << totals[INTRODUCED] << '\n'
           << "two_level_logical_error_rate=" << weighted_ler << " retained_probability=" << retained_probability
           << " low_gap_probability=" << low_probability << '\n'
           << "two_level_sample_mean_latency_us=" << (double(totals[HIGH_PM_NS]) + double(totals[LOW_LIBRA_NS])) / totals[SHOTS] / 1000.0 << '\n'
           << "two_level_probability_weighted_mean_latency_us=" << weighted_latency << '\n'
           << "sequential_probability_weighted_mean_latency_us=" << sequential_latency << '\n'
           << "Latency scope: full history; requested two-level mean uses first-decoder high-gap + Libra low-gap.\n"
           << "Sequential mean additionally includes first-decoder gap estimation on low-gap shots.\n"
           << "Histograms describe sampled shots; rare-event physical means use P(K) weighting without renormalizing.\n";
    if (conditional)
        report << "omitted_probability_bound=" << omitted_bound << " first_omitted_level=" << first_omitted
               << " (bounds missing logical-error contribution, not missing latency)\n";
    std::cout << report.str();
    if (print_histograms)
        print_histogram(std::cout, study.libra_low);
    // Keep all distributions in the saved report without flooding the console.
    print_histogram(report, study.pm_high);
    print_histogram(report, study.pm_low);
    print_histogram(report, study.libra_low);
    std::ofstream summary(output / "summary.txt");
    summary << "distance=" << distance << " rounds=" << rounds << " p=" << std::setprecision(17) << p
            << " seed=" << conf.seed << " first_decoder=" << first_decoder << " gap_metric=" << route_metric
            << " uf_weight_scale=" << uf_weight_scale << " gap_threshold_db=" << threshold_db << " method=" << conf.method
            << " memory_basis=" << (memory_x ? 'X' : 'Z') << " opposite_basis_detectors=true"
            << " num_detectors=" << dem.count_detectors() << " num_observables=" << dem.count_observables()
            << " ensemble_size=" << ensemble_size << " ensemble_seed=" << ensemble_seed
            << " libra_gap_threshold_db=" << libra_config.gap_threshold_db << '\n' << report.str();
    summary.close();
    if (!summary) throw std::runtime_error("Cannot write study summary");
    return 0;
}

} // anon

int
main(int argc, char** argv)
{
    int rank = 0, workers = 1;
#if defined(ENABLE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &workers);
#endif
    int status = 0;
    try { status = run(argc, argv, rank, workers); }
    catch (const std::exception& error)
    {
        if (rank == 0) std::cerr << "qudec_eval_gap: " << error.what() << '\n';
        status = 1;
#if defined(ENABLE_MPI)
        MPI_Abort(MPI_COMM_WORLD, status);
#endif
    }
#if defined(ENABLE_MPI)
    MPI_Finalize();
#endif
    return status;
}
