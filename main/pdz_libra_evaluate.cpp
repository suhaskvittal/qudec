/*
 *  contributor: OpenAI Codex
 *  purpose: Estimate the logical error rate of PyMatching -> PDZ -> Libra.
 * */

#include "circuit_generator.h"
#include "decoder/pdz_libra.h"
#include "sampler.h"

#include <stim/util_top/circuit_to_dem.h>

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{

struct Options
{
    uint32_t distance{5};
    uint32_t rounds{0};
    double physical_error_rate{1e-3};
    bool memory_x{false};
    ExperimentConfig experiment;
    pp::PDZ::config_type pdz;
    dec::Libra::Config libra;
};

void
print_usage(std::ostream& out)
{
    out << "usage: pdz_libra_evaluate DISTANCE [options]\n"
        << "  --rounds N                    default: distance\n"
        << "  --p P                         default: 0.001\n"
        << "  --memory-basis X|Z            default: Z\n"
        << "  --seed N                      default: 0\n"
        << "  --method monte_carlo|rare_event  default: monte_carlo\n"
        << "  --mc-max-samples N             default: 1000000\n"
        << "  --mc-stop-at-errors N          default: 25 (checked per batch)\n"
        << "  --aff-keep-threshold N         default: 2\n"
        << "  --libra-ensemble-size N        default: 60\n"
        << "  --rare-samples-per-level N     default: 10000\n"
        << "  --rare-max-errors-per-level N  default: 25\n"
        << "  --rare-max-level N             default: 128 (exclusive)\n"
        << "  --rare-relative-tail-tolerance P  default: 0.1 (0 disables early stop)\n"
        << "  --rare-tail-tolerance P        default: 1e-12\n"
        << "Opposing-basis detectors are enabled; PDZ uses same-basis scoring.\n";
}

uint32_t
positive_u32(const std::string& value)
{
    const auto number = std::stoll(value);
    if (number <= 0 || number > std::numeric_limits<uint32_t>::max())
        throw std::invalid_argument("expected a positive 32-bit integer: " + value);
    return static_cast<uint32_t>(number);
}

Options
parse_options(int argc, char** argv)
{
    Options out;
    if (argc < 2)
        throw std::invalid_argument("distance is required; use --help for usage");
    out.distance = positive_u32(argv[1]);
    for (int i = 2; i < argc; i++)
    {
        const std::string option = argv[i];
        if (++i >= argc)
            throw std::invalid_argument("missing value for " + option);
        const std::string value = argv[i];
        if (option == "--rounds") out.rounds = positive_u32(value);
        else if (option == "--p") out.physical_error_rate = std::stod(value);
        else if (option == "--seed") out.experiment.seed = std::stoll(value);
        else if (option == "--method") out.experiment.method = value;
        else if (option == "--mc-max-samples") out.experiment.monte_carlo.max_samples = positive_u32(value);
        else if (option == "--mc-stop-at-errors") out.experiment.monte_carlo.stop_at_error_count = positive_u32(value);
        else if (option == "--libra-ensemble-size") out.libra.ensemble_size = positive_u32(value);
        else if (option == "--rare-samples-per-level") out.experiment.rare_event.samples_per_level = positive_u32(value);
        else if (option == "--rare-max-errors-per-level") out.experiment.rare_event.max_errors_per_level = positive_u32(value);
        else if (option == "--rare-max-level") out.experiment.rare_event.max_level = positive_u32(value);
        else if (option == "--rare-relative-tail-tolerance") out.experiment.rare_event.relative_tail_tolerance = std::stod(value);
        else if (option == "--rare-tail-tolerance") out.experiment.rare_event.tail_tolerance = std::stod(value);
        else if (option == "--memory-basis")
        {
            if (value != "X" && value != "Z")
                throw std::invalid_argument("memory basis must be X or Z");
            out.memory_x = value == "X";
        }
        else if (option == "--aff-keep-threshold")
        {
            const auto threshold = std::stoll(value);
            if (threshold < 0 || threshold > pp::AffinityTable::MASK)
                throw std::invalid_argument("affinity threshold is outside the table range");
            out.pdz.aff_keep_threshold = static_cast<pp::PDZ::value_type>(threshold);
        }
        else throw std::invalid_argument("unknown option: " + option);
    }
    if (!std::isfinite(out.physical_error_rate) || out.physical_error_rate <= 0 || out.physical_error_rate >= 0.5)
        throw std::invalid_argument("physical error rate must be between 0 and 0.5");
    if (out.experiment.method != "monte_carlo" && out.experiment.method != "rare_event")
        throw std::invalid_argument("method must be monte_carlo or rare_event");
    if (out.experiment.seed < 0)
        throw std::invalid_argument("seed must be nonnegative");
    if (!out.rounds) out.rounds = out.distance;
    out.experiment.print_progress = true;
    return out;
}

int
run(int argc, char** argv, int rank, int workers)
{
    if (argc == 2 && std::string(argv[1]) == "--help")
    {
        if (!rank) print_usage(std::cout);
        return 0;
    }
    auto options = parse_options(argc, argv);
    options.libra.seed = (static_cast<uint64_t>(options.experiment.seed) + rank) ^ 0x4c49425241ULL;
    const auto circuit = sc_si1000(options.distance, options.rounds,
            options.physical_error_rate, options.memory_x, true);
    const auto dem = stim::circuit_to_dem(circuit, {.decompose_errors=true});
    if (!rank)
        std::cout << "decoder=uncorrelated_pymatching+pdz+libra distance=" << options.distance
                  << " rounds=" << options.rounds << " p=" << options.physical_error_rate
                  << " memory_basis=" << (options.memory_x ? 'X' : 'Z')
                  << " opposing_basis_detectors=true ignore_opposing_basis_detectors=true"
                  << " aff_keep_threshold=" << int(options.pdz.aff_keep_threshold)
                  << " ensemble_size=" << options.libra.ensemble_size
                  << " seed=" << options.experiment.seed << " workers=" << workers
                  << " method=" << options.experiment.method << '\n';

    dec::PDZLibra decoder(dem, options.pdz, options.libra);
    const double ler = estimate_logical_error_rate(dem, decoder, options.experiment);
    decoder.mpi_accumulate();
    if (!rank)
    {
        std::cout << std::setprecision(12) << "logical_error_rate=" << ler << '\n';
        decoder.print_stats(std::cout);
    }
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
        if (!rank) std::cerr << "pdz_libra_evaluate: " << error.what() << '\n';
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
