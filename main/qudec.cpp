/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#include "circuit_generator.h"
#include "decoder/surface_code.h"
#include "sampler.h"

#include <stim/gen/gen_surface_code.h>
#include <stim/util_top/circuit_to_dem.h>

#include <argparse/argparse.h>

#include <iostream>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int
main(int argc, char* argv[])
{
#if defined(ENABLE_MPI)
    MPI_Init(&argc, &argv);
#endif
    int world_rank{0};
#if defined(ENABLE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif

    /*
     * Simulation configuration:
     * */
    std::string decoder_name;
    int64_t d,
            r;
    double p;
    ExperimentConfig conf;

    /*
     * PyMatching parameters:
     * */
    bool enable_gap{false};

    ARGPARSE()
        .required("decoder-name", "Decoder to run", decoder_name)
        .required("code-distance", "Code distance of surface code", d)
        .optional("-p", "--physical-error-rate", "Physical error rate", p, 1e-3)
        .optional("-r", "--rounds", "Number of rounds (-1 = same as code distance)", r, -1)
        .optional("-v", "--verbose", "Verbosity level", conf.verbosity, 0)
        .optional("-m", "--method", "Sampler method: monte_carlo or rare_event", conf.method, "monte_carlo")
        .optional("-pp", "--print-progress", "Print simulation progress", conf.print_progress, false)
        .optional("-g", "--gap", "Enable complementary gap estimation (pymatching only)", enable_gap, false)

        .optional("", "--mc-max-samples", "Monte-carlo: max shots to sample", conf.monte_carlo.max_samples, 1000000)
        .optional("", "--mc-stop-at-errors", "Monte-carlo: stop after this many logical errors", conf.monte_carlo.stop_at_error_count, 25)

        .optional("", "--rare-samples-per-level", "Rare-event: samples per error level", conf.rare_event.samples_per_level, 10000)
        .optional("", "--rare-max-errors-per-level", "Rare-event: max errors per error level", conf.rare_event.max_errors_per_level, 25)

        .parse(argc, argv);

    if (r < 0)
        r = d;

    conf.rare_event.start_level = (d-1)/2 - 1;
    conf.rare_event.max_level = 128;

    auto circuit = si1000(d, r, p, false);

    // Convert to DEM with error decomposition required by PyMatching
    auto dem = stim::circuit_to_dem(circuit, {.decompose_errors = true});

    // Build decoder, run estimation, and report results.
    auto run = [&] (auto&& dec, const auto& error_callback)
    {
        double ler = estimate_logical_error_rate(dem, dec, conf, error_callback);
        dec.mpi_accumulate();
        if (world_rank == 0)
        {
            std::cout << "Logical error rate: " << ler << "\n";
            dec.print_stats(std::cout);
        }
    };

    if (decoder_name == "pymatching")
    {
        run(decoder::PyMatching(dem, enable_gap), [] (auto, auto, auto) {});
    }
    else if (decoder_name == "blossom5")
    {
        run(decoder::BlossomV(dem), [] (auto, auto, auto) {});
    }
    else
    {
        std::cerr << "unknown decoder name: " << decoder_name << _die{};
    }

#if defined(ENABLE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    return 0;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
