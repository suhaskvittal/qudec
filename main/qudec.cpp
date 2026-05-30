/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#include "sampler.h"
#include "decoder/surface_code.h"

#include <stim/gen/gen_surface_code.h>
#include <stim/util_top/circuit_to_dem.h>

#include <argparse/argparse.h>

#include <iostream>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

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
    int64_t d;
    EXPERIMENT_CONFIG conf;

    /*
     * CLUSTER-MATCH parameters:
     * */
    int64_t cm_astrea_hw_max;

    ARGPARSE()
        .required("code-distance", "Code distance of surface code", d)
        .optional("-v", "--verbose", "Verbosity level", conf.verbosity, 0)
        .optional("-s", "--samples", "Samples per error level", conf.samples_per_level, 10000)
        .optional("-pp", "--print-progress", "Print simulation progress", conf.print_progress, false)
        .optional("", "--cm-astrea-hw-max", "Max HW supported by Astrea decoder", cm_astrea_hw_max, 8)
        .parse(argc, argv);

    conf.start_level = (d-1)/2 - 1;
    conf.max_level = 128;

    // Generate d=11 rotated surface code memory-Z experiment (11 rounds, p=0.1%)
    stim::CircuitGenParameters params(d, d, "rotated_memory_z");
    double p = 1e-3;
    params.after_clifford_depolarization = p;
    params.before_round_data_depolarization = p;
    params.before_measure_flip_probability = p;
    params.after_reset_flip_probability = p;
    auto gen = stim::generate_surface_code_circuit(params);

    // Convert to DEM with error decomposition required by PyMatching
    auto dem = stim::circuit_to_dem(gen.circuit, {.decompose_errors = true});

    // Build decoder and run estimation
//  decoder::PYMATCHING dec(dem);
    decoder::CLUSTER_MATCH dec(dem, d, cm_astrea_hw_max, decoder::CLUSTER_MATCH::quantization_level::b16);

    double ler = estimate_logical_error_rate(dem, dec, conf);

    if (world_rank == 0)
    {
        std::cout << "Logical error rate: " << ler << "\n";
        dec.print_stats(std::cout);
    }

#if defined(ENABLE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    return 0;
}
