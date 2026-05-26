/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#include "sampler.h"
#include "decoder/surface_code.h"

#include "stim/gen/gen_surface_code.h"
#include "stim/util_top/circuit_to_dem.h"

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
    int d = atoi(argv[1]);

    // Generate d=11 rotated surface code memory-Z experiment (11 rounds, p=0.1%)
    stim::CircuitGenParameters params(d, d, "rotated_memory_z");
    params.after_clifford_depolarization = 0.001;
    auto gen = stim::generate_surface_code_circuit(params);

    // Convert to DEM with error decomposition required by PyMatching
    auto dem = stim::circuit_to_dem(gen.circuit, {.decompose_errors = true});

    // Build decoder and run estimation
//  decoder::PYMATCHING dec(dem);
    decoder::CLUSTER_MATCH dec(dem, d, 10, decoder::CLUSTER_MATCH::quantization_level::b32);

    EXPERIMENT_CONFIG conf{.samples_per_level=10000 };
    conf.start_level = (d-1)/2 - 1;
    conf.max_level = 128;
    conf.verbosity = 0;
    conf.print_progress = true;
    double ler = estimate_logical_error_rate(dem, dec, conf);

    std::cout << "Logical error rate: " << ler << "\n";
#if defined(ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
