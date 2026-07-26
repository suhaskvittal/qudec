/*
 *  author: Suhas Vittal
 *  date:   19 June 2026
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

    int64_t d;
    ExperimentConfig conf;

    bool enable_astrea{false};
    bool enable_filter{false};

    ARGPARSE()
        .required("code-distance", "Code distance of surface code", d)
        .optional("-s", "--samples", "Samples per error level", conf.samples_per_level, 10000)
        .optional("-pp", "--print-progress", "Print simulation progress", conf.print_progress, false)
        .optional("", "--astrea", "Enable Astrea HW emulation", enable_astrea, false)
        .optional("", "--filter", "Enable filter HW emulation", enable_filter, false)
        .parse(argc, argv);

    conf.start_level = (d-1)/2 - 1;
    conf.max_level = 128;

    stim::CircuitGenParameters params(d, d, "rotated_memory_z");
    double p = 1e-3;
    params.after_clifford_depolarization = p;
    params.before_round_data_depolarization = p;
    params.before_measure_flip_probability = p;
    params.after_reset_flip_probability = p;
    auto gen = stim::generate_surface_code_circuit(params);
    auto dem = stim::circuit_to_dem(gen.circuit, {.decompose_errors = true});

    uint8_t hw_emu_enable = decoder::ClusterMatch::hw_emu_flag::validate;
    if (enable_astrea) hw_emu_enable |= decoder::ClusterMatch::hw_emu_flag::astrea;
    if (enable_filter) hw_emu_enable |= decoder::ClusterMatch::hw_emu_flag::filter;

    decoder::ClusterMatch dec(dem, d, 6,
        decoder::ClusterMatch::quantization_level::b16,
        hw_emu_enable);

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
