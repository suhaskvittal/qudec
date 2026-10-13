#include "circuit_generator.h"
#include "decoder/libra.h"
#include "sampler.h"
#include <stim/util_top/circuit_to_dem.h>
#include <iostream>
#include <stdexcept>

namespace {
void run(int rank, int workers) {
    const auto circuit = sc_si1000(3, 3, 1e-3, false, true);
    const auto dem = stim::circuit_to_dem(circuit, {.decompose_errors = true});
    dec::Libra::Config cfg;
    cfg.ensemble_size = 2;
    cfg.seed = 20261001;
    cfg.validate = true;
    dec::Libra decoder(dem, cfg);
    ExperimentConfig experiment;
    experiment.monte_carlo.max_samples = 17;
    experiment.monte_carlo.stop_at_error_count = 1000;
    const double ler = estimate_logical_error_rate(dem, decoder, experiment);
    if (ler < 0 || ler > 1) throw std::runtime_error("Invalid Monte Carlo estimate");
    decoder.mpi_accumulate();
    if (decoder.statistics().shots != 17)
        throw std::runtime_error("Libra MPI statistics or global Monte Carlo budget is incorrect");

    // Every rank takes the fallback once, then check all counters are reduced.
    cfg.gap_threshold_db = 0;
    dec::Libra fallback(dem, cfg);
    dec::SyndromeType syndrome(dem.count_detectors()), truth(1);
    syndrome[0] = 1;
    const auto result = fallback.decode(syndrome, truth);
    if (result.flipped_obs[0] != result.matching_data.frame_flips[0])
        throw std::runtime_error("Generic result frame fields disagree");
    fallback.mpi_accumulate();
    const auto statistics = fallback.statistics();
    if (statistics.shots != uint64_t(workers) || statistics.fallback_shots != uint64_t(workers) ||
        statistics.ensemble_shots != 0)
        throw std::runtime_error("Libra fallback counter reduction is incorrect");
    if (!rank) std::cout << "Integrated Libra interface and MPI statistics passed\n";
}
}
int main(int argc, char** argv) {
    int rank = 0, workers = 1;
#if defined(ENABLE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &workers);
#endif
    int status = 0;
    try { run(rank, workers); }
    catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
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
