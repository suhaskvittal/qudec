/*
 *  contributor: OpenAI Codex
 * */

#include "sampler.h"
#include "decoder/types.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
struct AlwaysOne {
    uint64_t calls = 0;
    dec::Result decode(dec::SyndromeRef, dec::ObsRef) {
        ++calls;
        dec::Result result;
        result.flipped_obs[0] = 1;
        return result;
    }
};
void run(int rank) {
    // All channels flip L0, so error conditional on K is deterministic for
    // this decoder. K=0 is essential; starting at the mode biases the sum.
    const stim::DetectorErrorModel dem("error(0.1) D0 L0\nerror(0.4) D1 L0\n");
    AlwaysOne decoder;
    ExperimentConfig config;
    config.method = "rare_event";
    config.rare_event.max_level = 3;
    config.rare_event.samples_per_level = 17;
    config.rare_event.max_errors_per_level = 1000;
    const double estimate = estimate_logical_error_rate(dem, decoder, config);
    require(std::abs(estimate - 0.58) < 1e-13, "Full count sum disagrees with exact logical error rate");
    uint64_t calls = decoder.calls;
#if defined(ENABLE_MPI)
    MPI_Allreduce(MPI_IN_PLACE, &calls, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
    require(calls == 51, "Per-level global sample limit was not respected");
    config.rare_event.max_level = 2;
    bool rejected = false;
    try { estimate_logical_error_rate(dem, decoder, config); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "Unsafe truncated count range was accepted");
    config.rare_event.max_level = 3;
    config.rare_event.start_level = 1;
    rejected = false;
    try { estimate_logical_error_rate(dem, decoder, config); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "Skipping low-count levels was accepted");

    const stim::DetectorErrorModel unequal("error(0.1) D0 L0\nerror(0.4) D1 L0\nerror(0.7) D2 L0\n");
    Rng rng(4242 + rank);
    constexpr size_t trials = 30000;
    const auto [syndromes, truths] = generate_syndromes_with_k_errors(unequal, 1, trials, rng);
    std::array<size_t, 3> seen{};
    for (size_t i = 0; i < trials; ++i) {
        unsigned count = 0;
        for (size_t j = 0; j < 3; ++j) { count += bool(syndromes[i][j]); seen[j] += bool(syndromes[i][j]); }
        require(count == 1 && truths[i][0], "Public fixed-count generator returned inconsistent parity");
    }
    const std::array<double, 3> odds{1.0 / 9, 2.0 / 3, 7.0 / 3};
    const double total = odds[0] + odds[1] + odds[2];
    for (size_t j = 0; j < 3; ++j)
        require(std::abs(double(seen[j]) / trials - odds[j] / total) < 0.02,
                "Public generator does not use heterogeneous Bernoulli probabilities");

    ConditionalErrorSampler sampler(unequal, 4);
    const auto& pk = sampler.count_probabilities();
    for (size_t k = 0; k <= 4; ++k) {
        double exact_tail = 0;
        for (size_t j = k; j < pk.size(); ++j) exact_tail += pk[j];
        const double bound = sampler.tail_probability_bound(k);
        require(bound >= exact_tail - 1e-14 && bound <= 1,
                "Count-tail bound is smaller than the enumerated tail");
    }
    const stim::DetectorErrorModel deterministic("error(0) D0\nerror(1) D1 L0\n");
    const auto [ds, os] = generate_syndromes_with_k_errors(deterministic, 1, 8, rng);
    for (size_t i = 0; i < 8; ++i)
        require(!ds[i][0] && ds[i][1] && os[i][0], "Deterministic channels sampled incorrectly");
    const stim::DetectorErrorModel repeated("repeat 2 {\n error(0.2) D0 L0\n shift_detectors 1\n}\n");
    const auto [rs, ro] = generate_syndromes_with_k_errors(repeated, 2, 8, rng);
    for (size_t i = 0; i < 8; ++i)
        require(rs[i][0] && rs[i][1] && !ro[i][0], "Flattened instruction IDs or observable XOR is incorrect");
    const stim::DetectorErrorModel cancellation("error(0.2) D0 ^ D0 L0\n");
    const auto [cs, co] = generate_syndromes_with_k_errors(cancellation, 1, 8, rng);
    for (size_t i = 0; i < 8; ++i)
        require(!cs[i][0] && co[i][0], "Separated targets were not combined with XOR");
    require(compute_probability_polynomial(dem, 0).empty(), "Zero coefficient limit should return an empty polynomial");
    const stim::DetectorErrorModel empty;
    const auto empty_samples = generate_syndromes_with_k_errors(empty, 0, 8, rng);
    (void)empty_samples;
    if (!rank) std::cout << "Public conditional sampler, complete count sum, tail bounds, and MPI budgets passed\n";
}
}
int main(int argc, char** argv) {
    int rank = 0;
#if defined(ENABLE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    int status = 0;
    try { run(rank); }
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
