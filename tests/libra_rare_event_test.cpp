/*
 *  contributor: OpenAI Codex
 * */

#include "libra_rare_event.h"

#include <array>
#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <iostream>

void check(bool value) {
    if (!value) throw std::runtime_error("conditional sampler oracle failed");
}

int main()
{
    const stim::DetectorErrorModel dem(
        "error(0.1) D0 L0\n"
        "error(0.4) D1 L0\n"
        "error(0.7) D2 L0\n");
    LibraRareEventSampler sampler(dem, 4);
    const auto& pk = sampler.count_probabilities();

    const std::array<double, 4> expected{
        0.9 * 0.6 * 0.3,
        0.1 * 0.6 * 0.3 + 0.9 * 0.4 * 0.3 + 0.9 * 0.6 * 0.7,
        0.1 * 0.4 * 0.3 + 0.1 * 0.6 * 0.7 + 0.9 * 0.4 * 0.7,
        0.1 * 0.4 * 0.7};
    for (size_t k = 0; k < expected.size(); ++k)
        check(std::abs(pk[k] - expected[k]) < 1e-14);

    // At fixed K=1, instruction odds are p/(1-p). Count samples should follow
    // the normalized odds (1/9, 2/3, 7/3).
    const std::array<double, 3> odds{1.0 / 9.0, 2.0 / 3.0, 7.0 / 3.0};
    const double odds_sum = odds[0] + odds[1] + odds[2];
    std::array<uint64_t, 3> seen{};
    std::mt19937_64 rng(20261001);
    constexpr uint64_t trials = 300000;
    for (uint64_t i = 0; i < trials; ++i)
    {
        const auto selected = sampler.sample_errors(1, rng);
        check(selected.size() == 1 && selected[0] < 3);
        ++seen[selected[0]];
    }
    for (size_t i = 0; i < seen.size(); ++i)
    {
        const double observed = static_cast<double>(seen[i]) / trials;
        const double target = odds[i] / odds_sum;
        check(std::abs(observed - target) < 0.004);
    }

    for (size_t k = 0; k <= 3; ++k)
    {
        const auto selected = sampler.sample_errors(k, rng);
        check(selected.size() == k);
    }
    std::cout << "exact conditional sampler enumeration check passed\n";
}
