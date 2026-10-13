/*
 *  contributor: OpenAI Codex
 * */

#include "libra_local_splitting.h"
#include <array>
#include <iostream>

int main() {
    const stim::DetectorErrorModel dem(
        "error(0.02) D0 L0\n"
        "error(0.11) D0 D1\n"
        "error(0.23) D1 D2\n"
        "error(0.07) D2\n");
    LibraLocalSplittingModel model(dem);
    auto event = [](dec::SyndromeRef syn, dec::ObsRef) {
        return LibraSplittingModel::Labels{bool(syn[0]), bool(syn[0])};
    };
    constexpr double theta = 0.5;
    std::array<double, 4> expected{}, seen{};
    double event_mass = 0, original_mass = 0;
    for (unsigned mask = 0; mask < 16; ++mask) {
        if (bool(mask & 1) == bool(mask & 2))
            continue; // Exactly one D0 channel must be active.
        double tilted = 1, original = 1;
        for (unsigned j = 0; j < 4; ++j) {
            const double p = model.probabilities[j];
            const double odds = model.odds[j] * std::exp(theta);
            const double pt = odds / (1 + odds);
            tilted *= mask >> j & 1 ? pt : 1 - pt;
            original *= mask >> j & 1 ? p : 1 - p;
        }
        event_mass += tilted;
        original_mass += original;
        for (unsigned j = 0; j < 4; ++j)
            if (mask >> j & 1)
                expected[j] += tilted;
    }
    for (double& p : expected)
        p /= event_mass;
    auto state = model.make_state({0}, event);
    std::mt19937_64 rng(20261003);
    constexpr unsigned burnin = 20000, records = 1000000;
    double ratio_sum = 0;
    for (unsigned j = 0; j < burnin + records; ++j) {
        model.step(state, theta, rng, event);
        if (j % 1024 == 0) {
            std::array<bool, 3> parity{};
            bool observable = false;
            double active_odds = 0;
            for (size_t id : state.active) {
                if (state.position[id] < 0 || state.active[state.position[id]] != id)
                    throw std::runtime_error("local-swap active-position invariant failed");
                for (uint32_t d : model.patterns[id].detectors)
                    parity[d] = !parity[d];
                observable ^= model.patterns[id].observable;
                active_odds += model.odds[id];
            }
            for (unsigned d = 0; d < 3; ++d)
                if (parity[d] != bool(state.syndrome[d]))
                    throw std::runtime_error("local-swap detector parity invariant failed");
            if (observable != bool(state.truth[0]) ||
                std::abs(active_odds - state.active_odds) > 1e-9)
                throw std::runtime_error("local-swap truth/odds invariant failed");
        }
        if (j < burnin)
            continue;
        for (unsigned i = 0; i < 4; ++i)
            seen[i] += state.position[i] >= 0;
        ratio_sum += std::exp(model.log_ratio(theta, 0, state.active.size()));
    }
    for (unsigned i = 0; i < 4; ++i)
        if (std::abs(seen[i] / records - expected[i]) > 0.005)
            throw std::runtime_error("local-swap conditional distribution oracle failed");
    if (std::abs(ratio_sum / records - original_mass / event_mass) > 0.005)
        throw std::runtime_error("local-swap likelihood ratio oracle failed");
    std::cout << "local-swap exact conditional distribution and bridge ratio passed\n";
}
