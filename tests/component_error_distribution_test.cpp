// Implementation authored by OpenAI Codex.
#include "circuit_generator.h"
#include <stim/util_top/circuit_to_dem.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
void require(bool yes, const char* message) {
    if (!yes) throw std::runtime_error(message);
}
}
int main() {
    ComponentErrorDistribution distribution{0.7, 123};
    const stim::Circuit original("X_ERROR(0.001) 0 1\nTICK\nREPEAT 3 {\nX_ERROR(0.001) 0 1\nDEPOLARIZE2(0.002) 0 1\n}\n");
    require(apply_component_error_distribution(original, {}) == original, "Zero spread changed circuit");
    const auto noisy = apply_component_error_distribution(original, distribution);
    require(noisy == apply_component_error_distribution(original, distribution), "Seed not reproducible");
    require(noisy != apply_component_error_distribution(original, {0.7, 124}), "Seeds have same rates");
    const auto rescaled = apply_component_error_distribution(
        stim::Circuit("X_ERROR(0.002) 0 1\n"), distribution);
    require(std::abs(rescaled.operations[0].args[0] / noisy.operations[0].args[0] - 2) < 1e-12,
            "Changing nominal p resamples device disorder");
    std::vector<double> rates;
    noisy.for_each_operation([&](const stim::CircuitInstruction& op) {
        if (op.gate_type == stim::GateType::X_ERROR) rates.push_back(op.args[0]);
    });
    require(rates.size() == 8 && rates[0] != rates[1], "Qubits share component draw");
    for (size_t k = 2; k < rates.size(); ++k)
        require(rates[k] == rates[k%2], "Component draw changes across rounds");
    require(noisy.blocks.size() == 1, "REPEAT block flattened");
    require(apply_component_error_distribution(original.flattened(), distribution).flattened() == noisy.flattened(),
            "Rates depend on repeat expansion");
    long double mean = 0, log_mean = 0, log_square = 0;
    constexpr size_t n = 100000;
    for (size_t k = 0; k < n; ++k) {
        const double rate = distribution.probability(1e-6, k);
        const double log_ratio = std::log(rate / 1e-6);
        mean += rate;
        log_mean += log_ratio;
        log_square += log_ratio*log_ratio;
    }
    mean /= n; log_mean /= n; log_square /= n;
    require(std::abs(mean / 1e-6 - 1) < 0.015, "Lognormal arithmetic mean shifted");
    require(std::abs(log_mean + 0.245) < 0.015, "Lognormal mean incorrect");
    require(std::abs(log_square-log_mean*log_mean - 0.49) < 0.015, "Lognormal spread incorrect");
    for (double sigma : {-1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        bool failed = false;
        try { (void)apply_component_error_distribution(original, {sigma, 0}); }
        catch (const std::invalid_argument&) { failed = true; }
        require(failed, "Invalid sigma accepted");
    }
    require(distribution.probability(0, 0) == 0, "Zero physical noise became nonzero");
    for (size_t k = 0; k < 1000; ++k) {
        double rate = ComponentErrorDistribution{3, 1}.probability(0.5, k);
        require(rate >= 0 && rate <= 1, "Unphysical probability");
    }
    bool saturated_one = false, saturated_two = false;
    for (uint64_t seed = 0; seed < 100; ++seed) {
        const auto saturated = apply_component_error_distribution(
            stim::Circuit("R 0 1\nDEPOLARIZE1(0.5) 0\nDEPOLARIZE2(0.5) 0 1\nM 0\nDETECTOR rec[-1]\n"),
            {2, seed});
        for (const auto& op : saturated.operations) {
            if (op.gate_type == stim::GateType::DEPOLARIZE1) {
                require(op.args[0] <= 0.75, "Overmixing single-qubit channel");
                saturated_one |= op.args[0] == 0.75;
            }
            if (op.gate_type == stim::GateType::DEPOLARIZE2) {
                require(op.args[0] <= 0.9375, "Overmixing two-qubit channel");
                saturated_two |= op.args[0] == 0.9375;
            }
        }
        (void)stim::circuit_to_dem(saturated, {});
    }
    require(saturated_one && saturated_two, "High-spread channel saturation not exercised");
    for (bool memory_x : {false, true}) {
        const auto sc = sc_si1000(3, 3, 0.001, memory_x, true);
        require(sc_si1000(3, 3, 0.001, memory_x, true, distribution) ==
                apply_component_error_distribution(sc, distribution), "Surface generator ignores distribution");
        const auto toric = toric_si1000(3, 3, 0.001, memory_x, true);
        require(toric_si1000(3, 3, 0.001, memory_x, true, distribution) ==
                apply_component_error_distribution(toric, distribution), "Toric generator ignores distribution");
        for (auto color : {cc666_si1000, cc488_si1000}) {
            const auto c = color(3, 3, 0.001, memory_x, distribution);
            require(c == apply_component_error_distribution(color(3, 3, 0.001, memory_x, {}), distribution),
                    "Color generator ignores distribution");
            require(stim::circuit_to_dem(c, {}).count_errors() > 0, "Nonuniform color DEM failed");
        }
    }
    std::cout << "Component noise: fixed rates, seeds, zero spread, lognormal moments, and all SI1000 generators passed\n";
}
