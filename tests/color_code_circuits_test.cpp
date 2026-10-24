// Implementation authored by OpenAI Codex.
#include "circuit_generator.h"

#include <stim/simulators/frame_simulator_util.h>
#include <stim/simulators/tableau_simulator.h>
#include <stim/util_top/circuit_to_dem.h>

#include <cmath>
#include <functional>
#include <bit>
#include <iostream>
#include <iomanip>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
using Generator = std::function<stim::Circuit(uint32_t, uint32_t, double, bool, const ComponentErrorDistribution&)>;

auto detector_coordinates(const stim::Circuit& c) {
    std::set<uint64_t> indices;
    for (uint64_t k = 0; k < c.count_detectors(); ++k) indices.insert(k);
    return c.get_detector_coordinates(indices);
}

void check_geometry(const stim::Circuit& circuit, uint32_t n, uint32_t distance, bool oct) {
    std::vector<std::set<uint32_t>> supports;
    std::vector<int> colors;
    const auto coords = detector_coordinates(circuit);
    size_t detector = 0;
    circuit.for_each_operation([&](const stim::CircuitInstruction& op) {
        if (op.gate_type != stim::GateType::DETECTOR) return;
        const auto& pos = coords.at(detector++);
        require(pos.size() == 6 && pos[3] >= 0 && pos[3] <= 2 &&
                (pos[4] == 0 || pos[4] == 1) && (pos[5] == 0 || pos[5] == 1),
                "Invalid detector metadata");
        if (pos[5] == 1) require(op.targets.size() == 1, "Flag compares multiple rounds");
        // Final detectors have a data record for each support vertex.
        if (pos[2] != 1 || pos[5] != 0 || op.targets.size() < 4) return;
        std::set<uint32_t> support;
        for (auto t : op.targets)
            if (t.rec_offset() >= -(int)n) support.insert(n + t.rec_offset());
        require(support.size() == 4 || support.size() == (oct ? 8 : 6), "Wrong face weight");
        if (oct && pos[3] == 0) require(support.size() == 4, "Red 488 face is not square");
        supports.push_back(std::move(support));
        colors.push_back(pos[3]);
    });
    require(supports.size() == (n - 1) / 2, "Missing final check supports");
    std::vector<uint64_t> rows;
    for (size_t i = 0; i < supports.size(); ++i) {
        uint64_t row = 0;
        for (auto q : supports[i]) row |= uint64_t{1} << q;
        rows.push_back(row);
        for (size_t j = 0; j < i; ++j) {
            size_t overlap = 0;
            for (auto q : supports[i]) overlap += supports[j].contains(q);
            require(overlap % 2 == 0, "X and Z checks do not commute");
            if (colors[i] == colors[j]) require(overlap == 0, "Same-color faces overlap");
        }
    }
    // Independent face generators imply n - 2 rank(H) = 1 logical qubit.
    size_t rank = 0;
    for (uint32_t q = 0; q < n; ++q) {
        size_t pivot = rank;
        while (pivot < rows.size() && !(rows[pivot] >> q & 1)) ++pivot;
        if (pivot == rows.size()) continue;
        std::swap(rows[rank], rows[pivot]);
        for (size_t j = rank + 1; j < rows.size(); ++j)
            if (rows[j] >> q & 1) rows[j] ^= rows[rank];
        ++rank;
    }
    require(2 * rank + 1 == n, "Patch does not encode one logical qubit");
    // Exhaustively establish the CSS distance for the small patches. All-data
    // logical Z/X means a commuting odd-weight Pauli is a logical operator.
    if (distance <= 5) {
        uint32_t smallest = n;
        for (uint64_t mask = 1; mask < (uint64_t{1} << n); ++mask) {
            uint32_t weight = std::popcount(mask);
            if (weight % 2 == 0 || weight >= smallest) continue;
            bool commutes = true;
            for (auto row : rows) commutes &= (std::popcount(row & mask) % 2 == 0);
            if (commutes) smallest = weight;
        }
        require(smallest == distance, "Wrong code distance");
    }
}

void check_absolute_parities(const stim::Circuit& circuit) {
    // Detector samplers subtract a reference; explicitly catch constant -1
    // detector/observable parities as well as random detection events.
    auto reference = stim::TableauSimulator<64>::reference_sample_circuit(circuit);
    uint64_t count = 0;
    circuit.for_each_operation([&](const stim::CircuitInstruction& op) {
        count += op.count_measurement_results();
        if (op.gate_type != stim::GateType::DETECTOR &&
            op.gate_type != stim::GateType::OBSERVABLE_INCLUDE) return;
        bool parity = false;
        for (auto target : op.targets) parity ^= reference[count + target.rec_offset()];
        require(!parity, "Nonzero absolute noiseless parity");
    });
}

void check_schedule_and_noise(const stim::Circuit& circuit, double p) {
    std::set<uint32_t> busy;
    std::set<stim::GateType> channels;
    circuit.for_each_operation([&](const stim::CircuitInstruction& op) {
        if (op.gate_type == stim::GateType::TICK) busy.clear();
        if (op.gate_type == stim::GateType::R || op.gate_type == stim::GateType::H ||
            op.gate_type == stim::GateType::CX || op.gate_type == stim::GateType::M ||
            op.gate_type == stim::GateType::MR)
            for (auto t : op.targets)
                require(busy.insert(t.qubit_value()).second, "Gate collision within a tick");
        auto near = [](double x, double y) { return std::abs(x - y) < 1e-12; };
        if (op.gate_type == stim::GateType::DEPOLARIZE1) {
            require(near(op.args[0], p / 10) || near(op.args[0], p) || near(op.args[0], 2 * p),
                    "Wrong SI1000 single-qubit depolarization");
            channels.insert(op.gate_type);
        } else if (op.gate_type == stim::GateType::DEPOLARIZE2) {
            require(near(op.args[0], p), "Wrong SI1000 two-qubit depolarization");
            channels.insert(op.gate_type);
        } else if (op.gate_type == stim::GateType::X_ERROR) {
            require(near(op.args[0], 2 * p) || near(op.args[0], 5 * p), "Wrong reset/readout error rate");
            channels.insert(op.gate_type);
        }
    });
    require(channels.size() == 3, "Missing SI1000 noise channels");
}

void check_flag_independence(Generator generate, int basis) {
    auto original = generate(3, 2, 0, false, {});
    const auto coords = detector_coordinates(original);
    // Flip one flag's readout bit in round zero. An accidental consecutive-
    // round flag XOR would fire twice. Test X and Z extraction separately.
    uint64_t record = 0, count = 0, detector = 0;
    bool found = false;
    original.for_each_operation([&](const stim::CircuitInstruction& op) {
        count += op.count_measurement_results();
        if (op.gate_type != stim::GateType::DETECTOR) return;
        const auto& pos = coords.at(detector++);
        if (!found && pos[2] == 0 && pos[4] == basis && pos[5] == 1) {
            record = count + op.targets[0].rec_offset();
            found = true;
        }
    });
    require(found, "Missing flag detector");
    bool inserted = false;
    count = 0;
    stim::Circuit faulty;
    original.for_each_operation([&](const stim::CircuitInstruction& op) {
        if (!inserted && (op.gate_type == stim::GateType::M || op.gate_type == stim::GateType::MR))
            if (count <= record && record < count + op.targets.size()) {
                faulty.safe_append_ua("X_ERROR", {op.targets[record - count].qubit_value()}, 1);
                inserted = true;
            }
        count += op.count_measurement_results();
        faulty.safe_append(op);
    });
    std::mt19937_64 rng(19);
    auto [dets, obs] = stim::sample_batch_detection_events<64>(faulty, 64, rng);
    size_t fired = 0;
    for (size_t k = 0; k < original.count_detectors(); ++k) if (dets[k].not_zero()) {
        ++fired;
        require(coords.at(k)[5] == 1 && coords.at(k)[2] == 0 && coords.at(k)[4] == basis,
                "Flag fault leaked into another round/check");
    }
    require(inserted && fired == 1 && !obs.data.not_zero(), "Single flag readout fault should fire exactly once");
}

// Author: OpenAI Codex. Physical concurrency and resonator-wait budget.
void check_parallel_rounds(const stim::Circuit& c,uint32_t n,uint32_t rounds,double p) {
    size_t resets=0,transitions=0;
    std::vector<size_t> waits(n);
    c.for_each_operation([&](const auto& op) {
        if (op.gate_type==stim::GateType::R) {
            ++resets;
            require(op.targets.size()==2*n-1,"Initialization must reset data and all Bell ancillas together");
        }
        if (op.gate_type==stim::GateType::MR) {
            ++transitions;
            require(op.targets.size()==n-1,"Basis transition must measure/reset every color together");
            std::set<uint32_t> targets;
            for (auto t:op.targets) targets.insert(t.qubit_value());
            for (uint32_t q=n;q<2*n-1;q++) require(targets.contains(q),"Ancilla missing from joint basis transition");
        }
        if (op.gate_type==stim::GateType::DEPOLARIZE1 && std::abs(op.args[0]-2*p)<1e-12)
            for (auto t:op.targets) if (t.qubit_value()<n) ++waits[t.qubit_value()];
    });
    require(resets==1 && transitions==2*rounds,"Extra reset/readout subrounds in Bell extraction");
    for (auto count:waits) require(count==2*rounds,"Each data qubit must have exactly two 2p waits per round");
}
} // namespace

int main() {
    std::mt19937_64 rng(42);
    for (bool oct : {false, true}) {
        Generator generate = [oct](uint32_t d, uint32_t r, double p, bool x,
                                   const ComponentErrorDistribution& noise) {
            return oct ? cc488_si1000(d,r,p,x,true,noise) : cc666_si1000(d,r,p,x,false,noise);
        };
        check_flag_independence(generate, 0);
        check_flag_independence(generate, 1);
        for (uint32_t d : {3, 5, 7, 9, 11, 15}) for (bool x : {false, true})
            for (uint32_t rounds : {1, 2, 5}) {
                const uint32_t n = oct ? (d * d + 2 * d - 1) / 2 : (3 * d * d + 1) / 4;
                auto c = generate(d, rounds, 0, x, {});
                require(c.count_qubits() == 2 * n - 1, "Wrong qubit count");
                // Contributor: OpenAI Codex. There are (n-1)/2 faces; each
                // basis contributes one check/flag detector per face.
                require(c.count_detectors() == (n - 1) * (oct ? 4 * rounds : 3 * rounds + 1) / 2, "Wrong detector count");
                require(c.count_measurements() == 2 * (n - 1) * rounds + n, "Wrong measurement count");
                check_absolute_parities(c);
                auto [dets, obs] = stim::sample_batch_detection_events<64>(c, 1024, rng);
                require(!dets.data.not_zero() && !obs.data.not_zero(), "Noiseless detection/logical event");
                (void)stim::circuit_to_dem(c, {});
                if (rounds == 1 && d <= 7) check_geometry(c, n, d, oct);
                // No graphlike decomposition: color-code errors are hyperedges.
                auto noisy = generate(d, rounds, 0.001, x, {});
                check_schedule_and_noise(noisy, 0.001);
                check_parallel_rounds(noisy,n,rounds,0.001);
                auto dem = stim::circuit_to_dem(noisy, {});
                require(dem.count_errors() > 0, "No SI1000 errors in noisy model");
                // Larger honeycomb coordinates need explicit stream precision;
                // Stim's default six significant digits round sqrt(3) positions.
                std::ostringstream text;
                text<<std::setprecision(17)<<noisy;
                require(stim::Circuit(text.str()).approx_equals(noisy,1e-12), "Stim text round trip failed");
            }
        for (double p : {-0.1, 0.201, std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity()}) {
            bool rejected = false;
            try { (void)generate(3, 1, p, false, {}); } catch (const std::invalid_argument&) { rejected = true; }
            require(rejected, "Invalid noise rate accepted");
        }
        for (auto [d, r] : {std::pair{2u, 1u}, {4u, 1u}, {3u, 0u}}) {
            bool rejected = false;
            try { (void)generate(d, r, 0, false, {}); } catch (const std::invalid_argument&) { rejected = true; }
            require(rejected, "Invalid geometry accepted");
        }
        (void)generate(3, 1, 0.2, false, {});
    }
    std::cout << "Color code circuits: geometry, distance, flag independence, SI1000 DEMs, "
                 "joint basis transitions, two waits per round, and 73728 zero-noise shots passed\n";
}
