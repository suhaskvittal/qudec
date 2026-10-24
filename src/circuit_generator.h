/*
 *  author: Claude Sonnet 5.1
 *  contributor: OpenAI Codex
 *  date:   26 July 2026
 * */

#ifndef CIRCUIT_GENERATOR_h
#define CIRCUIT_GENERATOR_h

#include <stim.h>
#include "component_error_distribution.h"
#include "circuit_noise_policy.h"

#include <cstdint>
#include <type_traits>
#include <utility>

// Model-independent entry points. Supply a policy defining each operation and
// its injected errors. The model-named functions below are compatibility adapters.
stim::Circuit sc_circuit(uint32_t distance, uint32_t rounds,
    const CircuitNoisePolicy& policy, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {});
stim::Circuit toric_circuit(uint32_t distance, uint32_t rounds,
    const CircuitNoisePolicy& policy, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {});
stim::Circuit cc666_circuit(uint32_t distance, uint32_t rounds,
    const CircuitNoisePolicy& policy, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {});
stim::Circuit cc488_circuit(uint32_t distance, uint32_t rounds,
    const CircuitNoisePolicy& policy, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {});

// Author: OpenAI Codex. Superdense extraction, arXiv:2609.21376v3.
// Both bases are acquired per reset/measure cycle. 6.6.6: odd d>=3;
// brickwork 4.8.8: d=4m-1 only. Policies must implement gate sequences.
// CNOT is the default; true selects native CZ compilation.
// Detector coordinates: (x,y,round,color,basis,is_flag), flags color=-1.
stim::Circuit cc666_superdense(uint32_t distance, uint32_t rounds,
    const CircuitNoisePolicy& policy, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {}, bool native_cz = false);
stim::Circuit cc488_superdense(uint32_t distance, uint32_t rounds,
    const CircuitNoisePolicy& policy, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {}, uint32_t start_phase = 0,
    bool native_cz = false);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Generates a rotated surface-code syndrome-extraction circuit under the
 * `SI1000` circuit-level noise model:
 *      (1) single-qubit gates    -> DEPOLARIZE1(p/10) after the gate,
 *      (2) two-qubit gates       -> DEPOLARIZE2(p) after the gate,
 *      (3) reset gates           -> X_ERROR(2p) after the gate,
 *      (4) measurement gates     -> X_ERROR(5p) before the gate (and, for a
 *                                   pure measurement, DEPOLARIZE1(p) after),
 *      (5) idle during a CX layer -> DEPOLARIZE1(p),
 *      (6) idle during measurement -> DEPOLARIZE1(2p).
 *
 * By default, detection events use only checks in the memory basis: Z checks
 * for Z memory and X checks for X memory. With
 * `include_opposite_basis_detectors`, comparisons between successive interior
 * rounds also emit detectors for the other check basis. The first round and
 * final data readout retain only the memory-basis detectors.
 *
 * Arguments:
 *      `distance`    -- code distance of the rotated surface code,
 *      `rounds`      -- number of syndrome-extraction rounds,
 *      `p`           -- physical error rate,
 *      `is_memory_x` -- true for an X-memory experiment, false for Z-memory.
 *      `include_opposite_basis_detectors` -- include other-basis interior checks.
 * */
// OpenAI GPT-6-Sol: Expose optional interior opposite-basis detectors while preserving existing calls.
stim::Circuit sc_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x,
                        bool include_opposite_basis_detectors = false,
                        const ComponentErrorDistribution& noise_distribution = {});

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Generates a toric-code syndrome-extraction circuit (periodic boundaries in
 * both directions) under the same `SI1000` circuit-level noise model as
 * `sc_si1000` above:
 *      (1) single-qubit gates    -> DEPOLARIZE1(p/10) after the gate,
 *      (2) two-qubit gates       -> DEPOLARIZE2(p) after the gate,
 *      (3) reset gates           -> X_ERROR(2p) after the gate,
 *      (4) measurement gates     -> X_ERROR(5p) before the gate (and, for a
 *                                   pure measurement, DEPOLARIZE1(p) after),
 *      (5) idle during a CX layer -> DEPOLARIZE1(p),
 *      (6) idle during measurement -> DEPOLARIZE1(2p).
 *
 * `distance` gives the linear size L of the L x L toric lattice, for 2*L^2
 * data qubits and 2*L^2 ancilla qubits (L^2 X-check vertices, L^2 Z-check
 * plaquettes). The toric code has two logical qubits; both observables are
 * emitted (`OBSERVABLE_INCLUDE(0)` and `OBSERVABLE_INCLUDE(1)`). Detection
 * events are emitted for the memory-basis checks. Optionally, comparisons
 * between successive interior rounds also emit the other-basis checks.
 *
 * Arguments:
 *      `distance`    -- linear size L of the toric lattice (also its code distance),
 *      `rounds`      -- number of syndrome-extraction rounds,
 *      `p`           -- physical error rate,
 *      `is_memory_x` -- true for an X-memory experiment, false for Z-memory.
 *      `include_opposite_basis_detectors` -- include other-basis interior checks.
 * */
// OpenAI GPT-6-Sol: Keep toric SI1000's detector option consistent with the rotated code.
stim::Circuit toric_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x,
                           bool include_opposite_basis_detectors = false,
                           const ComponentErrorDistribution& noise_distribution = {});

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// Triangular, one-logical-qubit color-code memories. Odd distance >= 3,
// rounds >= 1, finite 0 <= p <= 0.2. The memory basis is deterministic at
// initialization/readout. Opposite-basis interior detectors are optional.
// Every flag is a SINGLE-round detector, in both bases, including round one.
// DETECTOR coordinates: (x, y, round, color, basis, is_flag), where
// color = 0/1/2 for red/green/blue; basis = 0/1 for X/Z.
// Bulk faces: cc666 has hexagons; cc488 has red squares and green/blue octagons.
// By default emits only memory-basis check detectors (plus flags in both
// bases). When enabled, also compares opposite-basis checks between rounds.
stim::Circuit cc666_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x,
                            bool include_opposite_basis_detectors = false,
                            const ComponentErrorDistribution& noise_distribution = {});
// Preserve the earlier five-argument call where the fifth argument is noise.
template <typename T>
    requires std::is_same_v<std::remove_cvref_t<T>, ComponentErrorDistribution>
inline stim::Circuit cc666_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x,
                                  T&& noise_distribution) {
    return cc666_si1000(distance,rounds,p,is_memory_x,false,
                        std::forward<T>(noise_distribution));
}
stim::Circuit cc488_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x,
                            bool include_opposite_basis_detectors = false,
                            const ComponentErrorDistribution& noise_distribution = {});
// Preserve the earlier five-argument call where the fifth argument is noise.
template <typename T>
    requires std::is_same_v<std::remove_cvref_t<T>, ComponentErrorDistribution>
inline stim::Circuit cc488_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x,
                                  T&& noise_distribution) {
    return cc488_si1000(distance,rounds,p,is_memory_x,false,
                        std::forward<T>(noise_distribution));
}

// Contributor: OpenAI Codex. Ideal-check color-code memories. Independent
// X(p) and Z(p) data errors precede each extraction round. Preparation,
// check gates and final data readout are perfect; phenomenological check
// outcomes flip with probability measurement_p. Neither model has flags.
// Code capacity always has one extraction round; phenomenological rounds>=1.
// Detector metadata matches SI1000: (x,y,time,color,basis,is_flag=0).
// Other-basis detectors compare interior rounds only when requested.
stim::Circuit cc666_code_capacity(uint32_t distance, double p, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {});
stim::Circuit cc488_code_capacity(uint32_t distance, double p, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {});
stim::Circuit cc666_phenomenological(uint32_t distance, uint32_t rounds,
    double data_p, double measurement_p, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {});
stim::Circuit cc488_phenomenological(uint32_t distance, uint32_t rounds,
    double data_p, double measurement_p, bool is_memory_x,
    bool include_opposite_basis_detectors = false,
    const ComponentErrorDistribution& noise_distribution = {});

#endif // CIRCUIT_GENERATOR_h
