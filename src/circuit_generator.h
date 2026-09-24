/*
 *  author: Claude Sonnet 5.1
 *  date:   26 July 2026
 * */

#ifndef CIRCUIT_GENERATOR_h
#define CIRCUIT_GENERATOR_h

#include <stim.h>

#include <cstdint>

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
                        bool include_opposite_basis_detectors = false);

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
                           bool include_opposite_basis_detectors = false);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#endif // CIRCUIT_GENERATOR_h
