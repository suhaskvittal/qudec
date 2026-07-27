/*
 *  author: Suhas Vittal
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
 * Detection events are only emitted for the checks relevant to the requested
 * memory experiment: Z checks (which detect X errors) for a Z-memory
 * experiment, and X checks for an X-memory experiment.
 *
 * Arguments:
 *      `distance`    -- code distance of the rotated surface code,
 *      `rounds`      -- number of syndrome-extraction rounds,
 *      `p`           -- physical error rate,
 *      `is_memory_x` -- true for an X-memory experiment, false for Z-memory.
 * */
stim::Circuit si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#endif // CIRCUIT_GENERATOR_h
