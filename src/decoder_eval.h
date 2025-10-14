/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#ifndef DECODER_EVAL_h
#define DECODER_EVAL_h

#include "decoder/common.h"

#include "stim/circuit/circuit.h"
#include "stim/mem/simd_bits.h"

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

using syndrome_type = stim::simd_bits<stim::MAX_BITWORD_WIDTH>;

struct DECODER_STATS
{
    using hw_histogram_type = std::array<uint64_t, 128>;

    uint64_t errors{0};
    uint64_t trials{0};
    uint64_t trivial_trials{0};
    uint64_t total_time_us{0};

    hw_histogram_type time_us_by_hamming_weight{};
    hw_histogram_type trials_by_hamming_weight{};
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

/*
 * Calls `IMPL::decode` which should take in a `std::vector<GRAPH_COMPONENT_ID>` of
 * detector indices that are flipped, and return a `DECODER_RESULT`
 *
 * This function manages any updates to `DECODER_STATS` during the call
 *
 * As clocks do introduce a substantial overhead to runtime (not during
 * `IMPL::decode` but moreso around it), `do_not_clock` can be set to `true`
 * to disable the timing of the call.
 * */
template <class IMPL> void decode(IMPL&,
                                    DECODER_STATS&,
                                    syndrome_type detector_flips,
                                    syndrome_type observable_flips,
                                    bool do_not_clock=false);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class IMPL> DECODER_STATS benchmark_decoder(const stim::Circuit&,
                                                        IMPL&,
                                                        uint64_t num_trials,
                                                        uint64_t batch_size=8192,
                                                        bool do_not_clock=false,
                                                        uint64_t seed=0);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#include "decoder_eval.tpp"

#endif
