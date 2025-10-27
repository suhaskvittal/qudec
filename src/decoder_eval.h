/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#ifndef DECODER_EVAL_h
#define DECODER_EVAL_h

#include "decoder/common.h"

#include <stim/circuit/circuit.h>

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

// Global debug configuration variable
extern bool GL_DEBUG_DECODER;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

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
 *
 * `GL_DEBUG_DECODER` can be set to `true` to enable debugging logical errors. We also
 * provide an `ERROR_CALLBACK` that can be used to provide more information (i.e., checking
 * the result against a reference decoder).
 *
 * `ERROR_CALLBACK` will be given 
 *      (1) the detectors, 
 *      (2) the true observable flips, 
 *      (3) the decoder's prediction,
 *      (4) a debug stream
 *
 *  `ERROR_CALLBACK` will be called before printing out the debug information for the decoder.
 *  If `ERROR_CALLBACK` returns false, then the debug information will not be printed. This is
 *  useful for limiting debug information for only correctable errors, for example.
 *
 *  Note that `ERROR_CALLBACK` is only called if `GL_DEBUG_DECODER` is true.
 * */

template <class IMPL, class ERROR_CALLBACK> 
void decode(IMPL&,
            DECODER_STATS&,
            syndrome_type detector_flips,
            syndrome_type observable_flips,
            const ERROR_CALLBACK&,
            bool do_not_clock=false);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class IMPL> DECODER_STATS benchmark_decoder(const stim::Circuit&,
                                                        IMPL&,
                                                        uint64_t num_trials,
                                                        uint64_t batch_size=8192,
                                                        bool do_not_clock=false,
                                                        uint64_t seed=0,
                                                        uint64_t stop_at_k_errors=10);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

/*
 * This is an expanded version of `benchmark_decoder`
 *
 * `ERROR_CALLBACK` should take in the detector and observable flips. It need not return anything.
 * */ 

template <class IMPL, class ERROR_CALLBACK> 
DECODER_STATS benchmark_decoder(const stim::Circuit&,
                                IMPL&,
                                uint64_t num_trials,
                                const ERROR_CALLBACK&,
                                uint64_t batch_size=8192,
                                bool do_not_clock=false,
                                uint64_t seed=0,
                                uint64_t stop_at_k_errors=10);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#include "decoder_eval.tpp"

#endif
