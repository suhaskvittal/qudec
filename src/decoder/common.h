/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#ifndef DECODER_COMMON_h
#define DECODER_COMMON_h

#include "decoding_graph.h"

#include <stim/mem/simd_bits.h>

#include <unordered_set>

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

using syndrome_type = stim::simd_bits<stim::MAX_BITWORD_WIDTH>;
using syndrome_ref = stim::simd_bits_range_ref<stim::MAX_BITWORD_WIDTH>;

constexpr size_t DEFAULT_OBS_BIT_WIDTH{256};

struct DECODER_RESULT
{
    syndrome_type flipped_observables{DEFAULT_OBS_BIT_WIDTH};
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // DECODER_COMMON_h
