/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#ifndef DECODER_COMMON_h
#define DECODER_COMMON_h

#include "globals.h"

#include <stim.h>

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

using syndrome_type = stim::simd_bits<64>;
using syndrome_ref = stim::simd_bits_range_ref<64>;
using obs_type = syndrome_type;
using obs_ref = syndrome_ref;

/*
 * `result_type` is a generic result output that must
 * be returned by a decoder. Feel free to add extra
 * fields if they are useful.
 * */
struct result_type
{
    obs_type flipped_obs{1};

    // Decoder specific parameters:
    int matching_weight{};
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder

#endif // DECODER_COMMON_h
