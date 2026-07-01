/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#ifndef DECODER_COMMON_h
#define DECODER_COMMON_h

#include "decoder/surface_code/matching_data.h"
#include "globals.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace decoder
{

/*
 * `result_type` is a generic result output that must
 * be returned by a decoder. Feel free to add extra
 * fields if they are useful.
 * */
struct result_type
{
    obs_type flipped_obs{1};

    // Decoder specific parameters:
    MATCHING_DATA matching_data;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder

#endif // DECODER_COMMON_h
