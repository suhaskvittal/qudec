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

namespace dec
{

/*
 * `result_type` is a generic result output that must
 * be returned by a decoder. Feel free to add extra
 * fields if they are useful.
 * */
struct result_type
{
    ObsType flipped_obs{1};

    // Decoder specific parameters:
    MatchingData matching_data;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec

#endif // DECODER_COMMON_h
