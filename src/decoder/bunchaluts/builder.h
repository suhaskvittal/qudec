/*
 *  author: Suhas Vittal
 *  date:   21 September 2026
 * */

#ifndef DECODER_BAL_BUILDER_h
#define DECODER_BAL_BUILDER_h

#include "decoder/bunchaluts/hypergraph.h"

namespace decoder
{
namespace bal
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

struct LUTKey
{
};

struct LUTEntry
{
    using frame_repr_type = std::unordered_set<size_t>; 

    /*
     * Set of frame flips for common case (most detectors)
     * and unique case (some detectors)
     * */
    struct
    {
        frame_repr_type common;
        std::unordered_map<id_type, frame_repr_type> unique;
    } frame_flips;

    double probability;
};

using LUT = std::unordered_map<LUTKey, LUTEntry>;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

enum class DemType { toric_xz, surface_xz };

void build_lut(DemType, 
                const stim::DetectorErrorModel&, 
                size_t distance, 
                double target_error_probability);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace bal
} // namespace decoder

#endif
