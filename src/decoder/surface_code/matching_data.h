/*
 *  author: Suhas Vittal
 *  date:   22 June 2026
 * */

#ifndef DECODER_SURFACE_CODE_MATCHING_DATA_h
#define DECODER_SURFACE_CODE_MATCHING_DATA_h

#include "decoder/types.h"

#include <cstdint>
#include <iosfwd>
#include <vector>

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * `MATCHING_DATA` contains information that is useful for debugging
 * matching decoders.
 * */
struct MATCHING_DATA
{
    using det_id_type = int64_t;
    using weight_type = uint64_t;
    struct assignment_type
    {
        det_id_type d1;
        det_id_type d2;
        double      pr;
        weight_type w_qu;
        obs_type    frame_flips;

        int matching_step{0};
        int cluster_id{0};
    };

    weight_type                  total_weight{0};
    double                       probability{1.0};
    obs_type                     frame_flips{1};
    std::vector<assignment_type> assignments;

    void add(const assignment_type&);
    void merge(const MATCHING_DATA&);
};

void matching_show_diff(std::ostream& ostrm, MATCHING_DATA, MATCHING_DATA, size_t num_observables);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder

#endif // DECODER_SURFACE_CODE_MATCHING_DATA_h
