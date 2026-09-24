/*
 *  author: Suhas Vittal
 *  date:   22 June 2026
 * */

#ifndef DECODER_SURFACE_CODE_MatchingData_h
#define DECODER_SURFACE_CODE_MatchingData_h

#include "decoder/types.h"

#include <cstdint>
#include <iosfwd>
#include <vector>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * `MatchingData` contains information that is useful for debugging
 * matching decoders.
 * */
struct MatchingData
{
    using det_id_type = int64_t;
    using weight_type = uint64_t;
    struct assignment_type
    {
        det_id_type d1;
        det_id_type d2;
        double      pr;
        weight_type w_qu;
        ObsType     frame_flips;

        int matching_step{0};
        int cluster_id{0};
    };

    weight_type                  total_weight{0};
    double                       probability{1.0};
    ObsType                      frame_flips{1};
    std::vector<assignment_type> assignments;

    double gap;

    void add(const assignment_type&);
    void merge(const MatchingData&);
};

void matching_show_diff(std::ostream& ostrm, MatchingData, MatchingData, size_t num_observables);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec

#endif // DECODER_SURFACE_CODE_MatchingData_h
