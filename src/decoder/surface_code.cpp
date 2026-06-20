/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#include "decoder/surface_code.h"

#include <pymatching/sparse_blossom/driver/user_graph.h>

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PYMATCHING::PYMATCHING(const stim::DetectorErrorModel& dem)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    mwpm_(pm::detector_error_model_to_mwpm(dem, pm::NUM_DISTINCT_WEIGHTS))
{}

result_type
PYMATCHING::decode(syndrome_ref syn)
{
    // Collect indices of fired detectors.
    std::vector<uint64_t> det_events;
    for (size_t i = 0; i < num_detectors; i++)
        if (syn[i]) 
            det_events.push_back(i);
    // Decode.
    result_type res;
    res.flipped_obs = obs_type(num_observables);
    pm::total_weight_int weight = 0;
    pm::decode_detection_events(mwpm_, det_events, res.flipped_obs.u8, weight, false);
    return res;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder
