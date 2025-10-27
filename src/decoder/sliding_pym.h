/*
 *  author: Suhas Vittal
 * */

#ifndef DECODER_SLIDING_PYM_h
#define DECODER_SLIDING_PYM_h

#include "decoder/common.h"
#include "decoder/surface_code.h"

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/search/search_detector_node.h"

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

/*
 * IMPORTANT:
 *
 * The provided stim::Circuit provides an error model for a subset of rounds. This error model
 * should have `window_size+1` rounds. This is because:
 *  (1) The first round of the first window is unaffected by CNOT/measurement errors from a prior round (as there
 *      is no previous round). So, the window will span detectors from the first round onward.
 *  (2) Remaining windows will span detectors starting from the second round onward.
 * */

class SLIDING_PYMATCHING
{
public:
    const size_t commit_size;
    const size_t window_size;
    const size_t detectors_per_round;
    const size_t total_rounds;
private:
    pm::Mwpm mwpm;
public:
    SLIDING_PYMATCHING(const stim::Circuit&, 
                            size_t commit_size, 
                            size_t window_size, 
                            size_t detectors_per_round,
                            size_t total_rounds);

    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>, std::ostream& debug_strm);
private:
    GRAPH_COMPONENT_ID min_detector_id(size_t round) const;
    GRAPH_COMPONENT_ID max_detector_id(size_t round) const;
    GRAPH_COMPONENT_ID max_commit_detector_id(size_t round) const;
    GRAPH_COMPONENT_ID detector_offset(size_t round) const;
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

namespace pm_ext  // "pymatching extension"
{

// for decoding syndromes with a sliding window decoder: only detection events in the commit region are committed
// to `syndrome_ref`
//
// returns a set of detector IDs that were committed -- these should not be passed in future calls
std::unordered_set<GRAPH_COMPONENT_ID> 
    decode_detection_events_in_commit_region(pm::Mwpm&,
                                                std::vector<uint64_t> detectors,
                                                uint64_t commit_region_end_id,
                                                syndrome_ref obs,
                                                std::ostream& debug_strm);

} // namespace pm_ext

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // DECODER_SLIDING_PYM_h
