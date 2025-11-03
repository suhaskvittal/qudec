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
    struct decode_options
    {
        // do not commit any boundary edges for detectors in this set:
        std::unordered_set<GRAPH_COMPONENT_ID> do_not_commit_boundary_edges_set{};

        bool do_not_commit_any_boundary_edges{false};
    };

    using window_bounds_type = std::tuple<GRAPH_COMPONENT_ID, GRAPH_COMPONENT_ID, GRAPH_COMPONENT_ID>;

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

    // generic decode function:
    void decode_and_update_inplace(syndrome_ref, syndrome_ref, std::ostream& debug_strm, decode_options);
private:
    void decode_window(syndrome_ref, syndrome_ref, window_bounds_type, std::ostream&, decode_options); 
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // DECODER_SLIDING_PYM_h
