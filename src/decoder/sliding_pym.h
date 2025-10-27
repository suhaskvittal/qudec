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

class SLIDING_PYMATCHING
{
public:
    const size_t commit_size;
    const size_t window_size;
    const size_t detectors_per_round;
    const size_t total_rounds;

    const bool pym_has_singleton_detectors_in_first_round;
private:
    pm::Mwpm mwpm;
public:
    SLIDING_PYMATCHING(const stim::Circuit&, 
                            size_t commit_size, 
                            size_t window_size, 
                            size_t detectors_per_round,
                            size_t total_rounds,
                            bool pym_has_single_detector_first_round);

    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>, std::ostream& debug_strm) const;
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

namespace pm_ext  // "pymatching extension"
{

// for decoding syndromes with a sliding window decoder: only detection events in the commit region are committed
// to `syndrome_ref`
//
// returns a set of detector IDs that were committed -- these should not be passed in future calls
template <class DET_ITER>
std::unordered_set<GRAPH_COMPONENT_ID> 
    decode_detection_events_in_commit_region(pm::Mwpm&,
                                                DET_ITER d_begin, 
                                                DET_ITER d_end,
                                                const std::unordered_set<GRAPH_COMPONENT_ID>& commit_region,
                                                syndrome_ref obs);

} // namespace pm_ext

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

/*
 * IMPLEMENTATION OF `decode_detection_events_in_commit_region`
 * */ 

namespace pm_ext
{

// Based on original PyMatching code by Oscar Higgott and Craig Gidney
// Modified for sliding window decoding with commit regions
template <class DET_ITER> std::unordered_set<GRAPH_COMPONENT_ID>
decode_detection_events_in_commit_region(pm::Mwpm& mwpm,
                                            DET_ITER d_begin,
                                            DET_ITER d_end,
                                            const std::unordered_set<GRAPH_COMPONENT_ID>& commit_region,
                                            syndrome_ref obs)
{
    // Convert detection events to vector format expected by PyMatching
    std::vector<uint64_t> detection_events(d_begin, d_end);

    // Get edges from the matching solution using the public API
    std::vector<int64_t> edges;
    pm::decode_detection_events_to_edges(mwpm, detection_events, edges);

    // Track which detectors were committed (had their observables applied)
    std::unordered_set<GRAPH_COMPONENT_ID> committed_detectors;

    // Process edges in pairs: [node1, node2, node1, node2, ...]
    for (size_t i = 0; i < edges.size(); i += 2) 
    {
        const int64_t node1 = edges[i];
        const int64_t node2 = edges[i + 1];

        // Only commit observables if at least one detector is in commit region
        const bool node1_in_commit = (node1 >= 0) && commit_region.count(static_cast<GRAPH_COMPONENT_ID>(node1));
        const bool node2_in_commit = (node2 >= 0) && commit_region.count(static_cast<GRAPH_COMPONENT_ID>(node2));

        if (!node1_in_commit && !node2_in_commit)
            continue;  // Skip edges entirely outside commit region

        if (node1 >= 0 && node2 >= 0) 
        {
            // Regular edge between two detectors
            const auto& detector_node = mwpm.search_flooder.graph.nodes[node1];
            const auto* neighbor_ptr = &mwpm.search_flooder.graph.nodes[node2];

            const size_t neighbor_idx = detector_node.index_of_neighbor(neighbor_ptr);
            const auto& obs_indices = detector_node.neighbor_observable_indices[neighbor_idx];

            // Apply observable flips
            for (const size_t obs_idx : obs_indices)
                obs[obs_idx] ^= 1;

            committed_detectors.insert(static_cast<GRAPH_COMPONENT_ID>(node1));
            committed_detectors.insert(static_cast<GRAPH_COMPONENT_ID>(node2));
        } 
        else if (node1 >= 0 && node2 == -1) 
        {
            // Boundary edge (node2 == -1 represents boundary)
            const auto& detector_node = mwpm.search_flooder.graph.nodes[node1];

            const size_t neighbor_idx = detector_node.index_of_neighbor(nullptr);
            const auto& obs_indices = detector_node.neighbor_observable_indices[neighbor_idx];

            // Apply observable flips
            for (const size_t obs_idx : obs_indices)
                obs[obs_idx] ^= 1;

            // Mark detector as committed if it's in commit region
            committed_detectors.insert(static_cast<GRAPH_COMPONENT_ID>(node1));
        }
    }

    return committed_detectors;
}

}   // namespace pm_ext

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // DECODER_SLIDING_PYM_h
