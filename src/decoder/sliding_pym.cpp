/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "decoder/sliding_pym.h"
#include "decoder/surface_code.h"

extern bool GL_DEBUG_DECODER;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

SLIDING_PYMATCHING::SLIDING_PYMATCHING(
        const stim::Circuit& circuit, 
        size_t _commit_size, 
        size_t _window_size,
        size_t _detectors_per_round,
        size_t _total_rounds)
    :commit_size(_commit_size),
    window_size(_window_size),
    detectors_per_round(_detectors_per_round),
    total_rounds(_total_rounds),
    mwpm{pymatching_create_mwpm_from_circuit(circuit, true)}
{
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

DECODER_RESULT
SLIDING_PYMATCHING::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm)
{
    DECODER_RESULT result;
    for (size_t r = 0; !dets.empty(); r += commit_size)
    {
        const GRAPH_COMPONENT_ID min_detector = min_detector_id(r),
                                 max_detector = max_detector_id(r),
                                 max_commit_detector = max_commit_detector_id(r),
                                 det_offset = detector_offset(r);
        const GRAPH_COMPONENT_ID total_offset{min_detector-det_offset};

        if (dets.front() < min_detector)
        {
            std::cerr << "dets.front() = " << dets.front() << ", min_detector = " << min_detector << "\n";
            throw std::runtime_error("SLIDING_PYMATCHING: dets.front() < min_detector");
        }

        if (dets.front() >= max_commit_detector)
            continue;

        auto d_begin = dets.begin();
        auto d_end = std::find_if(d_begin, dets.end(),
                [max_detector] (GRAPH_COMPONENT_ID d) { return d > max_detector; });

        std::vector<uint64_t> pm_dets;
        pm_dets.reserve(std::distance(d_begin, d_end));
        std::transform(d_begin, d_end, std::back_inserter(pm_dets),
                [total_offset] (GRAPH_COMPONENT_ID d) 
                { 
                    return static_cast<uint64_t>(d-total_offset);
                });

        const uint64_t commit_end = static_cast<uint64_t>(max_commit_detector-total_offset);

        if (GL_DEBUG_DECODER)
        {
            debug_strm << "round " << r << ":\n\tdets =";
            for (auto d : pm_dets)
                debug_strm << " " << d << " (" << (d+total_offset) << ")";
            debug_strm << "\ncommit region is until " << commit_end << " (i.e. " << (commit_end+total_offset) << ")";
            debug_strm << "\n";
        }

        auto done = pm_ext::decode_detection_events_in_commit_region(mwpm, 
                                                                    pm_dets,
                                                                    commit_end,
                                                                    result.flipped_observables,
                                                                    debug_strm);

        auto it = std::remove_if(dets.begin(), dets.end(),
                [&done, total_offset] (GRAPH_COMPONENT_ID d) 
                { 
                    return done.count(d-total_offset);
                });
        dets.erase(it, dets.end());
    }

    return result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

GRAPH_COMPONENT_ID
SLIDING_PYMATCHING::min_detector_id(size_t round) const
{
    return static_cast<GRAPH_COMPONENT_ID>(round * detectors_per_round);
}

GRAPH_COMPONENT_ID
SLIDING_PYMATCHING::max_detector_id(size_t round) const
{
    return static_cast<GRAPH_COMPONENT_ID>((round+window_size) * detectors_per_round);
}

GRAPH_COMPONENT_ID
SLIDING_PYMATCHING::max_commit_detector_id(size_t round) const
{
    return static_cast<GRAPH_COMPONENT_ID>((round+commit_size) * detectors_per_round);
}

GRAPH_COMPONENT_ID
SLIDING_PYMATCHING::detector_offset(size_t round) const
{
    return (round>0) ? detectors_per_round : 0;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

namespace pm_ext
{

std::unordered_set<GRAPH_COMPONENT_ID> 
decode_detection_events_in_commit_region(pm::Mwpm& mwpm,
                                        std::vector<uint64_t> detection_events,
                                        uint64_t commit_region_end_id,
                                        syndrome_ref obs,
                                        std::ostream& debug_strm)
{
    // Get edges from the matching solution using the public API
    std::vector<int64_t> edges;
    pm::decode_detection_events_to_edges(mwpm, detection_events, edges);

    // Track which detectors were committed (had their observables applied)
    std::unordered_set<GRAPH_COMPONENT_ID> committed_detectors;

    // if a matching is >1 weight and exits the commit region, we also want to commit it unless we hit the boundary.
    std::unordered_set<int64_t> ok_to_commit_outside;

    // track visited indices so we don't double count:
    std::unordered_set<size_t> visited_idx;

    // Process edges in pairs: [node1, node2, node1, node2, ...] -- make multiple passes if necessary
    size_t pass{0};
    bool any_commits;
    do
    {
        any_commits = false;

        debug_strm << "\tpass " << pass << ":\n";
        for (size_t i = 0; i < edges.size(); i += 2) 
        {
            if (visited_idx.count(i))
                continue;

            const int64_t node1 = edges[i];
            const int64_t node2 = edges[i + 1];

            // Only commit observables if at least one detector is in commit region
            const bool node1_in_commit = (node1 >= 0) && (node1 < commit_region_end_id || ok_to_commit_outside.count(node1));
            const bool node2_in_commit = (node2 >= 0) && (node2 < commit_region_end_id || ok_to_commit_outside.count(node2));

            if (!node1_in_commit && !node2_in_commit)
            {
                if (GL_DEBUG_DECODER)
                {
                    debug_strm << "\t\tskipping edge between " 
                        << node1 << " and " << node2 
                        << " (both outside commit region)\n";
                }
                continue;  // Skip edges entirely outside commit region
            }

            visited_idx.insert(i);

            // Regular edge between two detectors
            const auto& detector_node = mwpm.search_flooder.graph.nodes[node1];
            auto* neighbor_ptr = node2 >= 0 ? &mwpm.search_flooder.graph.nodes[node2] : nullptr;

            const size_t neighbor_idx = detector_node.index_of_neighbor(neighbor_ptr);
            const auto& obs_indices = detector_node.neighbor_observable_indices[neighbor_idx];

            // Apply observable flips
            for (const size_t obs_idx : obs_indices)
                obs[obs_idx] ^= 1;

            if (GL_DEBUG_DECODER)
            {
                debug_strm << "\t\tedge between " << node1 << " and " << node2 << ", flipped observables:";
                for (const size_t obs_idx : obs_indices)
                    debug_strm << " " << obs_idx;
                debug_strm << ", weight = " << detector_node.neighbor_weights[neighbor_idx] << "\n";
            }

            committed_detectors.insert(static_cast<GRAPH_COMPONENT_ID>(node1));
            ok_to_commit_outside.insert(static_cast<GRAPH_COMPONENT_ID>(node1));

            if (node2 >= 0)
            {
                committed_detectors.insert(static_cast<GRAPH_COMPONENT_ID>(node2));
                ok_to_commit_outside.insert(static_cast<GRAPH_COMPONENT_ID>(node2));
            }

            any_commits = true;
        }
        pass++;
    }
    while (any_commits);

    return committed_detectors;
}

}   // namespace pm_ext

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
