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

    // it is easier to work with the bit representation:
    syndrome_type syndrome(detectors_per_round * (total_rounds+1));
    syndrome.clear();
    
    for (auto d : dets)
        syndrome[d] = 1;

    size_t r{0};
    while (syndrome.popcnt() > 0)
    {
        if (GL_DEBUG_DECODER)
            debug_strm << "round " << r << ":\n";

        const GRAPH_COMPONENT_ID min_detector = r*detectors_per_round,
                                 max_detector = (r+window_size)*detectors_per_round,
                                 max_commit_detector = (r+commit_size)*detectors_per_round;
        decode_window(syndrome, result.flipped_observables, min_detector, max_detector, max_commit_detector, debug_strm);

        r += commit_size;
    }

    return result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void
SLIDING_PYMATCHING::decode_window(syndrome_ref syndrome, 
                                syndrome_ref obs,
                                GRAPH_COMPONENT_ID min, 
                                GRAPH_COMPONENT_ID window_max,
                                GRAPH_COMPONENT_ID commit_max,
                                std::ostream& debug_strm)
{
    std::vector<uint64_t> window_dets;

    const size_t offset = (min == 0) ? 0 : detectors_per_round;
    for (size_t i = min; i < window_max && i < syndrome.num_bits_padded(); i++)
    {
        if (syndrome[i])
            window_dets.push_back(i - min + offset);
    }

    if (window_dets.empty() || window_dets.front() >= commit_max)
        return;

    if (GL_DEBUG_DECODER)
    {
        debug_strm << "\t(min = " << min << ", max = " << window_max << ", commit_max = " << commit_max 
            << ") detectors in window:";
        for (auto d : window_dets)
            debug_strm << " " << (d-offset+min);
        debug_strm << "\n";
    }
    
    // run pymatching and get matched edges
    std::vector<int64_t> edges;
    pm::decode_detection_events_to_edges(mwpm, window_dets, edges);

    for (size_t i = 0; i < edges.size(); i += 2) 
    {
        const int64_t node1 = edges[i];
        const int64_t node2 = edges[i+1];

        const int64_t true_node1 = (node1 <= 0) ? node1 : node1 - offset + min,
                      true_node2 = (node2 <= 0) ? node2 : node2 - offset + min;

        // Only commit observables if at least one detector is in commit region
        const bool node1_in_commit = (true_node1 >= 0) && (true_node1 < commit_max);
        const bool node2_in_commit = (true_node2 >= 0) && (true_node2 < commit_max);

        if (!node1_in_commit && !node2_in_commit)
        {
            if (GL_DEBUG_DECODER)
            {
                debug_strm << "\tskipping edge between " 
                    << true_node1 << " and " << true_node2 
                    << " (both outside commit region)\n";
            }
            continue;  // Skip edges entirely outside commit region
        }

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
            debug_strm << "\tedge between " << true_node1 << " and " << true_node2 
                << ", flipped observables:";
            for (const size_t obs_idx : obs_indices)
                debug_strm << " " << obs_idx;
            debug_strm << ", weight = " << detector_node.neighbor_weights[neighbor_idx] << "\n";
        }

        syndrome[true_node1] ^= 1;
        if (true_node2 >= 0)
            syndrome[true_node2] ^= 1;
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
