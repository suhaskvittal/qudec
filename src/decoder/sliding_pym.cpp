/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "decoder/sliding_pym.h"
#include "decoder/surface_code.h"

#include <utility>

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

    decode_and_update_inplace(syndrome, result.flipped_observables, debug_strm, decode_options{});
    return result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void
SLIDING_PYMATCHING::decode_and_update_inplace(syndrome_ref syndrome, 
                                                syndrome_ref obs, 
                                                std::ostream& debug_strm, 
                                                decode_options opts)
{
    size_t r{0};
    while (r < total_rounds+1 && syndrome.popcnt() > 0)
    {
        if (GL_DEBUG_DECODER)
            debug_strm << "round " << r << ":\n";

        const GRAPH_COMPONENT_ID min_detector = r*detectors_per_round,
                                 max_detector = (r+window_size)*detectors_per_round,
                                 max_commit_detector = (r+commit_size)*detectors_per_round;

        window_bounds_type bounds{min_detector, max_detector, max_commit_detector};
        decode_window(syndrome, obs, bounds, debug_strm, opts);

        r += commit_size;
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void
SLIDING_PYMATCHING::decode_window(syndrome_ref syndrome, 
                                syndrome_ref obs,
                                window_bounds_type bounds,
                                std::ostream& debug_strm,
                                decode_options opts)
{
    std::vector<uint64_t> window_dets;

    auto [d_min, d_max, d_commit_max] = bounds;
    const size_t offset = (d_min == 0) ? 0 : detectors_per_round;

    // use `_true_id` to make the code less verbose
    auto _true_id = [d_min, offset] (const int64_t node) { return (node <= 0) ? node : node - offset + d_min; };

    for (size_t i = d_min; i < d_max && i < syndrome.num_bits_padded(); i++)
    {
        if (syndrome[i])
            window_dets.push_back(i-d_min+offset);
    }

    if (window_dets.empty() || _true_id(window_dets.front()) >= d_commit_max)
        return;

    if (GL_DEBUG_DECODER)
    {
        debug_strm << "\t(min = " << d_min << ", max = " << d_max << ", commit_max = " << d_commit_max 
            << ") detectors in window:";
        for (auto d : window_dets)
            debug_strm << " " << _true_id(d);
        debug_strm << "\n";
    }
    
    // run pymatching and get matched edges
    std::vector<int64_t> edges;
    pm::decode_detection_events_to_edges(mwpm, window_dets, edges);

    for (size_t i = 0; i < edges.size(); i += 2) 
    {
        int64_t node1 = edges[i];
        int64_t node2 = edges[i+1];

        if (node1 < 0)
            std::swap(node1, node2);

        const int64_t true_node1 = (node1 < 0) ? node1 : _true_id(node1),
                      true_node2 = (node2 < 0) ? node2 : _true_id(node2);

        // Only commit observables if at least one detector is in commit region
        const bool node1_in_commit = (true_node1 >= 0) && (true_node1 < d_commit_max);
        const bool node2_in_commit = (true_node2 >= 0) && (true_node2 < d_commit_max);

        if (!node1_in_commit && !node2_in_commit)
        {
            if (GL_DEBUG_DECODER)
            {
                debug_strm << "\tskipping edge between " << true_node1 << " and " << true_node2 
                    << " (both outside commit region)\n";
            }
            continue;  // Skip edges entirely outside commit region
        }

        if (opts.do_not_commit_boundary_edges_set.count(true_node1))
        {
            if (GL_DEBUG_DECODER)
            {
                debug_strm << "\tskipping edge between " << true_node1 << " and " << true_node2 
                    << " (touches boundary)\n";
            }
            continue;  // Skip edges touching boundary
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
