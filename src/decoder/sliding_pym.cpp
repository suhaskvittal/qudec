/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "decoder/sliding_pym.h"

extern bool GL_DEBUG_DECODER;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

SLIDING_PYMATCHING::SLIDING_PYMATCHING(
        const stim::Circuit& circuit, 
        size_t _commit_size, 
        size_t _window_size,
        size_t _detectors_per_round,
        size_t _total_rounds,
        bool _pym_has_singleton_detectors_in_first_round)
    :commit_size(_commit_size),
    window_size(_window_size),
    detectors_per_round(_detectors_per_round),
    total_rounds(_total_rounds),
    pym_has_singleton_detectors_in_first_round(_pym_has_singleton_detectors_in_first_round),
    pym{std::make_unique<PYMATCHING>(circuit)}
{
    if (total_rounds % commit_size != 0 || total_rounds % window_size != 0)
        throw std::runtime_error("SLIDING_PYMATCHING: total_rounds must be a multiple of commit_size and window_size");
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

DECODER_RESULT
SLIDING_PYMATCHING::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm) const
{
    DECODER_RESULT result;
    auto d_begin = dets.begin();

    for (size_t r = 0; r < total_rounds; r += commit_size)
    {
        const GRAPH_COMPONENT_ID max_detector{(r + window_size) * detectors_per_round},
                                 max_commit_detector{(r + commit_size) * detectors_per_round};

        if (d_begin == dets.end() || *d_begin > max_detector)
            continue;

        auto d_end = std::find_if(d_begin, dets.end(),
                [max_detector] (GRAPH_COMPONENT_ID d) { return d > max_detector; });

        std::vector<GRAPH_COMPONENT_ID> window_dets(d_begin, d_end);

        auto d_commit_end = std::find_if(d_begin, d_end,
                [max_commit_detector] (GRAPH_COMPONENT_ID d) { return d > max_commit_detector; });
        std::unordered_set<GRAPH_COMPONENT_ID> commit_region(d_begin, d_commit_end);

        auto compressed_edges = pym->decode_and_get_compressed_edges(window_dets, debug_strm);

        for (const auto& edge : compressed_edges)
            apply_compressed_edge_observables(edge, commit_region, result, debug_strm);

        d_begin = d_end;
    }

    return result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void
SLIDING_PYMATCHING::apply_compressed_edge_observables(const pm::CompressedEdge& edge,
                                                      const std::unordered_set<GRAPH_COMPONENT_ID>& commit_region,
                                                      DECODER_RESULT& result,
                                                      std::ostream& debug_strm) const
{
    // Check if this edge has at least one endpoint in commit region
    GRAPH_COMPONENT_ID det1 = static_cast<GRAPH_COMPONENT_ID>(edge.detector_node1->id_in_graph),
                       det2 = static_cast<GRAPH_COMPONENT_ID>(edge.detector_node2->id_in_graph);

    if (commit_region.count(det1) || commit_region.count(det2))
    {
        // Extract observables from compressed edge and apply them
        pm::total_weight_int weight = 0;
        std::vector<uint8_t> observables(pym->get_user_graph().get_num_observables(), 0);

        const_cast<pm::Mwpm&>(pym->mwpm).extract_paths_from_match_edges({edge}, observables.data(), weight);

        for (size_t i = 0; i < observables.size(); i++)
        {
            if (observables[i] & 1)
                result.flipped_observables.insert(static_cast<int64_t>(i));
        }

        if (GL_DEBUG_DECODER)
            debug_strm << "committing compressed edge (" << det1 << ", " << det2 << ")\n";
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
