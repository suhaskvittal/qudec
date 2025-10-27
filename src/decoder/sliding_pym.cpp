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
    for (size_t r = 0; r < total_rounds && !dets.empty(); r += commit_size)
    {
        const GRAPH_COMPONENT_ID max_detector{(r + window_size) * detectors_per_round},
                                 max_commit_detector{(r + commit_size) * detectors_per_round};

        if (dets.front() > max_commit_detector)
            continue;

        auto d_begin = dets.begin();
        auto d_end = std::find_if(d_begin, dets.end(),
                [max_detector] (GRAPH_COMPONENT_ID d) { return d > max_detector; });
        auto d_commit_end = std::find_if(d_begin, d_end,
                [max_commit_detector] (GRAPH_COMPONENT_ID d) { return d > max_commit_detector; });

        std::unordered_set<GRAPH_COMPONENT_ID> commit_region(d_begin, d_commit_end);
        auto done = pm_ext::decode_detection_events_in_commit_region(mwpm, 
                                                                    d_begin, 
                                                                    d_end, 
                                                                    commit_region,
                                                                    result.flipped_observables);

        dets.erase(d_begin, d_commit_end);  // can bulk delete these
        auto it = std::remove_if(dets.begin(), dets.end(),
                [&done] (GRAPH_COMPONENT_ID d) { return done.count(d); });
        dets.erase(it, dets.end());
    }

    return result;
}


/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
