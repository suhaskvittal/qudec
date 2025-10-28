/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "decoder/epr_pym.h"


/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

DECODER_RESULT
EPR_PYMATCHING::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm)
{
    DECODER_RESULT result;

    // split syndrome into inner and outer parts:
    std::vector<GRAPH_COMPONENT_ID> inner_dets,
                                    outer_dets;
    inner_dets.reserve(dets.size());
    outer_dets.reserve(dets.size());

    for (auto d : dets)
    {
        auto it = m_global_to_inner.find(d);
        if (it == m_global_to_inner.end())
            outer_dets.push_back(m_global_to_outer.at(d));
        else
            inner_dets.push_back(it->second);
    }

    auto remaining_inner_dets = decode_inner(inner_dets, result.flipped_observables, debug_strm);
    
    // add remaining detectors to `outer_dets` and decode:
    std::transform(remaining_inner_dets.begin(), remaining_inner_dets.end(), std::back_inserter(outer_dets),
                    [this] (GRAPH_COMPONENT_ID d) { return this->m_inner_to_outer.at(d); });

    [[ maybe_unused ]] pm::total_weight_int w{0};
    decode_outer(outer_dets, result.flipped_observables, debug_strm);

    return result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

std::vector<GRAPH_COMPONENT_ID>
EPR_PYMATCHING::decode_inner(std::vector<GRAPH_COMPONENT_ID> dets, syndrome_ref obs, std::ostream& debug_strm)
{
    std::vector<GRAPH_COMPONENT_ID> remaining_dets;
    for (size_t r = 0; !dets.empty(); r += inner_commit_size)
    {
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void
EPR_PYMATCHING::decode_outer(std::vector<GRAPH_COMPONENT_ID> dets, syndrome_ref obs, std::ostream& debug_strm)
{
    std::vector<uint64_t> detection_events(dets.size());
    std::copy(dets.begin(), dets.end(), detection_events.begin());

    pm::decode_detection_events(outer_mwpm, detection_events, obs.u8, w, false);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

