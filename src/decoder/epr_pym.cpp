/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "decoder/epr_pym.h"

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void
_init_detector_map_from_circuit(EPR_PYMATCHING::detector_map& m, const stim::Circuit& circuit)
{
    for (size_t i = 0; i < circuit.count_detectors(); i++)
    {
        const auto coords = circuit.coords_of_detector(i);
        if (EPR_PYMATCHING::GLOBAL_DETECTOR_COORD_IDX < coords.size())
        {
            const auto global_id = static_cast<GRAPH_COMPONENT_ID>(coords[EPR_PYMATCHING::GLOBAL_DETECTOR_COORD_IDX]);
            m[global_id] = static_cast<GRAPH_COMPONENT_ID>(i);
        }
    }
}

EPR_PYMATCHING::EPR_PYMATCHING(const stim::Circuit& inner, 
                                const stim::Circuit& outer,
                                size_t inner_commit_size,
                                size_t inner_window_size,
                                size_t inner_detectors_per_round,
                                size_t inner_total_rounds)
    :inner_circuit(inner),
    outer_circuit(outer),
    dec_inner{new SLIDING_PYMATCHING(inner, inner_commit_size, inner_window_size, inner_detectors_per_round, inner_total_rounds)},
    dec_outer{new PYMATCHING(outer)}
{
    // initialize detector maps:

    // Get detector counts for both circuits
    const size_t inner_detector_count = inner_circuit.count_detectors();
    const size_t outer_detector_count = outer_circuit.count_detectors();

    // Build global_to_inner map
    _init_detector_map_from_circuit(m_global_to_inner, inner_circuit);
    _init_detector_map_from_circuit(m_global_to_outer, outer_circuit);

    // Build inner_to_outer map by matching global detector IDs
    for (const auto& [global_id, inner_id] : m_global_to_inner) 
    {
        const auto outer_it = m_global_to_outer.find(global_id);
        if (outer_it != m_global_to_outer.end()) 
        {
            m_inner_to_outer[inner_id] = outer_it->second;

            // only add inner detectors that correspond to some outer detector to this set:
            do_not_commit_boundary_edges_set.insert(inner_id);
        }
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

DECODER_RESULT
EPR_PYMATCHING::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm)
{
    DECODER_RESULT result;

    // split syndrome into inner and outer parts:
    syndrome_type inner_s(inner_circuit.count_detectors()),
                  outer_s(outer_circuit.count_detectors());

    inner_s.clear();
    outer_s.clear();

    for (auto d : dets)
    {
        const auto inner_it = m_global_to_inner.find(d);
        if (inner_it != m_global_to_inner.end())
            inner_s[inner_it->second] ^= 1;
        else
            outer_s[m_global_to_outer.at(d)] ^= 1;
    }

    // decode inner part first:
    SLIDING_PYMATCHING::decode_options opts;
    opts.do_not_commit_boundary_edges_set = do_not_commit_boundary_edges_set;
    dec_inner->decode_and_update_inplace(inner_s, result.flipped_observables, debug_strm);

    // move r=maining syndrome bits from inner to outer (also make detector list)
    for (size_t i = 0; i < inner_s.num_bits_padded(); i++)
    {
        if (inner_s[i])
            outer_s[m_inner_to_outer.at(i)] ^= 1;
    }

    // convert to detector list:
    std::vector<GRAPH_COMPONENT_ID> outer_dets;
    for (size_t i = 0; i < outer_s.num_bits_padded(); i++)
    {
        if (outer_s[i])
            outer_dets.push_back(i);
    }

    // decode outer part:
    auto outer_result = dec_outer->decode(outer_dets, debug_strm);
    result.flipped_observables ^= outer_result.flipped_observables;

    return result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

