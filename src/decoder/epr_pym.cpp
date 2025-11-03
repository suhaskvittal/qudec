/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "decoder/epr_pym.h"

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

EPR_PYMATCHING::EPR_PYMATCHING(const stim::Circuit& global,
                                const stim::Circuit& inner, 
                                const stim::Circuit& outer,
                                size_t code_distance,
                                size_t _num_super_rounds,
                                size_t _num_sub_rounds_per_super_round)
    :global_circuit(global),
    inner_circuit(inner),
    outer_circuit(outer),
    num_super_rounds(_num_super_rounds),
    num_sub_rounds_per_super_round(_num_sub_rounds_per_super_round)
{
    // initialize detector map:

    read_first_round_of_detectors(outer,
            [this] (auto d, auto base, auto super_round, auto sub_round) 
            {
                if (this->m_detector_info.count(base))
                {
                    throw std::runtime_error("EPR_PYMATCHING: duplicate base detector: " + std::to_string(base));
                }
                this->m_detector_info[base].outer_id = d;
                this->outer_detectors_per_round++;
            });

    read_first_round_of_detectors(inner,
            [this] (auto d, auto base, auto super_round, auto sub_round) 
            {
                if (!this->m_detector_info.count(base))
                {
                    throw std::runtime_error("EPR_PYMATCHING: no detector in outer circuit for base: " 
                                                + std::to_string(base));
                }
                this->m_detector_info[base].inner_id = d;
                this->inner_detectors_per_round++;
            });

    // initialize remaining fields:
    total_detectors_per_super_round = inner_detectors_per_round*num_sub_rounds_per_super_round
                                        + outer_detectors_per_round;

    dec_inner = std::make_unique<SLIDING_PYMATCHING>(inner, 
                                                    code_distance,
                                                    2*code_distance,
                                                    inner_detectors_per_round,
                                                    num_sub_rounds_per_super_round+1);
    dec_outer = std::make_unique<PYMATCHING>(outer);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

DECODER_RESULT
EPR_PYMATCHING::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm)
{
    DECODER_RESULT result;

    syndrome_type s_outer(outer_detectors_per_round * (num_super_rounds+1));
    s_outer.clear();

    // decode the sub rounds of each super round in a first pass:
    // note that there are `num_super_rounds-1` epochs of sub-rounds:
    for (size_t sr = 0; sr < num_super_rounds-1; sr++)
    {
        syndrome_type s_inner(inner_detectors_per_round * (num_sub_rounds_per_super_round+1));
        s_inner.clear();

        size_t d_min = total_detectors_per_super_round * sr,
               d_max = total_detectors_per_super_round * (sr+1);

        auto d_begin = std::find_if(dets.begin(), dets.end(), [d_min] (auto d) { return d >= d_min; });
        auto d_end = std::find_if(d_begin, dets.end(), [d_max] (auto d) { return d >= d_max; });

        std::for_each(d_begin, d_end,
                [this, &s_inner] (auto d)
                {
                    auto idx = this->get_inner_syndrome_detector_idx(d);
                    if (idx.has_value())
                        s_inner[*idx] ^= 1;
                });
        
        SLIDING_PYMATCHING::decode_options opts;
        opts.do_not_commit_boundary_edges_set = do_not_commit_boundary_edges_set;
        dec_inner->decode_and_update_inplace(s_inner, result.flipped_observables, debug_strm, opts);

        // move remaining bits from inner to outer:
        size_t nonzero_bits = s_inner.popcount();
        for (auto it = d_begin; it != d_end && nonzero_bits > 0; ++it)
        {
            auto idx = get_inner_syndrome_detector_idx(*it);
            if (idx.has_value() && s_inner[*idx])
            {
                nonzero_bits--;

                // update `s_outer`:
                size_t outer_idx = m_detector_info.at(*it).outer_id + outer_detectors_per_round*sr;
                s_outer[outer_idx] ^= 1;
            }
        }
    }

    // convert to detector list:
    std::vector<GRAPH_COMPONENT_ID> outer_dets;
    for (size_t i = 0; i < s_outer.num_bits_padded(); i++)
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

std::optional<size_t>
EPR_PYMATCHING::get_inner_syndrome_detector_idx(size_t global_detector_idx)
{
    const auto coords = this->global_circuit.coords_of_detector(global_detector_idx);
    size_t base = static_cast<size_t>(coords[BASE_DETECTOR_IDX]),
           sub_round_idx = static_cast<size_t>(coords[SUB_ROUND_IDX]);

    if (m_detector_info.at(base).inner_id < 0)
        return std::nullopt;

    size_t idx = m_detector_info.at(base).inner_id + inner_detectors_per_round*sub_round_idx;
    return std::make_optional(idx);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

