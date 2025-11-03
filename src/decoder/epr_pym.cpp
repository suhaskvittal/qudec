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
    size_t inner_commit_size = code_distance;
    size_t inner_window_size = 2*inner_commit_size;
    size_t inner_total_rounds = (num_sub_rounds_per_super_round+1) * num_super_rounds + 1;

    dec_inner = std::make_unique<SLIDING_PYMATCHING>(inner, 
                                                    inner_commit_size,
                                                    inner_window_size,
                                                    inner_detectors_per_round,
                                                    inner_total_rounds);
    dec_outer = std::make_unique<PYMATCHING>(outer);

    std::cout << "EPR_PYMATCHING: initialized with " 
        << "inner_detectors_per_round = " << inner_detectors_per_round
        << ", outer_detectors_per_round = " << outer_detectors_per_round
        << ", total_detectors_per_super_round = " << total_detectors_per_super_round
        << ", inner decoder total rounds = " << inner_total_rounds
        << ", global total detectors = " << global_circuit.count_detectors()
        << "\n";
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

DECODER_RESULT
EPR_PYMATCHING::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm)
{
    DECODER_RESULT result;

    syndrome_type s_outer(outer_circuit.count_detectors());
    s_outer.clear();
        
    // initialize decoder optionsL:
    SLIDING_PYMATCHING::decode_options opts;
    opts.do_not_commit_any_boundary_edges = true;
    opts.do_not_commit_boundary_edges_set = do_not_commit_boundary_edges_set;

    if (GL_EPR_PYMATCHING_VERBOSE)
    {
        std::cout << "EPR_PYMATCHING: decode start... dets =";
        for (auto d : dets)
            std::cout << " " << d;
        std::cout << "\n";
    }

    std::stringstream inner_debug_strm, outer_debug_strm;

    // decode the sub rounds of each super round in a first pass:
    const size_t num_inner_bits = inner_detectors_per_round * ((num_sub_rounds_per_super_round+1)*num_super_rounds + 1);
    syndrome_type s_inner(num_inner_bits);
    s_inner.clear();

    if (GL_DEBUG_DECODER)
        debug_strm << "inner syndrome detectors (bit count = " << num_inner_bits << ") =";

    for (auto d : dets)
    {
        auto idx = get_inner_syndrome_detector_idx(d);
        if (idx.has_value())
        {
            s_inner[*idx] ^= 1;

            if (GL_DEBUG_DECODER)
                debug_strm << " " << d << "(" << *idx << ")";
        }
        else
        {
            s_outer[get_outer_syndrome_detector_idx(d)] ^= 1;
        }
    }

    if (GL_DEBUG_DECODER)
        debug_strm << "\ninner decoder call:\n";

    dec_inner->decode_and_update_inplace(s_inner, result.flipped_observables, inner_debug_strm, opts);
    concat_debug_strm(debug_strm, inner_debug_strm, 1);

    // move remaining bits from inner to outer:
    size_t nonzero_bits = s_inner.popcnt();
    for (size_t i = 0; i < dets.size() && nonzero_bits > 0; i++)
    {
        auto d = dets[i];
        auto idx = get_inner_syndrome_detector_idx(d);
        if (idx.has_value() && s_inner[*idx])
        {
            nonzero_bits--;

            // update `s_outer`:
            size_t outer_idx = get_outer_syndrome_detector_idx(d);
            s_outer[outer_idx] ^= 1;

            if (GL_DEBUG_DECODER)
            {
                debug_strm << "\tmoving bit " << d  << "(" << outer_idx << ")"
                    << " from inner to outer (parity = " << s_outer[outer_idx] << ")\n";
            }
        }
    }

    // convert to detector list:
    std::vector<GRAPH_COMPONENT_ID> outer_dets;
    for (size_t i = 0; i < s_outer.num_bits_padded(); i++)
    {
        if (s_outer[i])
            outer_dets.push_back(i);
    }

    if (GL_DEBUG_DECODER)
        debug_strm << "outer decoder call:\n";

    // decode outer part:
    auto outer_result = dec_outer->decode(outer_dets, outer_debug_strm);
    result.flipped_observables ^= outer_result.flipped_observables;

    concat_debug_strm(debug_strm, outer_debug_strm, 1);

    return result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

std::optional<size_t>
EPR_PYMATCHING::get_inner_syndrome_detector_idx(size_t global_detector_idx)
{
    const auto coords = this->global_circuit.coords_of_detector(global_detector_idx);
    size_t base = static_cast<size_t>(coords[BASE_DETECTOR_IDX]),
           super_round_idx = static_cast<size_t>(coords[SUPER_ROUND_IDX]),
           sub_round_idx = static_cast<size_t>(coords[SUB_ROUND_IDX]);

    if (!m_detector_info.count(base))
    {
        throw std::runtime_error("EPR_PYMATCHING: no detector in outer circuit for base: " 
                                + std::to_string(base));
    }

    if (m_detector_info.at(base).inner_id < 0)
        return std::nullopt;

    size_t overall_sub_round_idx = super_round_idx*(num_sub_rounds_per_super_round+1) + sub_round_idx;
    size_t idx = m_detector_info.at(base).inner_id + inner_detectors_per_round*overall_sub_round_idx;

    /*
    std::cout << "global_detector_idx = " << global_detector_idx
                << ", base = " << base
                << ", super_round_idx = " << super_round_idx
                << ", sub_round_idx = " << sub_round_idx
                << ", overall_sub_round_idx = " << overall_sub_round_idx
                << ", idx = " << idx 
                << "\n";
    */

    return std::make_optional(idx);
}

size_t
EPR_PYMATCHING::get_outer_syndrome_detector_idx(size_t global_detector_idx)
{
    const auto coords = this->global_circuit.coords_of_detector(global_detector_idx);
    size_t base = static_cast<size_t>(coords[BASE_DETECTOR_IDX]),
           super_round_idx = static_cast<size_t>(coords[SUPER_ROUND_IDX]);

    if (!m_detector_info.count(base))
    {
        throw std::runtime_error("EPR_PYMATCHING: no detector in outer circuit for base: " 
                                + std::to_string(base));
    }

    size_t idx = m_detector_info.at(base).outer_id + outer_detectors_per_round*super_round_idx;
    return idx;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void
concat_debug_strm(std::ostream& target, std::stringstream& source, size_t tab_count)
{
    std::string line;
    while (std::getline(source, line))
        target << std::string(tab_count, '\t') << line << "\n";
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

