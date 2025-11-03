/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#ifndef DECODER_EPR_PYM_h
#define DECODER_EPR_PYM_h

#include "decoder/sliding_pym.h"
#include "decoder/surface_code.h"

#include <memory>
#include <optional>

extern bool GL_EPR_PYMATCHING_VERBOSE;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

/*
 * Decoder for handling lattice surgery merges between two
 * hardware substrates with significantly different gate latencies.
 *
 * A decoder for the inner circuit is used to decode errors that
 * are contained within the faster substrate.
 *
 * A decoder for the outer circuit is used to decoder errors on
 * the slower substrate, or any errors that cross substrates.
 *
 * Any detectors that the first decoder maps to the boundary are
 * retried by the second decoder.
 *
 * The first decoder is a sliding window decoder. The second is not.
 * */

class EPR_PYMATCHING
{
public:
    constexpr static size_t BASE_DETECTOR_IDX{2};
    constexpr static size_t SUPER_ROUND_IDX{3};
    constexpr static size_t SUB_ROUND_IDX{4};

    struct detector_info
    {
        GRAPH_COMPONENT_ID inner_id{-1};
        GRAPH_COMPONENT_ID outer_id;
    };

    using detector_info_map_type = std::unordered_map<GRAPH_COMPONENT_ID, detector_info>;

    const stim::Circuit& global_circuit;
    const stim::Circuit& inner_circuit;
    const stim::Circuit& outer_circuit;
private:
    size_t num_super_rounds;
    size_t num_sub_rounds_per_super_round;
    
    size_t outer_detectors_per_round{};
    size_t inner_detectors_per_round{};
    size_t total_detectors_per_super_round;

    std::unique_ptr<SLIDING_PYMATCHING> dec_inner;
    std::unique_ptr<PYMATCHING>         dec_outer;

    detector_info_map_type m_detector_info;

    std::unordered_set<GRAPH_COMPONENT_ID> do_not_commit_boundary_edges_set;
public:
    EPR_PYMATCHING(const stim::Circuit& global,
                    const stim::Circuit& inner, 
                    const stim::Circuit& outer,
                    size_t code_distance,
                    size_t num_super_rounds,
                    size_t num_sub_rounds_per_super_round);

    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>, std::ostream& debug_strm);
private:
    std::optional<size_t> get_inner_syndrome_detector_idx(size_t global_detector_idx);
    size_t get_outer_syndrome_detector_idx(size_t global_detector_idx);
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void concat_debug_strm(std::ostream& target, std::stringstream& source, size_t tab_count);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

// helper function for initializing `m_detector_info`:
template <class UPDATE_CALLBACK> void
read_first_round_of_detectors(const stim::Circuit& circ, const UPDATE_CALLBACK& cb)
{
    // only need first round of syndromes from each circuit:
    for (size_t i = 0; i < circ.count_detectors(); i++)
    {
        const auto coords = circ.coords_of_detector(i);
        size_t overall_round = static_cast<size_t>(coords[1]),
               base = static_cast<size_t>(coords[EPR_PYMATCHING::BASE_DETECTOR_IDX]),
               super_round_idx = static_cast<size_t>(coords[EPR_PYMATCHING::SUPER_ROUND_IDX]),
               sub_round_idx = static_cast<size_t>(coords[EPR_PYMATCHING::SUB_ROUND_IDX]);
        
        if (overall_round > 0)
            continue;

        cb(i, base, super_round_idx, sub_round_idx);
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // DECODER_EPR_PYM_h
