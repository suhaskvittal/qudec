/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#ifndef DECODER_EPR_PYM_h
#define DECODER_EPR_PYM_h

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
    constexpr size_t GLOBAL_DETECTOR_COORD_IDX{2};

    using detector_map = std::unordered_map<GRAPH_COMPONENT_ID, GRAPH_COMPONENT_ID>;

    const stim::Circuit& inner_circuit;
    const stim::Circuit& outer_circuit;

    const size_t inner_detectors_per_round;
    const size_t inner_commit_size;
    const size_t inner_window_size;
    const size_t inner_total_rounds;
private:
    std::unique_ptr<SLIDING_PYMATCHING> d_inner;
    std::unique_ptr<PYMATCHING>         d_outer;

    detector_map m_inner_to_outer;
    detector_map m_global_to_inner;
    detector_map m_global_to_outer;
public:
    EPR_PYMATCHING(const stim::Circuit& inner, 
                    const stim::Circuit& outer,
                    size_t inner_detectors_per_round,
                    size_t inner_commit_size,
                    size_t inner_window_size);

    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>, std::ostream& debug_strm);
private:
    std::vector<GRAPH_COMPONENT_ID> decode_inner(std::vector<GRAPH_COMPONENT_ID>, syndrome_ref, std::ostream&); 
    void                            decode_outer(std::vector<GRAPH_COMPONENT_ID>, syndrome_ref, std::ostream&);
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // DECODER_EPR_PYM_h
