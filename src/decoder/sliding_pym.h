/*
 *  author: Suhas Vittal
 * */

#ifndef DECODER_SLIDING_PYM_h
#define DECODER_SLIDING_PYM_h

#include "decoder/common.h"
#include "decoder/surface_code.h"

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

class SLIDING_PYMATCHING
{
public:
    const size_t commit_size;
    const size_t window_size;
    const size_t detectors_per_round;
    const size_t total_rounds;

    const bool pym_has_singleton_detectors_in_first_round;
private:
    std::unique_ptr<PYMATCHING> pym;

    void apply_compressed_edge_observables(const pm::CompressedEdge&, const std::unordered_set<GRAPH_COMPONENT_ID>&, DECODER_RESULT&, std::ostream&) const;
public:
    SLIDING_PYMATCHING(const stim::Circuit&, 
                            size_t commit_size, 
                            size_t window_size, 
                            size_t detectors_per_round,
                            size_t total_rounds,
                            bool pym_has_single_detector_first_round);

    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>, std::ostream& debug_strm) const;
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // DECODER_SLIDING_PYM_h
