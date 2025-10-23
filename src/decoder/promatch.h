/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#ifndef DECODER_PROMATCH_h
#define DECODER_PROMATCH_h

#include "decoding_graph.h"
#include "decoder/common.h"

// Global debug configuration variable
extern bool GL_DEBUG_DECODER;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

struct PROMATCH_INFO
{
    using vertex_ptr = SC_DECODING_GRAPH::VERTEX*;

    size_t               det_idx;
    vertex_ptr           vertex;
    std::vector<size_t>  neighbors;
    int8_t               induced_degree{0};
};

/*
 * Implements a software implementation of Promatch (pre-decoding)
 * atop some lower-level decoder.
 * */

template <class LL_DECODER>
class PROMATCH_SW 
{
public:
private:
    std::unique_ptr<SC_DECODING_GRAPH> dg;
    LL_DECODER ll_decoder;

    const GRAPH_COMPONENT_ID boundary_index;
public:
    PROMATCH_SW(const stim::Circuit&, LL_DECODER&&);
    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>, std::ostream& debug_strm) const;
private:
    void initialize_induced_subgraph(std::vector<PROMATCH_INFO>&, const std::vector<GRAPH_COMPONENT_ID>&) const;
    bool promatch_step(std::vector<PROMATCH_INFO>&, syndrome_ref, std::ostream& debug_strm) const;
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#include "promatch.tpp"

#endif  // DECODER_PROMATCH_h
