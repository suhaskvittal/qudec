/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#ifndef DECODING_GRAPH_h
#define DECODING_GRAPH_h

#include "hypergraph.h"

#include <stim/dem/detector_error_model.h>

#include <cmath>

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

struct DETECTOR_DATA
{
    // make sure to use `RED = 1` in your stim circuits!!!!
    enum class COLOR { NONE, RED, GREEN, BLUE };

    COLOR              color{COLOR::NONE};
    bool               is_flag{false};
    bool               is_boundary{false};
};

struct DECODER_ERROR_DATA
{
    using quantized_weight_type = int16_t;

    double                      error_probability;
    quantized_weight_type       quantized_weight;
    std::unordered_set<int64_t> flipped_observables;
};

// `DG_TYPE` is a generic decoding graph type
// colorability here refers to the maximum order of the hypergraph,
// as this often determines the colorability of the graph
template <size_t COLORABILITY=2> using DG_TYPE = HYPERGRAPH<DETECTOR_DATA, DECODER_ERROR_DATA, COLORABILITY>;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

using SC_DECODING_GRAPH = DG_TYPE<2>;

SC_DECODING_GRAPH* read_surface_code_decoding_graph(const stim::DetectorErrorModel& dem);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class DG_PTR> void
quantize_all_edge_weights(DG_PTR& dg)
{
    for (auto* e : dg->get_edges())
    {
        double ep = e->data.error_probability;
        DECODER_ERROR_DATA::quantized_weight_type qw = std::round(-std::log(ep)*100);
        e->data.quantized_weight = qw;
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif
