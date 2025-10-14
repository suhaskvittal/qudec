/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#ifndef DECODER_SURFACE_CODE_h
#define DECODER_SURFACE_CODE_h

#include "decoding_graph.h"
#include "decoder/common.h"

#include <stim/circuit/circuit.h>
#include <stim/dem/detector_error_model.h>
#include <stim/util_top/circuit_to_dem.h>

/*
 * This file contains all decoders for the surface code.
 * The simple-to-implement decoders are implemented in `surface_code.cpp`
 *
 * If you write your own decoder, but want your own header/source file,
 * include the header here and add the source file to CMake.
 * */


/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

/*
 * This is a utility function for creating a 
 * decoding graph from a stim circuit.
 * */ 

SC_DECODING_GRAPH* read_surface_code_decoding_graph(const stim::Circuit&);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

struct BLOSSOM5
{
public:
    using weight_type = DECODER_ERROR_DATA::quantized_weight_type;
private:
    std::unique_ptr<SC_DECODING_GRAPH> dg;
public:
    BLOSSOM5(const stim::Circuit&);
    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>&&) const;
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // DECODER_SURFACE_CODE_h

