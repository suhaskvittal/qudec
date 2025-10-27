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

#include <pymatching/sparse_blossom/driver/user_graph.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <pymatching/sparse_blossom/matcher/mwpm.h>
#include <pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.h>

#include <iosfwd>

// Global debug configuration variable
extern bool GL_DEBUG_DECODER;

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

    GRAPH_COMPONENT_ID boundary_id;
public:
    BLOSSOM5(const stim::Circuit&);
    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>, std::ostream& debug_strm) const;
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

struct PYMATCHING
{
public:
    using weight_type = double;
private:
    pm::Mwpm mwpm;
    const size_t num_observables;
public:
    using compressed_edge_result = std::vector<uint64_t>;

    PYMATCHING(const stim::Circuit&);
    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>, std::ostream& debug_strm);
};

pm::Mwpm pymatching_create_mwpm_from_circuit(const stim::Circuit&, bool enable_search_flooder=false);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

// Other decoders:
#include "decoder/sliding_pym.h"

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // DECODER_SURFACE_CODE_h

