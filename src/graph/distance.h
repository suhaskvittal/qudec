/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#ifndef GRAPH_DISTANCE_h
#define GRAPH_DISTANCE_h

#include "hypergraph.h"

#include <vector>

namespace graph
{

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class WEIGHT_TYPE>
struct DIJKSTRA_RESULT
{
    std::vector<WEIGHT_TYPE>        dist;
    std::vector<GRAPH_COMPONENT_ID> prev;
};

/*
 * Precondition: `dijkstra` assumes that the vertex id's are contiguous and start from 0.
 * */

template <class WEIGHT_TYPE,
            class GRAPH_TYPE, 
            class WEIGHT_FUNCTION,
            class EARLY_TERM_ITER=std::vector<GRAPH_COMPONENT_ID>::const_iterator> 
DIJKSTRA_RESULT<WEIGHT_TYPE>
dijkstra(const GRAPH_TYPE&,
            GRAPH_COMPONENT_ID src,
            const WEIGHT_FUNCTION&,
            bool terminate_early=false,
            EARLY_TERM_ITER et_begin={},
            EARLY_TERM_ITER et_end={});

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

std::vector<GRAPH_COMPONENT_ID> dijkstra_path(const std::vector<GRAPH_COMPONENT_ID>& prev, 
                                                GRAPH_COMPONENT_ID src,
                                                GRAPH_COMPONENT_ID dst,
                                                bool reverse_ok=false);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

}  // namespace graph

#include "distance.tpp"

#endif
