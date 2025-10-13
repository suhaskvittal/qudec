/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#include "decoding_graph.h"
#include "io/dem.h"

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

SC_DECODING_GRAPH*
read_surface_code_decoding_graph(const stim::DetectorErrorModel& dem)
{
    constexpr GRAPH_COMPONENT_ID BOUNDARY_INDEX{-1};

    auto result = io::read_dem_block(dem);

    // add `+1` for the boundary index:
    SC_DECODING_GRAPH* gr = new SC_DECODING_GRAPH{result.detectors.size()+1, result.errors.size()};
    auto* boundary = gr->add_vertex(BOUNDARY_INDEX, {});

    // create vertices for each detector:
    for (const auto& d : result.detectors)
        gr->add_vertex(d, {});
    
    // create edges for each error:
    for (const auto& [dets, ed] : result.errors)
    {
        if (dets.size() > 2 || dets.empty())
        {
            throw std::runtime_error("SC_DECODING_GRAPH: got error with " + std::to_string(dets.size()) + " detectors");
        }
        else
        {
            // if `dets.size() == 2`, then both `boundary`  will be overwritten.
            // otherwise, `boundary` will be overwritten only once.
            std::vector<SC_DECODING_GRAPH::VERTEX*> vlist(2, boundary);
            std::transform(dets.begin(), dets.end(), vlist.begin(),
                            [gr] (int64_t d) { return gr->get_vertex(d); });
            gr->add_edge(vlist, ed);
        }
    }

    return gr;
}


/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
