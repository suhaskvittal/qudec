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
    auto result = io::read_dem_block(dem);

    // add `+1` for the boundary index:
    SC_DECODING_GRAPH* gr = new SC_DECODING_GRAPH{result.detectors.size()+1, result.errors.size()};

    // create vertices for each detector:
    for (const auto& [id, dd] : result.detectors)
        gr->add_vertex(static_cast<GRAPH_COMPONENT_ID>(id), dd);

    // add boundary:
    auto* boundary = gr->add_vertex(result.detectors.size(), {});
    boundary->data.is_boundary = true;
    
    // create edges for each error:
    for (const auto& [dets, ed] : result.errors)
    {
        if (dets.size() > 2 || dets.empty())
        {
            std::cerr << "error info:"
                    << "\n\terror prob = " << ed.error_probability
                    << "\n\tdetectors =";
            for (auto d : dets)
                std::cerr << " " << d;

            std::cerr << "\n";
                
            throw std::runtime_error("SC_DECODING_GRAPH: got error with " + std::to_string(dets.size()) + " detectors");
        }
        else
        {
            // if `dets.size() == 2`, then both `boundary`  will be overwritten.
            // otherwise, `boundary` will be overwritten only once.
            std::array<SC_DECODING_GRAPH::VERTEX*, 2> vlist{boundary, boundary};
            std::transform(dets.begin(), dets.end(), vlist.begin(),
                            [gr] (int64_t d) { return gr->get_vertex(d); });
            // first, check if `vlist` already has an edge:
            auto* e = gr->get_edge_and_fail_if_nonunique(vlist.begin(), vlist.end());
            if (e != nullptr)
            {
                // update edge probability:
                double p1 = e->data.error_probability,
                       p2 = ed.error_probability;
                e->data.error_probability = p1*(1-p2) + (1-p1)*p2;
            }
            gr->add_edge(vlist.begin(), vlist.end(), ed);
        }
    }

    return gr;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
