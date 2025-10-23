/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#include "decoding_graph.h"
#include "io/dem.h"

#include <stim/simulators/error_matcher.h>

#include <iostream>

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

bool
search_for_bad_dem_errors(const stim::DetectorErrorModel& dem, const stim::Circuit& circuit)
{
    bool found_bad_errors{false};

    stim::DetectorErrorModel bad_errors_dem;

    dem.iter_flatten_error_instructions([&](const stim::DemInstruction& error_inst)
    {
        bool has_detectors{false};
        bool has_observables{false};

        for (const auto& target : error_inst.target_data)
        {
            if (target.is_relative_detector_id())
                has_detectors = true;
            else if (target.is_observable_id())
                has_observables = true;
        }

        if (!has_detectors && has_observables)
        {
            found_bad_errors = true;

            bad_errors_dem.append_error_instruction(
                error_inst.arg_data[0],
                error_inst.target_data,
                ""
            );
        }
    });

    if (found_bad_errors)
    {
        std::cerr << "Found errors that only flip observables (no detectors):\n";

        auto explained_errors = stim::ErrorMatcher::explain_errors_from_circuit(
            circuit,
            &bad_errors_dem,
            true
        );

        for (const auto& explained_error : explained_errors)
            std::cerr << explained_error << "\n";
    }

    return found_bad_errors;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
