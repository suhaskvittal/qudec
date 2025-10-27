/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#include "decoder/surface_code.h"
#include "decoder/sliding_pym.h"
#include "graph/distance.h"

#include <iostream>
#include <mutex>
#include <thread>

#include <PerfectMatching.h>

extern bool GL_DEBUG_DECODER;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

SC_DECODING_GRAPH*
create_sc_decoding_graph_from_circuit(const stim::Circuit& circuit)
{
    constexpr stim::DemOptions SC_DEM_OPTS
    {
        true,  // decompose_errors
        true,  // flatten_loops
        false, // allow_gauge_detectors
        0.0,   // approximate_disjoint_errors_threshold
        false, // ignore_decomposition_failures
        false  // block_decomposition_from_introducing_remnant_edges
    };

    stim::DetectorErrorModel dem = stim::circuit_to_dem(circuit, SC_DEM_OPTS);

    if (search_for_bad_dem_errors(dem, circuit))
        throw std::runtime_error("SC_DECODING_GRAPH: found bad DEM errors");

    auto* dg = read_surface_code_decoding_graph(dem);
    quantize_all_edge_weights(dg);
    return dg;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

/*
 * Blossom-V Implementation:
 * */

BLOSSOM5::BLOSSOM5(const stim::Circuit& circuit)
    :dg{create_sc_decoding_graph_from_circuit(circuit)},
    boundary_id(dg->get_vertices().size()-1)
{
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

constexpr auto DIJKSTRA_WF = [] (const auto* e) { return e->data.quantized_weight; };

DECODER_RESULT
BLOSSOM5::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm) const
{
    // add boundary if odd number of dets:
    if (dets.size() & 1)
        dets.push_back(boundary_id);

    // initialize b5 object
    const size_t n = dets.size();
    const size_t m = (n*(n-1)) >> 1;
    b5::PerfectMatching pm(n, m);
    pm.options.verbose = false;

    std::vector<graph::DIJKSTRA_RESULT<weight_type>> dijkstra_results(n);

    for (size_t i = 0; i < n; i++)
    {
        auto et_begin = dets.begin()+i,
             et_end = dets.end();
        auto result = graph::dijkstra<weight_type>(*dg, dets[i], DIJKSTRA_WF, true, et_begin, et_end);
        for (size_t j = i+1; j < n; j++)
        {
            pm.AddEdge(i, j, result.dist[dets[j]]);

#if defined(DEBUG_DECODER)
            debug_strm << "added edge between " << dets[i] << " and " << dets[j] 
                        << " with weight " << result.dist[dets[j]] << "\n";
#endif
        }

        dijkstra_results[i] = std::move(result);
    }

    pm.Solve(); 

    // determine frame changes -- it is faster to just count the parity
    // of observable flips rather than modifying a set over and over again
    DECODER_RESULT result;
    for (size_t i = 0; i < n; i++)
    {
        size_t j = pm.GetMatch(i);
        if (j < i)  // avoid double counting
            continue;

        auto dijk_result = std::move(dijkstra_results[i]);

        GRAPH_COMPONENT_ID src_id = dets[i],
                           dst_id = dets[j];

        auto id_path = graph::dijkstra_path(dijk_result.prev, src_id, dst_id, true);
        std::vector<SC_DECODING_GRAPH::VERTEX*> vertex_path(id_path.size());
        std::transform(id_path.begin(), id_path.end(), vertex_path.begin(),
                        [this] (GRAPH_COMPONENT_ID id) { return dg->get_vertex(id); });

        [[ maybe_unused ]] std::unordered_map<int64_t, size_t> path_flips;
        for (auto it = vertex_path.begin(); it != vertex_path.end()-1; it++)
        {
            auto* e = dg->get_edge_and_fail_if_nonunique(it, it+2);
            for (auto obs_id : e->data.flipped_observables)
            {
                result.flipped_observables[obs_id] ^= 1;
                path_flips[obs_id] ^= 1;
            }
        }

#if defined (DEBUG_DECODER)
        debug_strm << "match between " << src_id << " and " << dst_id << ", flipped observables:";
        for (const auto& [x, flips] : path_flips)
            if (flips & 1)
                debug_strm << " " << x;
        debug_strm << "\n";
#endif
    }

    return result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

/*
 * PyMatching Implementation:
 * */

PYMATCHING::PYMATCHING(const stim::Circuit& circuit)
    :mwpm{pymatching_create_mwpm_from_circuit(circuit, GL_DEBUG_DECODER)},
    num_observables{circuit.count_observables()}
{}

DECODER_RESULT
PYMATCHING::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm)
{
    // Convert detector IDs to PyMatching format (uint64_t vector)
    std::vector<uint64_t> detection_events(dets.size());
    std::copy(dets.begin(), dets.end(), detection_events.begin());

    // Create observables array
    DECODER_RESULT result;

    // Perform matching using PyMatching's decode function
    if (GL_DEBUG_DECODER)
    {
        debug_strm << "pymatching verbose (not performant):\n";
        pm_ext::decode_detection_events_in_commit_region(mwpm, 
                                                        detection_events,
                                                        std::numeric_limits<uint64_t>::max(),
                                                        result.flipped_observables,
                                                        debug_strm);
    }
    else
    {
        pm::total_weight_int weight{0};
        pm::decode_detection_events(mwpm, detection_events, result.flipped_observables.u8, weight, false);
    }

    return result;
}

pm::Mwpm 
pymatching_create_mwpm_from_circuit(const stim::Circuit& circuit, bool enable_search_flooder)
{
    auto dem = stim::circuit_to_dem(circuit, {true, true, false, 0.0, false, false});
    return pm::detector_error_model_to_mwpm(dem, pm::NUM_DISTINCT_WEIGHTS, enable_search_flooder, false);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

