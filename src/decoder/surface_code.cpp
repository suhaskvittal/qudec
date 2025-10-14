/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#include "decoder/surface_code.h"
#include "graph/distance.h"

#include <mutex>
#include <thread>

#include <PerfectMatching.h>

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
    :dg{create_sc_decoding_graph_from_circuit(circuit)}
{}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

constexpr auto DIJKSTRA_WF = [] (const auto* e) { return e->data.quantized_weight; };

DECODER_RESULT
BLOSSOM5::decode(std::vector<GRAPH_COMPONENT_ID>&& dets) const
{
    // initialize b5 object
    const size_t n = dets.size();
    const size_t m = (n*(n-1)) >> 1;
    b5::PerfectMatching pm(n, m);

    std::vector<graph::DIJKSTRA_RESULT<weight_type>> dijkstra_results(n);

    for (size_t i = 0; i < n; i++)
    {
        auto et_begin = dets.begin()+i,
             et_end = dets.end();
        auto result = graph::dijkstra<weight_type>(*dg, dets[i], DIJKSTRA_WF, true, et_begin, et_end);
        for (size_t j = i+1; j < n; j++)
            pm.AddEdge(i, j, result.dist[dets[j]]);

        dijkstra_results[i] = std::move(result);
    }

    pm.Solve(); 

    // determine frame changes -- it is faster to just count the parity
    // of observable flips rather than modifying a set over and over again
    std::unordered_map<int64_t, size_t> observable_flips_by_id{};
    observable_flips_by_id.reserve(4);
    for (size_t i = 0; i < n; i++)
    {
        size_t j = pm.GetMatch(i);
        if (j < i)  // avoid double counting
            continue;

        auto result = std::move(dijkstra_results[i]);

        GRAPH_COMPONENT_ID src_id = dets[i],
                           dst_id = dets[j];

        auto id_path = graph::dijkstra_path(result.prev, src_id, dst_id, true);
        std::vector<SC_DECODING_GRAPH::VERTEX*> vertex_path(id_path.size());
        std::transform(id_path.begin(), id_path.end(), vertex_path.begin(),
                        [this] (GRAPH_COMPONENT_ID id) { return dg->get_vertex(id); });

        for (auto it = vertex_path.begin(); it != vertex_path.end()-1; it++)
        {
            auto* e = dg->get_edge_and_fail_if_nonunique(it, it+2);
            for (auto obs_id : e->data.flipped_observables)
                observable_flips_by_id[obs_id] ^= 1;
        }
    }

    std::unordered_set<int64_t> flipped_observables;
    flipped_observables.reserve(4);
    for (const auto& [x, flips] : observable_flips_by_id)
    {
        if (flips & 1)
            flipped_observables.insert(x);
    }

    return DECODER_RESULT{flipped_observables};
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

