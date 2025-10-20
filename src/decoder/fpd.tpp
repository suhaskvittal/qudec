/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "graph/distance.h"

#define TEMPL_PARAMS    template <class LL_DECODER>
#define TEMPL_CLASS     FPD<LL_DECODER>

extern SC_DECODING_GRAPH* create_sc_decoding_graph_from_circuit(const stim::Circuit&);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

TEMPL_PARAMS
TEMPL_CLASS::FPD(const stim::Circuit& circuit, LL_DECODER&& _ll_decoder, FPD_CONFIG conf)
    :conf(conf),
    ll_decoder(std::move(_ll_decoder))
{
    auto* dg = create_sc_decoding_graph_from_circuit(circuit);
    boundary_index = dg->get_vertices().size()-1;

    ec_cache.resize(dg->get_vertices().size());

    std::cout << "FPD: creating cache...\n";

    // build `ec_cache`
    for (GRAPH_COMPONENT_ID id = 0; id < dg->get_vertices().size(); id++)
        init_ec_cache_array(dg, id);

    std::cout << "\tdone\n";

    delete dg;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

TEMPL_PARAMS DECODER_RESULT
TEMPL_CLASS::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm) const
{
    if (dets.size() & 1)
        dets.push_back(boundary_index);

    // identify preferences:
    std::vector<preference_entry_type> prefs = compute_prefs(dets);
    syndrome_type flipped_obs{DEFAULT_OBS_BIT_WIDTH};

    // only operate on detectors with single preferences:
    for (size_t i = 0; i < dets.size(); i++)
    {
        const auto& pi = prefs[i];
        if (dets[i] < 0 || pi.count > 1 || pi.pref == NO_PREF)
            continue;

        // check that other detector also has a single preference and prefers `i` (consensual pairing)
        size_t j = pi.pref;
        const auto& pj = prefs[j];
        if (pj.count > 1 || pj.pref != i)
            continue;

        syndrome_ref obs = ec_cache.at(dets[i]).at(dets[j]).flipped_obs;

#if defined(DEBUG_DECODER)
        debug_strm << "FPD: pairing " << dets[i] << " and " << dets[j] << ", flipped observables:";
        for (size_t i = 0; i < obs.num_bits_padded(); i++)
        {
            if (obs[i])
                debug_strm << " " << i;
        }
        debug_strm << "\n";
#endif

        flipped_obs ^= obs;
        dets[i] = -1;
        dets[j] = -1;
    }

    std::vector<GRAPH_COMPONENT_ID> unmatched_dets;
    unmatched_dets.reserve(dets.size());
    std::copy_if(dets.begin(), dets.end(), std::back_inserter(unmatched_dets), 
                    [b=boundary_index] (auto id) { return id >= 0 && id != b; });

#if defined(DEBUG_DECODER)
    debug_strm << "FPD: unmatched detectors:";
    for (auto d : unmatched_dets)
        debug_strm << " " << d;
    debug_strm << "\n";
#endif

    if (unmatched_dets.empty())
    {
        DECODER_RESULT result;
        result.flipped_observables = std::move(flipped_obs);
        return result;
    }

    auto result = ll_decoder.decode(unmatched_dets, debug_strm);
    result.flipped_observables ^= flipped_obs;
    return result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

TEMPL_PARAMS void 
TEMPL_CLASS::init_ec_cache_array(SC_DECODING_GRAPH* dg, GRAPH_COMPONENT_ID base_id)
{
    using weight_type = int16_t;
    constexpr auto DIJKSTRA_WF = [] (const auto* e) { return e->data.quantized_weight; };

    auto& m = ec_cache[base_id];

    // execute dijkstra's on `dg`
    auto result = graph::dijkstra<weight_type>(*dg, base_id, DIJKSTRA_WF);
    for (GRAPH_COMPONENT_ID id = 0; id < dg->get_vertices().size(); id++)
    {
        if (id == base_id)
            continue;

        auto id_path = graph::dijkstra_path(result.prev, base_id, id, true);

        int8_t len = static_cast<int8_t>(id_path.size());
        if (len > conf.cache_chain_limit)
            continue;

        std::vector<SC_DECODING_GRAPH::VERTEX*> vertex_path(id_path.size());
        std::transform(id_path.begin(), id_path.end(), vertex_path.begin(),
                        [dg] (GRAPH_COMPONENT_ID id) { return dg->get_vertex(id); });

        syndrome_type flipped_obs{DEFAULT_OBS_BIT_WIDTH};
        for (auto it = vertex_path.begin(); it != vertex_path.end()-1; it++)
        {
            auto* e = dg->get_edge_and_fail_if_nonunique(it, it+2);
            for (auto obs_id : e->data.flipped_observables)
                flipped_obs[obs_id] ^= 1;
        }

        error_info e{len, std::move(flipped_obs)};
        m.insert({id, e});
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

TEMPL_PARAMS std::vector<typename TEMPL_CLASS::preference_entry_type>
TEMPL_CLASS::compute_prefs(const std::vector<GRAPH_COMPONENT_ID>& dets) const
{
    std::vector<preference_entry_type> prefs(dets.size());
    // identify the detector each detector ideally wants to pair up with:
    for (size_t i = 0; i < dets.size(); i++)
    {
        const auto& m = ec_cache[dets[i]];
        for (size_t j = 0; j < dets.size(); j++)
        {
            auto it = m.find(dets[j]);
            if (it == m.end())
                continue;

            if (it->second.length < prefs[i].length)
            {
                prefs[i].pref = j;
                prefs[i].length = it->second.length;
            }
        }

        if (prefs[i].pref != NO_PREF)
            prefs[prefs[i].pref].count++;
    }

    return prefs;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
