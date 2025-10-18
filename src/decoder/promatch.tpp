/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#define TEMPL_PARAMS    template <class LL_DECODER>
#define TEMPL_CLASS     PROMATCH_SW<LL_DECODER>

// in `surface_code.h`
extern SC_DECODING_GRAPH* create_sc_decoding_graph_from_circuit(const stim::Circuit&);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

TEMPL_PARAMS
TEMPL_CLASS::PROMATCH_SW(const stim::Circuit& circuit, LL_DECODER&& ll_decoder)
    :dg(create_sc_decoding_graph_from_circuit(circuit)),
    ll_decoder(std::move(ll_decoder)),
    boundary_index(dg->get_vertices().size()-1)
{
}

TEMPL_PARAMS DECODER_RESULT
TEMPL_CLASS::decode(std::vector<GRAPH_COMPONENT_ID> dets, std::ostream& debug_strm) const
{
    if (dets.size() & 1)
        dets.push_back(boundary_index);

    // first, compute the degree of each detector when only considering the subgraph
    // induced by the set of detectors `dets`:
    std::vector<PROMATCH_INFO> pm(dets.size());
    initialize_induced_subgraph(pm, dets);

    // try and match bits
    syndrome_type pm_obs_flips{DEFAULT_OBS_BIT_WIDTH};
    bool done{false};
    while (!done)
    {
#if defined (DEBUG_DECODER)
        debug_strm << "promatch: step\n";
#endif
        done = !promatch_step(pm, pm_obs_flips, debug_strm);
    }

    // delete matched bits from `dets`:
    std::vector<GRAPH_COMPONENT_ID> unmatched_dets;
    unmatched_dets.reserve(dets.size()/2);
    for (size_t i = 0; i < pm.size(); i++)
    {
        if (pm[i].induced_degree >= 0 && dets[i] != boundary_index)
            unmatched_dets.push_back(dets[i]);
    }

#if defined (DEBUG_DECODER)
    debug_strm << "promatch: unmatched detectors:";
    for (auto d : unmatched_dets)
        debug_strm << " " << d;
    debug_strm << "\n";
#endif
    
    if (unmatched_dets.empty())
    {
        DECODER_RESULT result;
        result.flipped_observables = std::move(pm_obs_flips);
        return result;
    }

    auto main_result = ll_decoder.decode(unmatched_dets, debug_strm);
    main_result.flipped_observables ^= pm_obs_flips;

    return main_result;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

TEMPL_PARAMS void
TEMPL_CLASS::initialize_induced_subgraph(std::vector<PROMATCH_INFO>& pm, const std::vector<GRAPH_COMPONENT_ID>& dets) const
{
    for (size_t i = 0; i < dets.size(); i++)
    {
        pm[i].det_idx = i;
        pm[i].vertex = dg->get_vertex(dets[i]);
        pm[i].neighbors.reserve(2);
    }

    // utility function for adding neighbors in promatch
    auto add_neighbor = [&pm] (size_t i, size_t j)
    {
        pm[i].neighbors.push_back(j);
        pm[j].neighbors.push_back(i);
        pm[i].induced_degree++;
        pm[j].induced_degree++;
    };

    // identify neighbors in the induced subgraph:
    for (size_t i = 0; i < dets.size(); i++)
    {
        auto* v = pm[i].vertex;
        for (size_t j = i+1; j < dets.size(); j++)
        {
            auto* w = pm[j].vertex;
            if (dg->get_edge_and_fail_if_nonunique(v, w) != nullptr)
                add_neighbor(i, j);
        }
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

TEMPL_PARAMS bool
TEMPL_CLASS::promatch_step(std::vector<PROMATCH_INFO>& pm, syndrome_ref flipped_obs, std::ostream& debug_strm) const
{
    bool any_deletion{false};
    for (size_t i = 0; i < pm.size(); i++)
    {
        // do not try anything if the detector is a singleton
        auto& pi = pm[i];
        if (pi.induced_degree <= 0)
            continue;

        for (size_t j : pi.neighbors)
        {
            if (i > j)
                continue;

            auto& pj = pm[j];
            if (pj.induced_degree <= 0)
                continue;

            // we want to check if the deletion of `(i,j)` causes a singleton to appear:
            std::vector<int8_t> updated_degrees(pm.size());
            std::transform(pm.begin(), pm.end(), updated_degrees.begin(), [] (const auto& p) { return p.induced_degree; });
            for (size_t k : pi.neighbors)
                updated_degrees[k]--;
            for (size_t k : pj.neighbors)
                updated_degrees[k]--;

            // `pi` and `pj` are ok to have 0, but the others, not so much -- set their degrees to 1
            // so they do not get matched
            updated_degrees[pi.det_idx] = -1;
            updated_degrees[pj.det_idx] = -1;

            auto it = std::find(updated_degrees.begin(), updated_degrees.end(), 0);
            if (it == updated_degrees.end())
            {
                // update induced degrees:
                for (size_t k : pi.neighbors)
                    pm[k].induced_degree--;
                for (size_t k : pj.neighbors)
                    pm[k].induced_degree--;

                pi.induced_degree = -1;
                pj.induced_degree = -1;

                // update decoder result:
                auto* e = dg->get_edge_and_fail_if_nonunique(pi.vertex, pj.vertex);
                for (auto obs_id : e->data.flipped_observables)
                    flipped_obs[obs_id] ^= 1;

                any_deletion = true;

#if defined (DEBUG_DECODER)
                debug_strm << "\tpromatch: matched " << pi.vertex->id << " and " << pj.vertex->id << ", obs flips:";
                for (auto obs_id : e->data.flipped_observables)
                    debug_strm << " " << obs_id;
                debug_strm << "\n";
#endif

                break;
            }
        }
    }

    return any_deletion;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#undef TEMPL_PARAMS
#undef TEMPL_CLASS
