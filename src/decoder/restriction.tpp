/*
 *  author: Suhas Vittal
 *  date:   16 October 2026
 * */

#include <cassert>
#include <iostream>

#define TEMPL_PARAMS <class M>
#define TEMPL_CLASS RestrictionDecoder<M>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS Result
TEMPL_CLASS::decode(SyndromeRef s, ObsRef true_obs) const
{
    detectors_.clear();
    orbits_.clear();
    // populate `detectors_`
    detectors_.reserve(s.popcnt());
    orbits_.reserve(s.popcnt());
    for (size_t i = 0; i < num_detectors; i++)
    {
        if (s[i])
        {
            detectors_.push_back( static_cast<det_id_type>(i) );
            if (hg_->v(i).color == restriction::C_RED)
                orbits_[i] = lift_orbit{.root=i};
        }
    }

    // initial `orbits_` for each red detector.

    // run matching procedures:
    matching(rg_matcher_, C_BLUE);
    matching(rb_matcher_, C_GREEN);

    Result out;
    out.flipped_obs = ObsType(num_observables);
    for (const auto& orb : orbits_)
        out.flipped_obs ^= lift(orb).flipped_obs;
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS void
TEMPL_CLASS::matching(ma_ptr& m, Color ignored_color) const
{
    
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS void
TEMPL_CLASS::build_connected_components(std::vector<restriction::Match> matchings)
{
    std::vector<restriction::ConnectedComponent> components;

    restriction::dfs::graph_type gr;

    std::unordered_map<det_id_type, hg::id_type> vertex_id_map;
    vertex_id_map.reserve(detectors_.size());
    for (auto d : detectors_)
        vertex_id_map[d] = gr.add_vertex(restriction::dfs::v_data_type{ d }); 
    for (Color c : {C_RED, C_GREEN, C_BLUE})
        vertex_id_map[boundary_id(c)] = gr.add_vertex(restriction::dfs::v_data_type{boundary_id(c)});

    // build graph from `m`
    for (const auto& m : matchings)
    {
        for (auto d : {m.x, m.y})
        {
            if (vertex_id_map.count(d) == 0)
            {
                std::cerr << "RestrictionDecoder::build_connected_components: got detector " << d 
                            << " without corresponding dfs graph vertex" << _die{};
            }
        }
        auto v = vertex_id_map[m.x],
             w = vertex_id_map[m.y];
        auto comm_edges = edges_incident_to({v,w});
        if (comm_edges.empty())
            gr.add_edge({v,w}, restriction::dfs::e_data_type{});
    }

    const auto br = vertex_id_map[boundary_id(C_RED)];
    if (gr.degree(br) == 0)
        return;

    // search for connected components via dfs
    std::vector<det_id_type> cc_dsts;
    cc_dsts.reserve(4);

    std::vector<hg::id_type> st{br};
    std::vector<std::optional<hg::id_type>> prev(vertex_id_map.size(), std::nullopt);
    prev[br] = br;
    while (st.size() > 0)
    {
        auto v = st.back();
        st.pop_back();
        for (const auto& [w, __unused_mul] : gr.adjacency(v))
        {
            if (prev[w].has_value())
                continue;
            prev[w] = v;
            if (is_boundary_id(gr.v(w).d))
            {
                cc_dsts.push_back(w);
                continue;
            }
        }
    }

    // resolve components using `prev`
    for (auto v : cc_dsts)
    {
        std::vector<det_id_type> path;
        path.reserve(8);
        while (prev[v] != v)
        {
            path.push_back(gr.v(v).d);
            v = prev[v];
        }
        path.push_back(v);

        Color c = complementary_color(hg_->v(path.front()).color, hg_->v(path.back()).color);
        components.push_back( restriction::ConnectedComponent{path, c} );
    }

    return components;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec
