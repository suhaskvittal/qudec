/*
 *  author: Suhas Vittal
 *  date:   21 September 2026
 * */

#include <algorithm>
#include <deque>
#include <iostream>
#include <unordered_set>

#define TEMPL_PARAM template <class V, class E, size_t K>
#define TEMPL_CLASS Hypergraph<V,E,K>

namespace decoder
{
namespace bal
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM typename TEMPL_CLASS::id_type
TEMPL_CLASS::add_vertex(V d)
{
    const id_type x = static_cast<id_type>(vertex_count());
    vertex_data_.push_back(d);
    incident_edges_.push_back({});
    if constexpr (K == 2)
        adj_.push_back(si_adj_elem_type{});
    else
        adj_.add_vertex(hy_adj_vdata_type{});
    vertex_count_++;
    return x;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM typename TEMPL_CLASS::id_type
TEMPL_CLASS::add_edge(std::vector<id_type> inc, E d)
{
    // check that all elements of `inc` are valid vertices:
    validate_vertex_list("add_edge()", inc);
    // canonicalize: an edge's identity does not depend on the order or
    // repetition of its endpoints.
    std::sort(inc.begin(), inc.end());
    inc.erase(std::unique(inc.begin(), inc.end()), inc.end());
    // check that order is ok as well:
    if (inc.size() > K)
        std::cerr << "Hypergraph::add_edge(): tried to add edge with order > max order (" << K << ")" << _die{};
    // finally, check that the edge is unique:
    if (edges_incident_to(inc).size() > 0)
        std::cerr << "Hypergraph::add_edge(): edge is non-unique" << _die{};

    // create edge support:
    edge_support_type supp{};
    supp.fill(INV);
    std::copy(inc.begin(), inc.end(), supp.begin());

    const id_type x = static_cast<id_type>(edge_count());
    edge_data_.push_back(d);
    edge_support_.push_back(supp);
    edge_count_++;

    if constexpr (K == 2)
    {
        if (inc.size() == 2)
        {
            id_type v = inc[0],
                    w = inc[1];
            adj_[v].insert({w, x});
            adj_[w].insert({v, x});
        }
    }
    else
    {
        for (size_t i = 0; i < inc.size(); i++)
        {
            id_type v = inc[i];
            for (size_t j = i+1; j < inc.size(); j++)
            {
                id_type w = inc[j];
                // two possibiltiies: (v,w) already in `adj_` or it is not:
                auto edges = adj_.edges_incident_to({v,w});
                id_type y;
                if (edges.empty())  // create new edge:
                    y = adj_.add_edge({v,w}, hy_adj_edata_type{});
                else
                    y = edges[0];
                adj_.e(y).incident.push_back(x);
            }
        }
    }

    for (auto v : inc)
        incident_edges_[v].push_back(x);

    return x;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM template <class CALLBACK> void
TEMPL_CLASS::for_each_edge_incident_to(this auto& G, std::vector<id_type> inc, const CALLBACK& cb)
{
    // validate that all of `inc` is valid:
    G.validate_vertex_list("for_each_edge_incident_to()", inc);

    if (inc.empty())
        std::cerr << "Hypergraph::for_each_edge_incident_to(): vertex list is empty" << _die{};

    // simple case:
    if (inc.size() > K)
        std::cerr << "Hypergraph::for_each_edge_incident_to(): tried to query order > max order (" << K << ")" << _die{};

    if (inc.size() == 1)
    {
        for (auto e : G.incident_edges_[inc[0]])
            cb(e);
        return;
    }

    if (inc.size() == 2)
    {
        id_type v = inc[0],
                w = inc[1];
        if constexpr (K == 2)
        {
            auto it = G.adj_[v].find(w);
            if (it != G.adj_[v].end())
                cb(it->second);
        }
        else
        {
            auto edges = G.adj_.edges_incident_to(inc);
            if (edges.empty())
                return;
            for (auto e : G.adj_.e(edges[0]).incident)
                cb(e);
        }
        return;
    }

    // K > 2: basic idea is to get an initial list using the interaction graph.
    // Then, we can refine this list iteratively.

    if constexpr (K != 2)
    {
        std::vector<id_type> out = G.adj_.edges_incident_to({inc[0], inc[1]});
        auto it = std::remove_if(out.begin(), out.end(),
                            [&G, &inc] (auto e)
                            {
                                const auto& _supp = G.support(e);
                                std::unordered_set<id_type> supp(_supp.begin(), _supp.end());
                                const bool all_present = std::all_of(inc.begin()+2, inc.end(),
                                                                [&supp] (id_type x) { return supp.count(x) > 0; });
                                return !all_present;
                            });
        out.erase(it, out.end());
        for (auto e : out)
            cb(e);
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM std::vector<typename TEMPL_CLASS::id_type>
TEMPL_CLASS::edges_incident_to(std::vector<id_type> inc) const
{
    std::vector<id_type> out;
    for_each_edge_incident_to(inc, [&out] (auto e) { out.push_back(e); });
    return out;
}

TEMPL_PARAM std::vector<std::pair<typename TEMPL_CLASS::id_type, size_t>>
TEMPL_CLASS::adjacency(id_type v) const
{
    using entry_type = std::pair<id_type, size_t>;

    std::unordered_map<id_type, size_t> mmap;
    mmap.reserve(degree(v) * K);

    // initialize mmap
    for (auto e : incident_edges_[v])
        for (auto w : support(e))
            if (w != v && w != INV)
                mmap[w]++;

    // copy mmap contents to `out`
    std::vector<entry_type> out;
    out.reserve(mmap.size());
    for (const auto& kv : mmap)
        out.emplace_back(kv.first, kv.second);
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM void
TEMPL_CLASS::validate_vertex_list(std::string_view caller_id, const std::vector<id_type>& vlist) const
{
    for (id_type x : vlist)
        if (x < 0 || static_cast<size_t>(x) >= vertex_count())
            std::cerr << "Hypergraph::" << caller_id << ": vertex \"" << x << "\" not in hypergraph (N = " << N() << ")" << _die{};
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace hg
{

template <class G> ErrorCube
extract_error_cube(const G& gr, id_type base, size_t k)
{
    ErrorCube cu{ .base=base, .radius=k };
    cu.nodes.reserve(32);
    cu.edges.reserve(512);
    cu.dist_map.reserve(32);

    cu.dist_map[base] = ErrorCube::dist_type{};

    std::unordered_set<id_type> visited_edges;
    visited_edges.reserve(512);

    std::deque<id_type> bfsq{ base }; 
    while (bfsq.size() > 0)
    {
        id_type v = std::move(bfsq.front());
        bfsq.pop_front();

        auto d_v = cu.dist_map[v];
        cu.nodes.push_back(v);

        if (d_v.edges.size() >= k)
            continue;

        gr.for_each_edge_incident_to({v},
                [&] (auto e)
                {
                    if (visited_edges.count(e))
                        return;
                    auto supp = gr.support(e);

                    // add `e` to `cu.edges` since this is a potential edge
                    // in the graph
                    cu.edges.push_back(e);
                    visited_edges.insert(e);

                    // now handle traversal
                    auto d_w = d_v;
                    d_w.edges.push_back(e);
                    d_w.pr *= gr.e(e).error_probability;

                    for (auto w : supp)
                    {
                        auto d_w_it = cu.dist_map.find(w);
                        if (d_w_it != cu.dist_map.end())
                        {
                            if (d_w.pr > d_w_it->second.pr)
                                d_w_it->second = d_w;        
                            continue;
                        }
                        else
                        {
                            cu.dist_map.insert({w, d_w});
                        }
                        bfsq.push_back(w);
                    }
                });
    }
    return cu;
}

} // namespace hg

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace bal
} // namespace decoder

#undef TEMPL_PARAM
#undef TEMPL_CLASS
