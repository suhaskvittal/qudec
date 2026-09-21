/*
 *  author: Suhas Vittal
 *  date:   21 September 2026
 * */

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
    [[ maybe_unused ]] const id_type x = vertex_count();
    vertex_data_.push_back(d);
    incident_edges.push_back({});
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
    // check that order is ok as well:
    if (inc.size() > K)
        std::cerr << "Hypergraph::add_edge(): tried to add edge with order > max order (" << K << ")" << _die{};
    // finally, check that the edge is unique:
    if (edges_incident_to(inc).size() > 0)
        std::cerr << "Hypergraph::add_edge(): edge is non-unique" << _die{};

    // create edge support:
    edge_support_type supp{};
    supp.fill(id_type{-1});
    std::copy(inc.begin(), inc.end(), supp.begin());

    [[ maybe_unused ]] const id_type x = edge_count();
    edge_data_.push_back(d);
    edge_support_.push_back(supp);
    edge_count_++;
    
    if constexpr (K == 2)
    {
        id_type v = inc[0],
                w = inc[1];
        adj_[v].insert({w, x});
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
                auto edges = edges_incident_to({v,w});
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
        incident_edges_.push_back(x);

    return x;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM template <class CALLBACK> void
TEMPL_CLASS::for_each_edge_incident_to(this auto& G, std::vector<id_type> inc, const CALLBACK& cb)
{
    // validate that all of `inc` is valid:
    G.validate_vertex_list("for_each_edge_incident_to()", inc);
    
    // simple case:
    if (inc.size() > K)
        return;

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
            cb(adj_[v][w]);
        }
        else
        {
            auto x = G.adj_.edges_incident_to(inc)[0];
            for (auto e : G.adj_.e(x).incident)
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
                            [&G, &inc, &check_for_membership] (auto e)
                            {
                                const auto& _supp = G.support(e);
                                std::unordered_set<id_type> supp(_supp.begin(), _supp.end());
                                const bool all_present = std::all_of(inc.begin()+2, inc.end(),
                                                                [&supp] (id_type x) { return supp.count(x); });
                            });
        for (auto e : out)
            cb(e);
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM std::vector<typename TEMPL_CLASS::id_type>
TEMPL_CLASS::edges_incident_to(std::vector<id_type> inc)
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
    mmap.reserve(degree(v));

    // initialize mmap
    for (auto e : incident_edges_[v])
        for (auto w : support(e))
            if (w != v)
                mmap[w]++;

    // copy mmap contents to `out` -- both value types should be `entry_type`
    std::vector<entry_type> out(mmap.size());
    std::move(mmap.begin(), mmap.end(), out.begin());
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM void
TEMPL_CLASS::validate_vertex_list(std::string_view caller_id, const std::vector<id_type>& vlist)
{
    for (id_type x : vlist)
        if (x >= vertex_count())
            std::cerr << "Hypergraph::" << caller_id << ": vertex \"" << x << "\" not in hypergraph (N = " << N() << ")" << _die{};
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace bal
} // namespace decoder

#undef TEMPL_PARAM
#undef TEMPL_CLASS
