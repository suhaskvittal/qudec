/*
 *  author: Suhas Vittal
 *  date:   9 October 2025
 */

#include <algorithm>
#include <array>
#include <stdexcept>

#define TEMPL_PARAMS    template <class VERTEX_DATA_TYPE, class EDGE_DATA_TYPE, size_t MAX_ORDER>
#define TEMPL_CLASS     HYPERGRAPH<VERTEX_DATA_TYPE, EDGE_DATA_TYPE, MAX_ORDER>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS
TEMPL_CLASS::HYPERGRAPH(size_t reserve_vertices, size_t reserve_edges)
{
    vertices_.reserve(reserve_vertices);
    edges_.reserve(reserve_edges);
}

TEMPL_PARAMS
TEMPL_CLASS::~HYPERGRAPH()
{
    for (auto* v : vertices_)
        delete v;
    for (auto* e : edges_)
        delete e;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS typename TEMPL_CLASS::VERTEX*
TEMPL_CLASS::add_vertex(id_type id, VERTEX_DATA_TYPE data)
{
    if (vertex_id_map_.find(id) != vertex_id_map_.end())
        throw std::runtime_error("vertex already exists");

    VERTEX* v = new VERTEX{id, data};
    vertices_.push_back(v);
    vertex_id_map_[id] = v;
    return v;
}

TEMPL_PARAMS template <class ITER> typename TEMPL_CLASS::EDGE*
TEMPL_CLASS::add_edge(ITER v_begin, ITER v_end, EDGE_DATA_TYPE data)
{
    size_t order = std::distance(v_begin, v_end);
    if (order < 2)
        throw std::runtime_error("edge must have at least 2 vertices");
    if (order > MAX_ORDER)
        throw std::runtime_error("vertex list size exceeds MAX_ORDER");

    typename EDGE::vertex_list_type vertex_list;
    std::copy(v_begin, v_end, vertex_list.begin());
    EDGE* e = new EDGE{vertex_list, order, data};
    edges_.push_back(e);

    if constexpr (MAX_ORDER == 2)
    {
        VERTEX* v = vertex_list[0];
        VERTEX* w = vertex_list[1];
        adjacency_[v].emplace_back(w, e);
        adjacency_[w].emplace_back(v, e);
    }
    else
    {
        for (size_t i = 0; i < e->order; i++)
        {
            VERTEX* v = e->vertices[i];
            auto& v_adj = adjacency_[v];
            for (size_t j = i+1; j < e->order; j++)
            {
                VERTEX* w = e->vertices[j];
                auto it = std::find_if(v_adj.begin(), v_adj.end(),
                                        [w] (const auto& entry) { return entry.first == w; });
                if (it == v_adj.end())
                    v_adj.emplace_back(w, std::vector<EDGE*>{e});
                else
                    it->second.push_back(e);
            }
        }
    }

    return e;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS typename TEMPL_CLASS::VERTEX*
TEMPL_CLASS::get_vertex(id_type id) const
{
    auto it = vertex_id_map_.find(id);
    return it == vertex_id_map_.end() ? nullptr : it->second;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS void
TEMPL_CLASS::remove_vertex(VERTEX* v)
{
    auto it = std::find(vertices_.begin(), vertices_.end(), v);
    if (it == vertices_.end())
        throw std::runtime_error("vertex not found");
    vertices_.erase(it);

    // remove all edges that contain this vertex -- do a copy as we will modify the original
    // entry
    adjacency_list v_adj = adjacency_[v];
    if constexpr (MAX_ORDER == 2)
    {
        // max_order 2 is simple -- can remove as we see the edges
        for (auto& [w, e] : v_adj)
            remove_edge(e);
    }
    else
    {
        // a bit harder when max_order > 2 as we need to worry about duplicates
        std::unordered_set<EDGE*> incident_edges;  // there may be duplicates, so use a set
        incident_edges.reserve(v_adj.size());
        for (auto& [w, e_list] : v_adj)
        {
            for (auto* e : e_list)
                incident_edges.insert(e);
        }
    
        for (EDGE* e : incident_edges)
            remove_edge(e);
    }

    // remove the vertex from the incidence map
    vertex_id_map_.erase(v->id);
    adjacency_.erase(v);
    delete v;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS template <class ITER> typename TEMPL_CLASS::EDGE*
TEMPL_CLASS::get_edge_and_fail_if_nonunique(ITER v_begin, ITER v_end)
{
    std::vector<EDGE*> edges = get_all_incident_edges(v_begin, v_end);
    if (edges.size() > 1)
        throw std::runtime_error("non-unique edge");
    return edges.empty() ? nullptr : edges[0];
}

TEMPL_PARAMS template <class ITER> std::vector<typename TEMPL_CLASS::EDGE*>
TEMPL_CLASS::get_all_incident_edges(ITER v_begin, ITER v_end)
{
    const size_t v_count = std::distance(v_begin, v_end);
    
    // if `v_count > MAX_ORDER`, then there are no edges that are incident on all vertices
    if (v_count > MAX_ORDER)
        return {};

    std::vector<EDGE*> common_edges;
    common_edges.reserve(4);

    if (v_begin == v_end)
        throw std::runtime_error("empty vertex list");

    VERTEX* v0 = *v_begin;
    const auto& v0_adj = adjacency_[v0];
    if (v_count > 2)
    {
        if constexpr (MAX_ORDER == 2)
        {
            // we just need this to avoid compile errors in this block
            // -- note that we will never reach this code because of the `v_count > MAX_ORDER`
            // check above
            return {};
        }
        else
        {
            VERTEX* v1 = *(v_begin+1);
            
            // we only need to check `v1`'s entry in `adjacency_[v0]` since at least these
            // edges are incident on both `v0` and `v1`
            auto adj_it = std::find_if(v0_adj.begin(), v0_adj.end(),
                                        [v1] (const auto& entry) { return entry.first == v1; });
            const auto& [__unused_w, e_list] = *adj_it;

            // now, make sure `v_begin+2 ... v_end-1` are also incident on all edges in `e_list`
            std::copy_if(e_list.begin(), e_list.end(), std::back_inserter(common_edges),
                        [v_begin, v_end] (EDGE* e) 
                        {
                            return std::all_of(v_begin+2, v_end,
                                        [e] (VERTEX* v) 
                                        { 
                                            return std::find(e->vertices.begin(), e->vertices.end(), v) 
                                                                    != e->vertices.end();
                                        });
                        });
        }
    }
    if (v_count == 2)
    {
        VERTEX* v1 = *(v_begin+1);
        auto it = std::find_if(v0_adj.begin(), v0_adj.end(),
                                [v1] (const auto& entry) { return entry.first == v1; });
        if (it != v0_adj.end())
        {
            if constexpr (MAX_ORDER == 2)
                common_edges.push_back(it->second);
            else
                common_edges = it->second;
        }
    }
    else if (v_count == 1)
    {
        // copy all edges that are incident on v0
        for (const auto& [__unused_w, e_singleton_or_list] : v0_adj)
        {
            if constexpr (MAX_ORDER == 2)
                common_edges.push_back(e_singleton_or_list);
            else
                common_edges.insert(common_edges.end(), e_singleton_or_list.begin(), e_singleton_or_list.end());
        }
    }

    return common_edges;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS void
TEMPL_CLASS::remove_edge(EDGE* e)
{
    auto it = std::find(edges_.begin(), edges_.end(), e);
    if (it == edges_.end())
        throw std::runtime_error("edge not found");
    edges_.erase(it);

    for (auto* v : e->vertices)
    {
        const auto& v_adj = adjacency_[v];
        if constexpr (MAX_ORDER == 2)
        {
            // max_order == 2 : just find the entry and delete it
            auto it = std::find_if(v_adj.begin(), v_adj.end(),
                                    [e] (const auto& entry) { return entry.second == e; });
            if (it != v_adj.end())
                v_adj.erase(it);
        }
        else
        {
            // more involved for max_order > 2: need to find the entry containing the list,
            // and handle the case where the list is empty after deleting `e`
            auto it = std::remove_if(v_adj.begin(), v_adj.end(),
                                    [e] (auto& p)
                                    {
                                        auto& e_list = p.second;
                                        // we will mark `e_list` for removal if it is empty after removing `e`
                                        auto it = std::find(e_list.begin(), e_list.end(), e);
                                        if (it != e_list.end())
                                            e_list.erase(it);
                                        return e_list.empty();
                                    });
            v_adj.erase(it, v_adj.end());
        }
    }

    delete e;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
