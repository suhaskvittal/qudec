/*
 *  author: Suhas Vittal
 *  date:   9 October 2025
 */

#define TEMPL_PARAMS    template <class V_t, class E_t, bool IM>
#define TEMPL_CLASS     HYPERGRAPH<V_t, E_t, IM>

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

TEMPL_PARAMS typename TEMPL_CLASS::vertex_ptr
TEMPL_CLASS::add_vertex(id_type id, VERTEX_DATA_TYPE data)
{
    if constexpr (IM)
        throw std::runtime_error("cannot modify immutable hypergraph");
    
    vertex_ptr v{new vertex_type{id, data}};
    vertices_.push_back(v);
    vertex_id_map_[id] = v;
    return v;
}

TEMPL_PARAMS typename TEMPL_CLASS::edge_ptr
TEMPL_CLASS::add_edge(id_type id, std::vector<vertex_ptr> vertex_list, EDGE_DATA_TYPE data)
{
    if constexpr (IM)
        throw std::runtime_error("cannot modify immutable hypergraph");

    edge_ptr e{new edge_type{id, vertex_list, data}};
    edges_.push_back(e);

    for (auto* v : vertex_list)
        incidence_map_[v].push_back(e);

    edge_id_map_[id] = e;
    return e;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS typename TEMPL_CLASS::vertex_ptr
TEMPL_CLASS::get_vertex(id_type id) const
{
    return vertex_id_map_.at(id);
}

TEMPL_PARAMS typename TEMPL_CLASS::edge_ptr
TEMPL_CLASS::get_edge(id_type id) const
{
    return edge_id_map_.at(id);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS void
TEMPL_CLASS::remove_vertex(vertex_ptr v)
{
    if constexpr (IM)
        throw std::runtime_error("cannot modify immutable hypergraph");

    auto it = std::find(vertices_.begin(), vertices_.end(), v);
    if (it == vertices_.end())
        throw std::runtime_error("vertex not found");
    vertices_.erase(it);

    // remove all edges that contain this vertex -- do a copy as we will modify the original
    // entry
    incidence_array_type incidence_list(incidence_map_[v]);
    for (auto* e : incidence_list)
        remove_edge(e);

    // remove the vertex from the incidence map
    vertex_id_map_.erase(v->id);
    incidence_map_.erase(v);
    delete v;
}

TEMPL_PARAMS void
TEMPL_CLASS::remove_edge(edge_ptr e)
{
    if constexpr (IM)
        throw std::runtime_error("cannot modify immutable hypergraph");

    auto it = std::find(edges_.begin(), edges_.end(), e);
    if (it == edges_.end())
        throw std::runtime_error("edge not found");
    edges_.erase(it);

    for (auto* v : e->vertices)
    {
        auto& inc = incidence_map_[v];
        auto it = std::find(inc.begin(), inc.end(), e);
        inc.erase(it);
    }

    edge_id_map_.erase(e->id);
    delete e;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAMS std::vector<typename TEMPL_CLASS::edge_ptr>
TEMPL_CLASS::get_common_edges(std::vector<vertex_ptr> vertex_list)
{
    std::vector<edge_ptr> common_edges;

    vertex_ptr v0 = vertex_list[0];
    for (edge_ptr e : incidence_map_[v0])
    {
        bool has_all = std::all_of(vertex_list.begin()+1, vertex_list.end(),
                                [e] (vertex_ptr v) 
                                { 
                                    return std::find(e->vertices.begin(), e->vertices.end(), v) != e->vertices.end();
                                });
        if (has_all)
            common_edges.push_back(e);
    }
    return common_edges;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
