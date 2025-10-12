/*
 *  author: Suhas Vittal
 *  date:   9 October 2025
 */

#ifndef HYPERGRAPH_h
#define HYPERGRAPH_h

#include <cstdint>
#include <unordered_map>
#include <vector>

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class VERTEX_DATA_TYPE, class EDGE_DATA_TYPE, bool IMMUTABLE=false>
class HYPERGRAPH
{
public:
    using id_type = uint32_t;

    struct vertex_type
    {
        id_type          id;
        VERTEX_DATA_TYPE data;
    };
    
    using vertex_ptr = vertex_type*;

    struct edge_type
    {
        id_type                 id;
        std::vector<vertex_ptr> vertices{};
        EDGE_DATA_TYPE          data;
    };

    using edge_ptr = edge_type*;

    constexpr static id_type INVALID{-1};
protected:
    using incidence_array_type = std::vector<edge_ptr>;

    std::vector<vertex_ptr> vertices_;
    std::vector<edge_ptr>   edges_;

    std::unordered_map<id_type, vertex_ptr> vertex_id_map_;
    std::unordered_map<id_type, edge_ptr>   edge_id_map_;

    std::unordered_map<vertex_ptr, incidence_array_type> incidence_map_;
public:
    HYPERGRAPH(size_t reserve_vertices=1024, size_t reserve_edges=4096);
    ~HYPERGRAPH();

    HYPERGRAPH<VERTEX_DATA_TYPE, EDGE_DATA_TYPE, true> immutable_copy() const;

    vertex_ptr add_vertex(id_type, VERTEX_DATA_TYPE);
    edge_ptr   add_edge(id_type, std::vector<vertex_ptr>, EDGE_DATA_TYPE);

    vertex_ptr get_vertex(id_type) const;
    edge_ptr   get_edge(id_type) const;

    void remove_vertex(vertex_ptr);
    void remove_edge(edge_ptr);

    std::vector<edge_ptr> get_common_edges(std::vector<vertex_ptr>);

    const std::vector<vertex_ptr>& get_vertices() const { return vertices_; }
    const std::vector<edge_ptr>&   get_edges() const { return edges_; }
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#include "hypergraph.tpp"

#endif  // HYPERGRAPH_h
