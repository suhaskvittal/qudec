/*
 *  author: Suhas Vittal
 *  date:   9 October 2025
 */

#ifndef HYPERGRAPH_h
#define HYPERGRAPH_h

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

using GRAPH_COMPONENT_ID = int32_t;

template <class VERTEX_DATA_TYPE, class EDGE_DATA_TYPE, size_t MAX_ORDER=3>
class HYPERGRAPH
{
public:
    static_assert(MAX_ORDER > 1, "MAX_ORDER must be greater than 1");

    using id_type = GRAPH_COMPONENT_ID;

    struct VERTEX
    {
        id_type          id;
        VERTEX_DATA_TYPE data;
    };

    struct EDGE
    {
        using vertex_list_type = std::array<VERTEX*, MAX_ORDER>;

        vertex_list_type  vertices{};
        size_t            order;
        EDGE_DATA_TYPE    data;
    };

    using adjacency_list_entry = std::pair<
                                        VERTEX*,
                                        std::conditional_t<MAX_ORDER == 2, VERTEX*, std::vector<VERTEX*>>
                                        >;    
    using adjacency_list = std::vector<adjacency_list_entry>;
protected:
    std::vector<VERTEX*> vertices_;
    std::vector<EDGE*>   edges_;

    std::unordered_map<id_type, VERTEX*> vertex_id_map_;
    std::unordered_map<VERTEX*, adjacency_list> adjacency_;
public:
    HYPERGRAPH(size_t reserve_vertices=1024, size_t reserve_edges=4096);
    ~HYPERGRAPH();

    VERTEX*                     add_vertex(id_type, VERTEX_DATA_TYPE);
    template <class ITER> EDGE* add_edge(ITER v_begin, ITER v_end, EDGE_DATA_TYPE);

    VERTEX* get_vertex(id_type) const;
    void remove_vertex(VERTEX*);

    template <class ITER> EDGE*               get_edge_and_fail_if_nonunique(ITER v_begin, ITER v_end);
    template <class ITER> std::vector<EDGE*>  get_all_incident_edges(ITER v_begin, ITER v_end);
    void remove_edge(EDGE*);

    const std::vector<VERTEX*>& get_vertices() const { return vertices_; }
    const std::vector<EDGE*>&   get_edges() const { return edges_; }
    const adjacency_list&       get_adjacency_list(VERTEX* v) const { return adjacency_.at(v); }
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#include "hypergraph.tpp"

#endif  // HYPERGRAPH_h
