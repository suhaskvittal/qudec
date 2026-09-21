/*
 *  author: Suhas Vittal
 *  date:   21 September 2026
 * */

#ifndef DECODER_BAL_HYPERGRAPH_h
#define DECODER_BAL_HYPERGRAPH_h

#include <array>
#include <cstdint>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace decoder
{
namespace bal
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * `VDATA` = metadata for each vertex
 * `EDATA` = metadata for each hyperedge
 * `K` = hypergraph max order
 * */
template <class VDATA, class EDATA, size_t K>
class Hypergraph
{
public:
    using id_type = int32_t;
    using edge_support_type = std::array<id_type, K>;
    
    constexpr static id_type INV{-1};

    /*
     * When `Hypergraph` has `K = 2`, then we can track adjacency
     * using a hashmap that maps vertices to the corresponding edge.
     *  -- for example:
     *          v1: [ {v2, A}, {v4, B} ]
     *
     * When `K > 2` however, then this does not work as two vertices
     * may share multiple hyperedges. Then, we must use an interaction
     * graph.
     * */
private:
    using si_adj_elem_type = std::unordered_map<id_type, id_type>;
    using si_adj_type = std::vector<si_adj_elem_type>;

    struct hy_adj_vdata_type { };
    struct hy_adj_edata_type { std::vector<id_type> incident{}; };
    using hy_adj_type = Hypergraph<hy_adj_vdata_type, hy_adj_edata_type, 2>;
public:
    using adjacency_type = std::conditional_t< (K==2), si_adj_type, hy_adj_type >;
private:
    std::vector<VDATA> vertex_data_;
    std::vector<EDATA> edge_data_;

    /*
     * Adjacency structures:
     * */
    adjacency_type adj_;
    std::vector<std::vector<id_type>> incident_edges_;  // indexed by vertex id
    std::vector<edge_support_type> edge_support_;       // indexed by edge id

    size_t vertex_count_{0},
           edge_count_{0};
public:
    Hypergraph() =default;
    Hypergraph(const Hypergraph&) =default;

    /*
     * Graph updates:
     * */
    id_type add_vertex(VDATA);
    id_type add_edge(std::vector<id_type>, EDATA);

    /*
     * `for_each_edge_incident_to()` calls the given callback for
     * each edge that is incident to *all* specified vertices. The
     * callback is given the edge id. Data like edge support or metadata
     * can be obtained using the accessors below.
     *
     * When the number of input vertices is 1 or 2, this function is fastest and is
     * O(1). Otherwise, the runtime is linear in the degree.
     * */
    template <class CALLBACK>
    void for_each_edge_incident_to(this auto&, std::vector<id_type>, const CALLBACK&);

    /*
     * `edges_incident_to()` calls `for_each_edge_incident_to()`.
     * */
    std::vector<id_type> edges_indicent_to(std::vector<id_type>) const;

    /*
     * Returns all neighbors of the given vertex and their multiplicity (second argument in pair).
     * Note that multiplicity is always 1 for 2-graphs.
     * */
    std::vector<std::pair<id_type, size_t>> adjacency(id_type) const;

    /*
     * One-liners:
     * */
    size_t vertex_count() const { return vertex_count_; }
    size_t edge_count() const { return edge_count_; }

    size_t N() const { return vertex_count(); }
    size_t M() const { return edge_count(); }

    auto& v(this auto& g, id_type x) { return g.vertex_data_[x]; }
    auto& e(this auto& g, id_type x) { return g.edge_data_[x]; }

    size_t degree(id_type v) const { return incident_edges_[v].size(); }
    size_t order(id_type e) const
        { 
            auto begin = support(e).begin(), end = support(e).end();
            auto it = std::find(begin, end, INV);
            return std::distance(begin, it);
        }

    const edge_support_type& support(id_type x) const { return edge_support_[x]; }
private:
    void validate_vertex_list(std::string_view caller_id, const std::vector<id_type>&);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace bal
} // namespace decoder

#include "hypergraph.tpp"

#endif
