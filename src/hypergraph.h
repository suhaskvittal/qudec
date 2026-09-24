/*
 *  author: OpenAI GPT-6-Luna
 *  purpose: Define the reusable Hypergraph template in the global namespace.
 * */

#ifndef HYPERGRAPH_H
#define HYPERGRAPH_H

#include "globals.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

/*
 * `VDATA` = metadata for each vertex
 * `EDATA` = metadata for each hyperedge
 * `K` = hypergraph max order
 */
template <class VDATA, class EDATA, size_t K>
class Hypergraph
{
public:
    using id_type = int32_t;
    using edge_support_type = std::array<id_type, K>;
    constexpr static id_type INV{-1};

private:
    using si_adj_elem_type = std::unordered_map<id_type, id_type>;
    using si_adj_type = std::vector<si_adj_elem_type>;

    struct hy_adj_vdata_type { };
    struct hy_adj_edata_type { std::vector<id_type> incident{}; };
    using hy_adj_type = Hypergraph<hy_adj_vdata_type, hy_adj_edata_type, 2>;

public:
    using adjacency_type = std::conditional_t<(K == 2), si_adj_type, hy_adj_type>;

private:
    std::vector<VDATA> vertex_data_;
    std::vector<EDATA> edge_data_;
    adjacency_type adj_;
    std::vector<std::vector<id_type>> incident_edges_;
    std::vector<edge_support_type> edge_support_;
    size_t vertex_count_{0}, edge_count_{0};

public:
    Hypergraph() = default;
    Hypergraph(const Hypergraph&) = default;
    Hypergraph(Hypergraph&&) = default;
    Hypergraph& operator=(const Hypergraph&) = default;
    Hypergraph& operator=(Hypergraph&&) = default;

    id_type add_vertex(VDATA);
    id_type add_edge(std::vector<id_type>, EDATA);

    template <class CALLBACK>
    void for_each_edge_incident_to(this auto&, std::vector<id_type>, const CALLBACK&);

    std::vector<id_type> edges_incident_to(std::vector<id_type>) const;
    std::vector<std::pair<id_type, size_t>> adjacency(id_type) const;

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
    constexpr static size_t max_order() { return K; }

private:
    void validate_vertex_list(std::string_view caller_id, const std::vector<id_type>&) const;
};

#include "hypergraph.tpp"

#endif // HYPERGRAPH_H
