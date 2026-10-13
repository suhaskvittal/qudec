#ifndef QUDEC_LIBRA_ASSIGNMENT_PYMATCHING_H
#define QUDEC_LIBRA_ASSIGNMENT_PYMATCHING_H

#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace dec::libra_assignment {

struct Pair {
    uint32_t a;
    uint32_t b;
    double weight;
};

struct Result {
    // -1 means the independent single channel; otherwise the paired vertex.
    std::vector<int32_t> mate;
    std::vector<size_t> pair_indices;
    int64_t quantized_cost = 0;
};

inline int64_t quantize(double weight) {
    constexpr double scale = 1e6;
    if (!std::isfinite(weight) || weight < 0 ||
        weight > double(std::numeric_limits<int32_t>::max()) / scale)
        throw std::invalid_argument("Invalid Libra assignment weight");
    return std::llround(weight * scale);
}

// Solve the assignment used by Libra: every matched CSS edge must be covered
// exactly once, either by its single channel or by a correlated pair channel.
//
// Vertex i is fired. A boundary edge costs C+s_i, and a pair edge costs
// 2C+p_ij, where C=1+max(s_i), using the original 1e-6 integer quantization.
// An odd-parity solution with degree >=3 at v has at least two pair edges
// (v,a),(v,b). Remove these and toggle the boundary edges at a and b.
// This preserves all parities and lowers cost by at least 2C-s_a-s_b >0.
// Thus every optimum has degree exactly one at every vertex. The n*C term
// is then constant, so the original assignment objective is preserved.
//
// Blocks of at most eight vertices use exact subset dynamic programming to
// avoid matching-engine setup on these tiny assignments. Larger blocks use
// PyMatching. All optimal paths are direct edges. Its compressed endpoint API
// therefore supplies the assignment without a SearchFlooder or Dijkstra path
// reconstruction. The private dummy vertices and zero-weight dummy clique of
// the earlier Blossom V formulation are unnecessary.
inline Result solve(const std::vector<double>& singles, const std::vector<Pair>& pairs) {
    const size_t n = singles.size();
    if (n > size_t(std::numeric_limits<int32_t>::max()))
        throw std::invalid_argument("Libra assignment vertex index overflow");
    Result result;
    result.mate.assign(n, -1);
    if (!n) {
        if (!pairs.empty())
            throw std::invalid_argument("Pairs in an empty Libra assignment");
        return result;
    }

    std::vector<int64_t> single_costs;
    single_costs.reserve(n);
    int64_t largest_single = 0;
    for (double w : singles) {
        const int64_t cost = quantize(w);
        single_costs.push_back(cost);
        largest_single = std::max(largest_single, cost);
    }
    const int64_t penalty = largest_single + 1;
    auto even_edge_weight = [](int64_t cost) {
        // Even integer weights make all Sparse Blossom collision times integral.
        if (cost < 0 || cost > std::numeric_limits<int32_t>::max() / 2)
            throw std::overflow_error(
                "Libra assignment exceeds PyMatching's 32-bit edge range at 1e-6 precision");
        return static_cast<pm::signed_weight_int>(2 * cost);
    };

    struct BestPair {int64_t cost; size_t index;};
    std::map<std::pair<uint32_t, uint32_t>, BestPair> best_pairs;
    for (size_t index = 0; index < pairs.size(); ++index) {
        const auto& pair = pairs[index];
        if (pair.a >= n || pair.b >= n || pair.a == pair.b)
            throw std::invalid_argument("Invalid Libra assignment pair");
        const int64_t cost = quantize(pair.weight);
        if (cost >= single_costs[pair.a] + single_costs[pair.b])
            continue; // Two independent singles are at least as cheap.
        const auto endpoints = std::minmax(pair.a, pair.b);
        auto [it, inserted] = best_pairs.emplace(
            std::pair<uint32_t, uint32_t>{endpoints.first, endpoints.second},
            BestPair{cost, index});
        if (!inserted && cost < it->second.cost)
            it->second = {cost, index};
    }

    if (n <= 8) {
        // best[mask] is the exact minimum cover of the indicated vertices.
        // Cover the first vertex by its single or by any available pair.
        // Removing those vertices gives a smaller mask, already evaluated.
        constexpr int64_t absent = std::numeric_limits<int64_t>::max() / 4;
        std::array<std::array<int64_t, 8>, 8> pair_costs;
        std::array<std::array<size_t, 8>, 8> pair_indices;
        for (auto& row : pair_costs) row.fill(absent);
        for (const auto& [endpoints, pair] : best_pairs) {
            pair_costs[endpoints.first][endpoints.second] = pair.cost;
            pair_indices[endpoints.first][endpoints.second] = pair.index;
        }
        std::array<int64_t, 256> best;
        std::array<int8_t, 256> choice;
        best[0] = 0;
        const uint32_t full_mask = (uint32_t{1} << n) - 1;
        for (uint32_t mask = 1; mask <= full_mask; ++mask) {
            const uint32_t a = std::countr_zero(mask);
            const uint32_t remaining = mask ^ (uint32_t{1} << a);
            best[mask] = single_costs[a] + best[remaining];
            choice[mask] = -1;
            for (uint32_t b = a + 1; b < n; ++b) {
                if (!(remaining & (uint32_t{1} << b))) continue;
                const int64_t cost = pair_costs[a][b] + best[remaining ^ (uint32_t{1} << b)];
                if (cost < best[mask]) {
                    best[mask] = cost;
                    choice[mask] = b;
                }
            }
        }
        result.quantized_cost = best[full_mask];
        for (uint32_t mask = full_mask; mask;) {
            const uint32_t a = std::countr_zero(mask);
            const int32_t b = choice[mask];
            mask ^= uint32_t{1} << a;
            if (b < 0) continue;
            mask ^= uint32_t{1} << b;
            result.mate[a] = b;
            result.mate[b] = a;
            result.pair_indices.push_back(pair_indices[a][b]);
        }
        return result;
    }

    pm::MatchingGraph graph(n, 0, 1.0);
    for (size_t i = 0; i < n; ++i)
        graph.add_boundary_edge(i, even_edge_weight(penalty + single_costs[i]), {});
    for (const auto& [endpoints, pair] : best_pairs)
        graph.add_edge(endpoints.first, endpoints.second,
                       even_edge_weight(2 * penalty + pair.cost), {});

    pm::Mwpm matcher{pm::GraphFlooder(std::move(graph))};
    std::vector<uint64_t> events(n);
    std::iota(events.begin(), events.end(), uint64_t{0});
    pm::decode_detection_events_to_match_edges(matcher, events);

    const auto* nodes = matcher.flooder.graph.nodes.data();
    std::vector<bool> covered(n, false);
    for (const auto& edge : matcher.flooder.match_edges) {
        if (!edge.loc_from)
            throw std::runtime_error("Libra assignment returned an empty match origin");
        const size_t a = edge.loc_from - nodes;
        if (a >= n || covered[a])
            throw std::runtime_error("Libra assignment covers a vertex twice");
        covered[a] = true;
        if (!edge.loc_to) {
            result.quantized_cost += single_costs[a];
            continue;
        }
        const size_t b = edge.loc_to - nodes;
        if (b >= n || covered[b])
            throw std::runtime_error("Libra assignment covers a vertex twice");
        const std::pair<uint32_t, uint32_t> endpoints{
            static_cast<uint32_t>(std::min(a, b)),
            static_cast<uint32_t>(std::max(a, b))};
        auto found = best_pairs.find({endpoints.first, endpoints.second});
        if (found == best_pairs.end())
            throw std::runtime_error("Libra assignment returned a nonadjacent pair");
        covered[b] = true;
        result.mate[a] = b;
        result.mate[b] = a;
        result.pair_indices.push_back(found->second.index);
        result.quantized_cost += found->second.cost;
    }
    if (std::find(covered.begin(), covered.end(), false) != covered.end())
        throw std::runtime_error("Libra assignment left a vertex uncovered");
    return result;
}

} // namespace dec::libra_assignment

#endif
