/* Contributor: OpenAI Codex. Seed conditioning and canonical-path filtering. */
#include "decoder/seeded.h"
#include <pymatching/sparse_blossom/driver/user_graph.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <queue>
#include <stdexcept>
#include <tuple>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec
{
namespace
{
using Clock = std::chrono::steady_clock;
using Weight = int64_t;
constexpr Weight infinity = std::numeric_limits<Weight>::max() / 8;

uint64_t elapsed_ns(Clock::time_point start)
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
}
}

CorrelatedPyMatching::CorrelatedPyMatching(const stim::DetectorErrorModel& dem)
    : num_detectors(dem.count_detectors()), num_observables(dem.count_observables()),
      matcher_(pm::detector_error_model_to_mwpm(dem, pm::NUM_DISTINCT_WEIGHTS, true, true))
{
    if (num_observables > 64)
        throw std::invalid_argument("Seeded correlated PyMatching supports at most 64 observables");
}

CorrelatedPyMatching::~CorrelatedPyMatching()
{
    restore();
}

void CorrelatedPyMatching::prepare(const std::vector<int64_t>& selected_edges)
{
    if (prepared_ || !matcher_.flooder.graph.previous_weights.empty() ||
        !matcher_.search_flooder.graph.previous_weights.empty())
        throw std::logic_error("Correlated PyMatching graph is already conditioned");
    prepared_ = true;
    try
    {
        matcher_.flooder.graph.reweight_for_edges(selected_edges);
        matcher_.search_flooder.graph.reweight_for_edges(selected_edges);
    }
    catch (...)
    {
        restore();
        throw;
    }
}

Result CorrelatedPyMatching::decode_prepared(SyndromeRef syndrome, ObsRef)
{
    if (!prepared_)
        throw std::logic_error("Seed correction must condition PyMatching before its final pass");
    const auto start = Clock::now();
    std::vector<uint64_t> events;
    events.reserve(syndrome.popcnt());
    for (size_t i = 0; i < num_detectors; ++i)
        if (syndrome[i]) events.push_back(i);
    const auto decoded = pm::decode_detection_events_for_up_to_64_observables(matcher_, events, false);
    Result result;
    result.flipped_obs = ObsType(num_observables);
    for (size_t i = 0; i < num_observables; ++i)
        result.flipped_obs[i] = (decoded.obs_mask >> i) & 1;
    result.matching_data.frame_flips = result.flipped_obs;
    result.decode_latency_ns = elapsed_ns(start);
    ++calls_;
    total_ns_ += result.decode_latency_ns;
    return result;
}

void CorrelatedPyMatching::restore() noexcept
{
    if (!prepared_) return;
    matcher_.flooder.graph.undo_reweights();
    matcher_.search_flooder.graph.undo_reweights();
    prepared_ = false;
}

void CorrelatedPyMatching::print_stats(std::ostream& out) const
{
    out << "Seeded correlated PyMatching final pass: calls=" << calls_
        << " average_us=" << (calls_ ? double(total_ns_) / calls_ / 1000.0 : 0) << '\n';
}

void CorrelatedPyMatching::mpi_accumulate()
{
#if defined(ENABLE_MPI)
    uint64_t counters[]{calls_, total_ns_};
    MPI_Allreduce(MPI_IN_PLACE, counters, 2, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    calls_ = counters[0];
    total_ns_ = counters[1];
#endif
}

struct SeededPostprocessor::Impl
{
    struct Edge
    {
        uint32_t u, v;
        Weight original;
        pm::weight_int* live;
        uint64_t observable_mask;
    };
    struct Row
    {
        std::vector<Weight> distance;
        std::vector<uint32_t> parent;
        std::vector<uint32_t> parent_edge;
    };

    const size_t detectors;
    const size_t observables;
    std::vector<Edge> edges;
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> adjacency;
    std::map<std::tuple<int64_t, int64_t, uint64_t>, uint32_t> edge_ids;
    std::vector<uint32_t> component;
    std::vector<std::optional<Row>> rows;

    Impl(pm::Mwpm& matcher, size_t n, size_t obs)
        : detectors(n), observables(obs), adjacency(n + 2), component(n, UINT32_MAX), rows(n)
    {
        if (obs != 1 || matcher.search_flooder.graph.nodes.size() != n)
            throw std::invalid_argument("Seed postprocessing requires a one-observable surface-code graph");
        auto& nodes = matcher.search_flooder.graph.nodes;
        for (size_t u = 0; u < n; ++u)
            for (size_t j = 0; j < nodes[u].neighbors.size(); ++j)
            {
                const auto* neighbor = nodes[u].neighbors[j];
                const int64_t other = neighbor ? neighbor - nodes.data() : -1;
                if (other >= 0 && u > static_cast<size_t>(other)) continue;
                uint64_t mask = 0;
                for (size_t k : nodes[u].neighbor_observable_indices[j]) mask ^= uint64_t(1) << k;
                const uint32_t target = other >= 0 ? static_cast<uint32_t>(other) :
                    static_cast<uint32_t>(n + (mask & 1));
                const uint32_t id = edges.size();
                edges.push_back({static_cast<uint32_t>(u), target,
                                 nodes[u].neighbor_weights[j], &nodes[u].neighbor_weights[j], mask});
                adjacency[u].push_back({target, id});
                // Boundary vertices terminate paths; their adjacency is used
                // only for reachability, never for expansion in Dijkstra.
                if (target < n) adjacency[target].push_back({static_cast<uint32_t>(u), id});
                const auto key = std::tuple<int64_t, int64_t, uint64_t>{
                    static_cast<int64_t>(u), other, mask};
                if (!edge_ids.emplace(key, id).second)
                    throw std::invalid_argument("Ambiguous parallel physical edges in seed graph");
            }
        uint32_t component_id = 0;
        for (size_t root = 0; root < n; ++root)
        {
            if (component[root] != UINT32_MAX) continue;
            std::vector<uint32_t> stack{static_cast<uint32_t>(root)};
            component[root] = component_id;
            while (!stack.empty())
            {
                const uint32_t u = stack.back();
                stack.pop_back();
                for (auto [v, id] : adjacency[u])
                    if (v < n && component[v] == UINT32_MAX)
                    {
                        component[v] = component_id;
                        stack.push_back(v);
                    }
            }
            ++component_id;
        }
    }

    const Row& row(uint32_t source)
    {
        if (rows[source]) return *rows[source];
        Row result;
        result.distance.assign(detectors + 2, infinity);
        result.parent.assign(detectors + 2, UINT32_MAX);
        result.parent_edge.assign(detectors + 2, UINT32_MAX);
        using Item = std::pair<Weight, uint32_t>;
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
        result.distance[source] = 0;
        queue.push({0, source});
        while (!queue.empty())
        {
            const auto [distance, u] = queue.top();
            queue.pop();
            if (distance != result.distance[u] || u >= detectors) continue;
            for (auto [v, id] : adjacency[u])
            {
                if (v < detectors && component[v] != component[source]) continue;
                const Weight next = distance + edges[id].original;
                if (next < result.distance[v])
                {
                    result.distance[v] = next;
                    result.parent[v] = u;
                    result.parent_edge[v] = id;
                    queue.push({next, v});
                }
            }
        }
        rows[source] = std::move(result);
        return *rows[source];
    }

    Weight path_cost(uint32_t source, uint32_t target)
    {
        const Row& route = row(source);
        if (route.distance[target] == infinity) return infinity;
        Weight cost = 0;
        for (uint32_t v = target; v != source; v = route.parent[v])
        {
            if (route.parent_edge[v] == UINT32_MAX)
                throw std::logic_error("Incomplete canonical seed path");
            cost += *edges[route.parent_edge[v]].live;
        }
        return cost;
    }

    uint32_t selected_edge_id(const MatchingData::assignment_type::correction_edge_type& selected) const
    {
        int64_t u = selected.d1, v = selected.d2;
        if (v >= 0 && u > v) std::swap(u, v);
        const auto it = edge_ids.find({u, v, selected.observable_mask});
        if (it == edge_ids.end())
            throw std::invalid_argument("Seed correction edge is absent from the conditioned matching graph");
        return it->second;
    }

    SeededPostprocessor::Output apply(SyndromeRef syndrome, const MatchingData& matching)
    {
        SeededPostprocessor::Output out{SyndromeType(detectors), ObsType(observables)};
        out.remaining = syndrome;
        std::vector<uint32_t> events;
        for (size_t d = 0; d < detectors; ++d)
            if (syndrome[d]) events.push_back(d);
        for (const auto& pair : matching.assignments)
        {
            const int64_t a = pair.d1 >= 0 ? pair.d1 : pair.d2;
            const int64_t b = pair.d1 >= 0 ? pair.d2 : pair.d1;
            if (a < 0 || static_cast<size_t>(a) >= detectors ||
                (b >= 0 && static_cast<size_t>(b) >= detectors))
                throw std::invalid_argument("Seed matched endpoint is outside the graph");
            Weight own = 0;
            bool changed = pair.correction_edges.empty();
            for (const auto& selected : pair.correction_edges)
            {
                const auto& edge = edges[selected_edge_id(selected)];
                own += *edge.live;
                changed |= *edge.live != edge.original;
            }
            if (changed) continue;
            bool unsafe = false;
            for (uint32_t competitor : events)
            {
                if (competitor == a || competitor == b || component[competitor] != component[a])
                    continue;
                if (path_cost(a, competitor) < own ||
                    (b >= 0 && path_cost(b, competitor) < own))
                {
                    unsafe = true;
                    break;
                }
            }
            if (b >= 0)
                for (int64_t endpoint : {a, b})
                    for (size_t boundary = 0; boundary < 2; ++boundary)
                        unsafe |= path_cost(endpoint, detectors + boundary) < own;
            else
                unsafe |= path_cost(a, detectors + !bool(pair.frame_flips[0])) < own;
            if (unsafe) continue;
            out.remaining[a] = false;
            if (b >= 0) out.remaining[b] = false;
            out.fixed_frame ^= pair.frame_flips;
            ++out.frozen_pairs;
        }
        return out;
    }
};

SeededPostprocessor::SeededPostprocessor(pm::Mwpm& matcher, size_t detectors, size_t observables)
    : impl_(std::make_unique<Impl>(matcher, detectors, observables))
{}

SeededPostprocessor::~SeededPostprocessor() = default;

SeededPostprocessor::Output SeededPostprocessor::apply(SyndromeRef syndrome, const MatchingData& matching)
{
    return impl_->apply(syndrome, matching);
}

} // namespace dec
