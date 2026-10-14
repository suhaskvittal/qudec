/* Contributor: OpenAI Codex. */
#include "decoder/cluster_gap.h"
#include <stim.h>

namespace dec
{
void
validate_cluster_gap_dem(const stim::DetectorErrorModel& dem)
{
    if (dem.count_observables() != 1)
        throw std::invalid_argument("Cluster gap supports exactly one planar memory observable");
    dem.iter_flatten_error_instructions([&](const stim::DemInstruction& instruction)
    {
        if (instruction.arg_data[0] > 0.5)
            throw std::invalid_argument("Cluster gap currently requires nonnegative log-odds weights (p <= 1/2)");
        instruction.for_separated_targets([&](const auto& group)
        {
            size_t detectors = 0, observables = 0;
            for (const auto& target : group)
            {
                detectors += target.is_relative_detector_id();
                observables += target.is_observable_id();
            }
            if (detectors > 2 || (observables && (detectors != 1 || observables != 1)))
                throw std::invalid_argument("Cluster gap requires graphlike components with logical labels only on boundary edges");
        });
    });
}

ClusterGapGraph::ClusterGapGraph(size_t num_detectors)
    :adjacency(num_detectors + 2), relevant(num_detectors + 2, false),
     boundary_zero(num_detectors), boundary_one(num_detectors + 1)
{
}

void
ClusterGapGraph::add_edge(size_t u, size_t v, double length)
{
    if (!std::isfinite(length) || length < 0)
        throw std::invalid_argument("Cluster gap requires nonnegative finite edge weights");
    const size_t id = edges.size();
    edges.push_back({u, v, length});
    adjacency[u].push_back({v, id});
    adjacency[v].push_back({u, id});
}

void
ClusterGapGraph::finish()
{
    // Do not connect the opposing basis through the shared zero boundary.
    std::vector<size_t> stack{boundary_one};
    relevant[boundary_one] = true;
    while (!stack.empty())
    {
        const size_t u = stack.back();
        stack.pop_back();
        if (u == boundary_zero)
            continue;
        for (const auto& [v, id] : adjacency[u])
            if (!relevant[v])
            {
                relevant[v] = true;
                stack.push_back(v);
            }
    }
    if (!relevant[boundary_zero])
        throw std::invalid_argument("Cluster gap requires a path between inequivalent boundaries");
}

const std::vector<double>&
ClusterGapGraph::ball_coverage(const std::vector<std::pair<size_t, double>>& sources)
{
    distance_.assign(adjacency.size(), std::numeric_limits<double>::infinity());
    auto& distance = distance_;
    heap_.clear();
    for (const auto& [v, radius] : sources)
        if (relevant[v] && -radius < distance[v])
        {
            distance[v] = -radius;
            push(-radius, v);
        }
    while (!heap_.empty())
    {
        const auto [d, u] = pop();
        if (d != distance[u])
            continue;
        if (d >= 0)
            break;
        for (const auto& [v, id] : adjacency[u])
            if (relevant[v] && d + edges[id].length < distance[v])
            {
                distance[v] = d + edges[id].length;
                push(distance[v], v);
            }
    }
    auto& coverage = reset_coverage();
    for (size_t i = 0; i < edges.size(); i++)
    {
        const auto& edge = edges[i];
        coverage[i] = std::min(edge.length,
                std::max(0.0, -distance[edge.u]) + std::max(0.0, -distance[edge.v]));
    }
    return coverage;
}

std::vector<double>&
ClusterGapGraph::reset_coverage()
{
    coverage_.assign(edges.size(), 0.0);
    return coverage_;
}

double
ClusterGapGraph::gap(const std::vector<double>& coverage)
{
    shortest_.assign(adjacency.size(), std::numeric_limits<double>::infinity());
    auto& distance = shortest_;
    heap_.clear();
    distance[boundary_one] = 0;
    push(0, boundary_one);
    while (!heap_.empty())
    {
        const auto [d, u] = pop();
        if (d != distance[u])
            continue;
        if (u == boundary_zero)
            return d;
        for (const auto& [v, id] : adjacency[u])
        {
            if (!relevant[v])
                continue;
            const double remaining = std::max(0.0, edges[id].length - coverage[id]);
            if (d + remaining < distance[v])
            {
                distance[v] = d + remaining;
                push(distance[v], v);
            }
        }
    }
    return std::numeric_limits<double>::infinity();
}

void
ClusterGapGraph::push(double distance, size_t vertex)
{
    heap_.push_back({distance, vertex});
    std::push_heap(heap_.begin(), heap_.end(), std::greater<Item>{});
}

ClusterGapGraph::Item
ClusterGapGraph::pop()
{
    std::pop_heap(heap_.begin(), heap_.end(), std::greater<Item>{});
    const auto item = heap_.back();
    heap_.pop_back();
    return item;
}

} // namespace dec
