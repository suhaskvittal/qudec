/*
 *  contributor: OpenAI Codex
 * */

#include "decoder/uf.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <optional>
#include <stdexcept>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

// DEM targets describe parity: duplicate targets cancel in pairs.
void
cancel_duplicate_targets(std::vector<size_t>& targets)
{
    std::sort(targets.begin(), targets.end());
    size_t write{0};
    for (size_t first = 0; first < targets.size(); )
    {
        size_t next = first + 1;
        while (next < targets.size() && targets[next] == targets[first])
            next++;
        if ((next - first) % 2)
            targets[write++] = targets[first];
        first = next;
    }
    targets.resize(write);
}

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

UnionFind::UnionFind(const stim::DetectorErrorModel& dem)
    :UnionFind(dem, config_type{})
{
}

UnionFind::UnionFind(const stim::DetectorErrorModel& dem, config_type conf)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    conf(conf),
    cluster_gap_graph_(conf.enable_cluster_gap ? num_detectors : 0),
    num_vertices_(num_detectors),
    baseline_syndrome_(num_detectors),
    baseline_observables_(num_observables)
{
    if (conf.enable_cluster_gap)
        validate_cluster_gap_dem(dem);
    if (!conf.weight_scale)
        throw std::invalid_argument("UnionFind weight_scale must be positive");

    // Combine identical independent channels by their odd-parity probability.
    // Different observable actions remain separate parallel edges.
    using channel_type = std::pair<std::vector<size_t>, std::vector<size_t>>;
    std::map<channel_type, double> channels;
    dem.iter_flatten_error_instructions(
            [&] (const stim::DemInstruction& instruction)
            {
                const double p = instruction.arg_data[0];
                if (!std::isfinite(p) || p < 0 || p > 1)
                    throw std::invalid_argument("UnionFind requires probabilities in [0, 1]");
                instruction.for_separated_targets(
                        [&] (const auto& group)
                        {
                            channel_type channel;
                            for (const auto& target : group)
                            {
                                if (target.is_relative_detector_id())
                                    channel.first.push_back(target.val());
                                else if (target.is_observable_id())
                                    channel.second.push_back(target.val());
                            }
                            cancel_duplicate_targets(channel.first);
                            cancel_duplicate_targets(channel.second);
                            if (channel.first.size() > 2)
                                throw std::invalid_argument(
                                        "UnionFind requires graphlike DEM components; decompose errors first");
                            double& merged = channels[channel];
                            merged = merged * (1 - p) + (1 - merged) * p;
                        });
            });

    for (const auto& [channel, p] : channels)
    {
        if (p == 0)
            continue;
        if (conf.report_matching_data && p > 0.5)
            throw std::invalid_argument("UnionFind matched-endpoint reporting requires channel probabilities <= 1/2");

        // For p > 1/2, decode the complementary error variable with probability
        // 1-p and apply this channel as a fixed baseline correction.
        if (p > 0.5)
        {
            for (size_t d : channel.first)
                baseline_syndrome_[d] ^= true;
            for (size_t obs : channel.second)
                baseline_observables_[obs] ^= true;
        }
        if (channel.first.empty() || p == 1)
            continue;

        const double q = std::min(p, 1 - p);
        const double weight = std::log1p(-q) - std::log(q);
        const long double scaled = static_cast<long double>(weight) * conf.weight_scale;
        // Length is doubled so one unit of time is a half-edge growth step.
        if (scaled >= static_cast<long double>(std::numeric_limits<uint64_t>::max() / 2 - 1))
            throw std::overflow_error("UnionFind quantized edge weight exceeds uint64_t range");
        uint64_t length = static_cast<uint64_t>(std::floor(scaled + 0.5L));
        if (weight > 0)
            length = std::max<uint64_t>(length, 1);

        const size_t u = channel.first[0];
        // Give each boundary edge its own unconstrained terminal. Sharing a
        // terminal would merge unrelated clusters through the physical boundary.
        const size_t v = channel.first.size() == 2 ? channel.first[1] : num_vertices_++;
        edges_.push_back({u, v, 2 * length, channel.second});
        if (conf.enable_cluster_gap)
        {
            const size_t terminal = channel.second.empty() ? cluster_gap_graph_.boundary_zero
                                                           : cluster_gap_graph_.boundary_one;
            cluster_gap_edge_ids_.push_back(cluster_gap_graph_.edges.size());
            cluster_gap_graph_.add_edge(u, channel.first.size() == 2 ? v : terminal, 2 * length);
        }
    }

    if (conf.enable_cluster_gap)
        cluster_gap_graph_.finish();

    parent_.resize(num_vertices_);
    cluster_size_.resize(num_vertices_);
    parity_.resize(num_vertices_);
    boundary_.resize(num_vertices_);
    growth_.resize(edges_.size());
    forest_edges_.reserve(num_vertices_);
}

size_t
UnionFind::find_root(size_t vertex)
{
    while (parent_[vertex] != vertex)
    {
        parent_[vertex] = parent_[parent_[vertex]];
        vertex = parent_[vertex];
    }
    return vertex;
}

void
UnionFind::merge_clusters(size_t edge_id)
{
    const auto& edge = edges_[edge_id];
    size_t a = find_root(edge.u);
    size_t b = find_root(edge.v);
    if (a == b)
        return;
    if (cluster_size_[a] < cluster_size_[b])
        std::swap(a, b);
    parent_[b] = a;
    cluster_size_[a] += cluster_size_[b];
    parity_[a] ^= parity_[b];
    boundary_[a] |= boundary_[b];
    forest_edges_.push_back(edge_id);
}

void
UnionFind::grow_clusters()
{
    // Zero-weight edges form the initial erasure; merge them before growth.
    for (size_t i = 0; i < edges_.size(); i++)
        if (!edges_[i].length)
            merge_clusters(i);

    std::vector<uint8_t> rates(edges_.size());
    while (true)
    {
        bool any_active{false};
        for (size_t v = 0; v < num_vertices_; v++)
            if (parent_[v] == v && parity_[v] && !boundary_[v])
            {
                any_active = true;
                break;
            }
        if (!any_active)
            return;

        uint64_t step = std::numeric_limits<uint64_t>::max();
        for (size_t i = 0; i < edges_.size(); i++)
        {
            const auto& edge = edges_[i];
            const size_t a = find_root(edge.u);
            const size_t b = find_root(edge.v);
            rates[i] = 0;
            if (a == b)
                continue;
            const unsigned rate = (parity_[a] && !boundary_[a])
                                + (parity_[b] && !boundary_[b]);
            rates[i] = rate;
            if (rate)
            {
                const uint64_t remaining = edge.length - growth_[i];
                const uint64_t time = remaining / rate + (remaining % rate != 0);
                step = std::min(step, time);
            }
        }
        if (step == std::numeric_limits<uint64_t>::max())
            throw std::runtime_error("UnionFind syndrome has an odd component without a reachable boundary");

        // All rates are determined before any merge. Grow simultaneously, then
        // merge every saturated edge, including ties, in deterministic order.
        for (size_t i = 0; i < edges_.size(); i++)
        {
            const unsigned rate = rates[i];
            if (!rate)
                continue;
            const uint64_t remaining = edges_[i].length - growth_[i];
            const uint64_t time = remaining / rate + (remaining % rate != 0);
            if (step >= time)
                growth_[i] = edges_[i].length;
            else
                growth_[i] += step * rate;
        }
        for (size_t i = 0; i < edges_.size(); i++)
            if (rates[i] && growth_[i] == edges_[i].length)
                merge_clusters(i);
    }
}

void
UnionFind::peel_forest(std::vector<uint8_t>& syndrome, ObsRef correction,
                       std::vector<size_t>* selected_edges) const
{
    using neighbor_type = std::pair<size_t, size_t>; // Vertex, edge ID.
    std::vector<std::vector<neighbor_type>> neighbors(num_vertices_);
    for (size_t id : forest_edges_)
    {
        const auto& edge = edges_[id];
        neighbors[edge.u].push_back({edge.v, id});
        neighbors[edge.v].push_back({edge.u, id});
    }

    const size_t unseen = std::numeric_limits<size_t>::max();
    std::vector<size_t> tree_parent(num_vertices_, unseen);
    std::vector<size_t> parent_edge(num_vertices_, unseen);
    std::vector<size_t> order;
    order.reserve(num_vertices_);

    // Visit terminals first, so any component touching the boundary is rooted
    // there. The final residual parity is allowed only at such a root.
    std::vector<size_t> roots;
    roots.reserve(num_vertices_);
    for (size_t v = num_detectors; v < num_vertices_; v++)
        roots.push_back(v);
    for (size_t v = 0; v < num_detectors; v++)
        roots.push_back(v);
    for (size_t root : roots)
    {
        if (tree_parent[root] != unseen)
            continue;
        order.clear();
        tree_parent[root] = root;
        order.push_back(root);
        for (size_t j = 0; j < order.size(); j++)
        {
            const size_t v = order[j];
            for (const auto& [next, id] : neighbors[v])
                if (tree_parent[next] == unseen)
                {
                    tree_parent[next] = v;
                    parent_edge[next] = id;
                    order.push_back(next);
                }
        }
        for (size_t j = order.size(); j > 1; j--)
        {
            const size_t v = order[j - 1];
            if (!syndrome[v])
                continue;
            const auto& edge = edges_[parent_edge[v]];
            if (selected_edges)
                selected_edges->push_back(parent_edge[v]);
            for (size_t obs : edge.observables)
                correction[obs] ^= true;
            syndrome[tree_parent[v]] ^= true;
            syndrome[v] = 0;
        }
        if (root < num_detectors && syndrome[root])
            throw std::logic_error("UnionFind peeling left an unsatisfied detector");
    }
}

Result
UnionFind::decode(SyndromeRef syndrome, ObsRef /* actual_observables */)
{
    const auto start = std::chrono::steady_clock::now();
    if (syndrome.num_bits_padded() < num_detectors)
        throw std::invalid_argument("UnionFind syndrome is shorter than the DEM detector count");

    std::iota(parent_.begin(), parent_.end(), 0);
    std::fill(cluster_size_.begin(), cluster_size_.end(), 1);
    std::fill(parity_.begin(), parity_.end(), 0);
    std::fill(boundary_.begin(), boundary_.end(), 0);
    std::fill(growth_.begin(), growth_.end(), 0);
    forest_edges_.clear();

    std::vector<uint8_t> remaining(num_vertices_, 0);
    for (size_t v = 0; v < num_detectors; v++)
        remaining[v] = parity_[v] = bool(syndrome[v]) ^ bool(baseline_syndrome_[v]);
    for (size_t v = num_detectors; v < num_vertices_; v++)
        boundary_[v] = 1;

    grow_clusters();
    Result result;
    if (conf.enable_cluster_gap)
    {
        auto& coverage = cluster_gap_graph_.reset_coverage();
        for (size_t i = 0; i < edges_.size(); i++)
            coverage[cluster_gap_edge_ids_[i]] = growth_[i];
        result.cluster_gap_db = cluster_gap_graph_.gap(coverage)
                * 10.0 / (2.0 * double(conf.weight_scale) * std::log(10.0));
        result.has_cluster_gap = true;
        s_cluster_gap.add(result.cluster_gap_db);
    }
    result.flipped_obs = baseline_observables_;
    std::vector<size_t> selected_edges;
    peel_forest(remaining, result.flipped_obs, conf.report_matching_data ? &selected_edges : nullptr);
    if (conf.report_matching_data)
        report_matching(syndrome, selected_edges, result);
    result.matching_data.frame_flips = result.flipped_obs;
    result.matching_data.gap = 0;

    const auto elapsed = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - start).count());
    stats_.timed_shots++;
    stats_.decode_total_ns += elapsed;
    stats_.decode_min_ns = std::min(stats_.decode_min_ns, elapsed);
    stats_.decode_max_ns = std::max(stats_.decode_max_ns, elapsed);
    result.decode_latency_ns = elapsed;
    return result;
}

void
UnionFind::report_matching(SyndromeRef syndrome, const std::vector<size_t>& selected,
                           Result& result) const
{
    // The correction is a subforest. Pair paths at their meeting vertices,
    // carrying a single unpaired detector and its observable parity upward.
    // Different valid path decompositions can exist at branching vertices;
    // this deterministic decomposition preserves the chosen UF correction.
    using neighbor_type = std::pair<size_t, size_t>;
    struct path_type { int64_t detector; ObsType frame; };
    std::vector<std::vector<neighbor_type>> neighbors(num_vertices_);
    for (size_t id : selected)
    {
        const auto& edge = edges_[id];
        neighbors[edge.u].push_back({edge.v, id});
        neighbors[edge.v].push_back({edge.u, id});
    }
    std::vector<std::optional<path_type>> pending(num_vertices_);
    const size_t unseen = std::numeric_limits<size_t>::max();
    std::vector<size_t> parent(num_vertices_, unseen), parent_edge(num_vertices_, unseen);
    std::vector<size_t> roots, order;
    for (size_t v = num_detectors; v < num_vertices_; ++v) roots.push_back(v);
    for (size_t v = 0; v < num_detectors; ++v) roots.push_back(v);

    const auto consume = [&](size_t v, path_type path)
    {
        if (v >= num_detectors)
            result.matching_data.add({.d1=path.detector, .d2=-1, .pr=1, .w_qu=0,
                                      .frame_flips=std::move(path.frame)});
        else if (pending[v])
        {
            path.frame ^= pending[v]->frame;
            result.matching_data.add({.d1=pending[v]->detector, .d2=path.detector,
                                      .pr=1, .w_qu=0, .frame_flips=std::move(path.frame)});
            pending[v].reset();
        }
        else
            pending[v] = std::move(path);
    };
    for (size_t root : roots)
    {
        if (parent[root] != unseen) continue;
        order.clear();
        parent[root] = root;
        order.push_back(root);
        for (size_t i = 0; i < order.size(); ++i)
            for (const auto& [next, id] : neighbors[order[i]])
                if (parent[next] == unseen)
                {
                    parent[next] = order[i];
                    parent_edge[next] = id;
                    order.push_back(next);
                }
        for (size_t i = order.size(); i > 0; --i)
        {
            const size_t v = order[i - 1];
            if (v < num_detectors && syndrome[v])
                consume(v, {static_cast<int64_t>(v), ObsType(num_observables)});
            if (!pending[v]) continue;
            if (v == root)
                throw std::logic_error("UnionFind path decomposition left an unmatched detector");
            auto path = std::move(*pending[v]);
            pending[v].reset();
            for (size_t obs : edges_[parent_edge[v]].observables)
                path.frame[obs] ^= true;
            consume(parent[v], std::move(path));
        }
    }
}

UnionFind::Statistics
UnionFind::statistics() const
{
    auto stats = stats_;
    if (!stats.timed_shots)
        stats.decode_min_ns = 0;
    return stats;
}

void
UnionFind::print_stats(std::ostream& ostrm) const
{
    const auto stats = statistics();
    const double average_us = stats.timed_shots
        ? double(stats.decode_total_ns) / double(stats.timed_shots) / 1000.0 : 0.0;
    ostrm << "UnionFind decode latency (us): timed_shots=" << stats.timed_shots
          << " average=" << average_us
          << " min=" << double(stats.decode_min_ns) / 1000.0
          << " max=" << double(stats.decode_max_ns) / 1000.0
          << " includes_empty_syndromes=true scope=full_history\n";
    if (conf.enable_cluster_gap)
    {
        ostrm << "Cluster gap (dB): count=" << stats.timed_shots;
        if (stats.timed_shots)
            ostrm << " " << s_cluster_gap.to_string_some();
        ostrm << "\n";
        if (conf.print_cluster_gap_histogram)
            ostrm << s_cluster_gap.to_string_full() << "\n";
    }
}

void
UnionFind::mpi_accumulate()
{
    s_cluster_gap.mpi_accumulate();
#if defined(ENABLE_MPI)
    uint64_t sums[]{stats_.timed_shots, stats_.decode_total_ns};
    MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    stats_.timed_shots = sums[0];
    stats_.decode_total_ns = sums[1];
    MPI_Allreduce(MPI_IN_PLACE, &stats_.decode_min_ns, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &stats_.decode_max_ns, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
#endif
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec
