/*
 *  contributor: OpenAI Codex
 *  date:        13 October 2026
 *  purpose:     Implement Libra matching synthesis with PyMatching.
 * */

#include "decoder/libra.h"
#include "decoder/libra/assignment.h"
#include "globals.h"

#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <pymatching/sparse_blossom/driver/user_graph.h>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

using IdList = std::vector<uint32_t>;
using Clock = std::chrono::steady_clock;

uint64_t
elapsed_nanoseconds(Clock::time_point start)
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
}

struct DisjointSet
{
    std::vector<uint32_t> parent;

    explicit DisjointSet(size_t size)
        :parent(size)
    {
        std::iota(parent.begin(), parent.end(), 0);
    }

    uint32_t
    find(uint32_t vertex)
    {
        return parent[vertex] == vertex ? vertex : parent[vertex] = find(parent[vertex]);
    }

    void
    unite(uint32_t first, uint32_t second)
    {
        first = find(first);
        second = find(second);
        if (first != second)
            parent[second] = first;
    }
};

// Error configurations and cycles are sorted sets of channel ids. Applying a
// cycle toggles membership, so their composition is symmetric difference.
IdList
symmetric_difference(const IdList& first, const IdList& second)
{
    IdList result;
    result.reserve(first.size() + second.size());
    std::set_symmetric_difference(first.begin(), first.end(), second.begin(), second.end(),
                                 std::back_inserter(result));
    return result;
}

// Repeated occurrences cancel in pairs, as detector and error parities are binary.
IdList
cancel_duplicate_ids(IdList ids)
{
    std::sort(ids.begin(), ids.end());
    IdList result;
    for (size_t first = 0; first < ids.size();)
    {
        size_t next = first + 1;
        while (next < ids.size() && ids[next] == ids[first])
            ++next;
        if ((next - first) & 1)
            result.push_back(ids[first]);
        first = next;
    }
    return result;
}

uint64_t
edge_key(uint32_t first, uint32_t second)
{
    if (first > second)
        std::swap(first, second);
    return (uint64_t(first) << 32) | second;
}

double
log_sum_exp(double first, double second)
{
    if (first < second)
        std::swap(first, second);
    return first + std::log1p(std::exp(second - first));
}

} // anonymous namespace

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

struct Libra::Impl
{
    // A hyperedge is one independent error channel, possibly containing graph
    // edges in both CSS bases. Its weight is the log odds ln((1-p)/p).
    struct Hyperedge
    {
        IdList detectors;
        bool flips_observable;
        double probability;
        double weight;
        IdList graph_edge_ids;
    };

    struct GraphEdge
    {
        uint32_t first_detector;
        uint32_t second_detector; // UINT32_MAX represents the matching boundary.
        bool flips_observable;
        uint32_t single_channel_id = UINT32_MAX;
        IdList correlated_channel_ids;
    };

    struct Candidate
    {
        IdList error_ids;
        bool flips_observable = false;
        double weight = 0;
    };

    struct Cycle
    {
        IdList error_ids;
        double weight_delta;
        bool flips_observable = false;
    };

    struct EnsembleMember
    {
        std::unique_ptr<pm::Mwpm> matcher;
        std::vector<double> perturbed_weights;
    };

    Config config;
    size_t detector_count;
    std::vector<Hyperedge> hyperedges;
    std::vector<GraphEdge> graph_edges;
    std::unordered_map<uint64_t, uint32_t> edge_by_endpoints;
    stim::DetectorErrorModel canonical_dem;
    std::unique_ptr<pm::Mwpm> complementary_matcher;
    std::unique_ptr<pm::Mwpm> fallback_matcher;
    std::vector<EnsembleMember> ensemble_members;
    Statistics stats;
    size_t correlation_rule_count = 0;

    // Sparse scratch tables avoid clearing every detector/edge on each decode.
    // Only entries used by the current assignment or difference are reset.
    std::vector<int32_t> selected_edge_position;
    std::vector<int32_t> detector_owner;
    IdList touched_detectors;

    Impl(const stim::DetectorErrorModel& dem, Config conf)
        :config(conf),
        detector_count(dem.count_detectors())
    {
        if (dem.count_observables() != 1)
            throw std::invalid_argument("Libra requires exactly one logical observable");
        if (detector_count >= UINT32_MAX - 1)
            throw std::invalid_argument("Libra detector index overflow");
        if (!std::isfinite(conf.gap_threshold_db) || conf.gap_threshold_db < 0 ||
                conf.positive_cycle_capacity > 30)
            throw std::invalid_argument("Invalid Libra gap/cycle capacity");

        DisjointSet basis(detector_count);
        discover_css_components(dem, basis);
        const uint32_t logical_basis = find_logical_basis(dem, basis);
        build_canonical_model(dem, basis, logical_basis);
        build_matchers();

        ensemble_members.resize(config.ensemble_size);
        selected_edge_position.resize(graph_edges.size(), -1);
        detector_owner.resize(detector_count, -1);
    }

    void
    discover_css_components(const stim::DetectorErrorModel& dem, DisjointSet& basis)
    {
        dem.iter_flatten_error_instructions(
                [&] (const auto& instruction)
                {
                    instruction.for_separated_targets(
                            [&] (const auto& group)
                            {
                                IdList detectors;
                                for (auto target : group)
                                    if (target.is_relative_detector_id())
                                        detectors.push_back(target.val());
                                if (detectors.size() > 2)
                                    throw std::invalid_argument("Libra requires a decomposed DEM");
                                if (detectors.size() == 2)
                                    basis.unite(detectors[0], detectors[1]);
                            });
                });
    }

    uint32_t
    find_logical_basis(const stim::DetectorErrorModel& dem, DisjointSet& basis)
    {
        uint32_t logical_basis = UINT32_MAX;
        dem.iter_flatten_error_instructions(
                [&] (const auto& instruction)
                {
                    instruction.for_separated_targets(
                            [&] (const auto& group)
                            {
                                IdList detectors;
                                bool flips_observable = false;
                                for (auto target : group)
                                {
                                    if (target.is_relative_detector_id())
                                        detectors.push_back(target.val());
                                    else if (target.is_observable_id())
                                        flips_observable = !flips_observable;
                                }
                                if (!flips_observable)
                                    return;
                                if (detectors.size() != 1)
                                    throw std::invalid_argument(
                                        "Libra complementary matching requires observables on boundary edges");

                                const uint32_t component = basis.find(detectors[0]);
                                if (logical_basis != UINT32_MAX && logical_basis != component)
                                    throw std::invalid_argument(
                                        "Logical observable crosses multiple CSS components");
                                logical_basis = component;
                            });
                });
        if (logical_basis == UINT32_MAX)
            throw std::invalid_argument("Missing observable-crossing errors");
        return logical_basis;
    }

    void
    build_canonical_model(const stim::DetectorErrorModel& dem, DisjointSet& basis,
                          uint32_t logical_basis)
    {
        // Stim separator choices can contain redundant pieces. Recover each
        // channel's total parity, then split it by the CSS graph components.
        // Channels with identical action combine by odd-parity probability.
        std::map<std::pair<IdList, bool>, double> merged_channels;
        dem.iter_flatten_error_instructions(
                [&] (const auto& instruction)
                {
                    IdList detectors;
                    bool flips_observable = false;
                    for (auto target : instruction.target_data)
                    {
                        if (target.is_relative_detector_id())
                            detectors.push_back(target.val());
                        else if (target.is_observable_id())
                            flips_observable = !flips_observable;
                    }
                    detectors = cancel_duplicate_ids(std::move(detectors));
                    const double probability = instruction.arg_data[0];
                    if (probability == 0)
                        return;
                    if (!(probability > 0 && probability < 0.5))
                        throw std::invalid_argument("Libra requires 0 < error probabilities < 1/2");
                    if (detectors.empty())
                    {
                        if (flips_observable)
                            throw std::invalid_argument("Undetectable logical channel is unsupported");
                        return;
                    }
                    auto& merged = merged_channels[{detectors, flips_observable}];
                    merged = merged + probability - 2 * merged * probability;
                });

        for (auto& [pattern, probability] : merged_channels)
        {
            Hyperedge channel{pattern.first, pattern.second, probability,
                              std::log((1 - probability) / probability), {}};
            std::map<uint32_t, IdList> basis_components;
            for (uint32_t detector : channel.detectors)
                basis_components[basis.find(detector)].push_back(detector);
            if (basis_components.size() > 2)
                throw std::invalid_argument("Hyperedge spans more than two CSS components");

            std::vector<stim::DemTarget> targets;
            for (auto& [component, detectors] : basis_components)
            {
                if (detectors.size() > 2)
                    throw std::invalid_argument("Error has more than two detectors in one CSS basis");
                const bool flips_observable = channel.flips_observable && component == logical_basis;
                if (flips_observable && detectors.size() != 1)
                    throw std::invalid_argument("Canonical observable component is not a boundary edge");

                const uint32_t first = detectors[0];
                const uint32_t second = detectors.size() == 2 ? detectors[1] : UINT32_MAX;
                auto [position, inserted] = edge_by_endpoints.emplace(
                    edge_key(first, second), graph_edges.size());
                if (inserted)
                    graph_edges.push_back({first, second, flips_observable, UINT32_MAX, {}});
                else if (graph_edges[position->second].flips_observable != flips_observable)
                    throw std::invalid_argument("Parallel graph edges have different logical actions");
                channel.graph_edge_ids.push_back(position->second);

                if (!targets.empty())
                    targets.push_back(stim::DemTarget::separator());
                for (uint32_t detector : detectors)
                    targets.push_back(stim::DemTarget::relative_detector_id(detector));
                if (flips_observable)
                    targets.push_back(stim::DemTarget::observable_id(0));
            }
            if (channel.flips_observable && !basis_components.count(logical_basis))
                throw std::invalid_argument("Logical error has no logical-basis component");
            canonical_dem.append_error_instruction(probability, targets, "");

            const uint32_t channel_id = hyperedges.size();
            hyperedges.push_back(std::move(channel));
            if (hyperedges.back().graph_edge_ids.size() == 1)
            {
                auto& edge = graph_edges[hyperedges.back().graph_edge_ids[0]];
                if (edge.single_channel_id == UINT32_MAX ||
                        hyperedges[channel_id].weight < hyperedges[edge.single_channel_id].weight)
                    edge.single_channel_id = channel_id;
            }
            else
            {
                for (uint32_t edge_id : hyperedges.back().graph_edge_ids)
                    graph_edges[edge_id].correlated_channel_ids.push_back(channel_id);
            }
        }
        for (const auto& edge : graph_edges)
            if (edge.single_channel_id == UINT32_MAX)
                throw std::invalid_argument("A CSS component has no independent single-edge channel");
    }

    void
    build_matchers()
    {
        // Fold L0 into an explicit detector. Its unfired/fired state constrains
        // the two logical classes, while the other matching boundary remains.
        stim::DetectorErrorModel folded_dem;
        canonical_dem.iter_flatten_error_instructions(
                [&] (const auto& instruction)
                {
                    std::vector<stim::DemTarget> targets;
                    for (auto target : instruction.target_data)
                        targets.push_back(target.is_observable_id()
                            ? stim::DemTarget::relative_detector_id(detector_count) : target);
                    folded_dem.append_error_instruction(instruction.arg_data[0], targets, "");
                });
        complementary_matcher = std::make_unique<pm::Mwpm>(
            pm::detector_error_model_to_mwpm(folded_dem, pm::NUM_DISTINCT_WEIGHTS, true, true));
        fallback_matcher = std::make_unique<pm::Mwpm>(
            pm::detector_error_model_to_mwpm(canonical_dem, pm::NUM_DISTINCT_WEIGHTS, true, true));
        for (const auto& node : complementary_matcher->flooder.graph.nodes)
            for (const auto& implied_weights : node.neighbor_implied_weights)
                correlation_rule_count += implied_weights.size();
    }

    void
    prepare_member(size_t member_index)
    {
        auto& member = ensemble_members[member_index];
        if (member.matcher)
            return;
        std::mt19937_64 rng(config.seed + 0x9e3779b97f4a7c15ULL * (member_index + 1));
        std::normal_distribution<double> normal(
            0, (member_index & 1) ? std::log(4.) : std::log(2.));
        stim::DetectorErrorModel perturbed_dem;
        canonical_dem.iter_flatten_error_instructions(
                [&] (const auto& instruction)
                {
                    double probability = std::min(
                        0.499999999, instruction.arg_data[0] * std::exp(normal(rng)));
                    probability = std::max(probability, std::numeric_limits<double>::min());
                    perturbed_dem.append_error_instruction(probability, instruction.target_data, "");
                    member.perturbed_weights.push_back(std::log((1 - probability) / probability));
                });
        member.matcher = std::make_unique<pm::Mwpm>(
            pm::detector_error_model_to_mwpm(perturbed_dem, pm::NUM_DISTINCT_WEIGHTS, true, true));
    }

    Candidate
    decode_candidate(pm::Mwpm& matcher, const std::vector<uint64_t>& events, bool folded,
                     const std::vector<double>* perturbed_weights = nullptr)
    {
        std::vector<int64_t> matched_endpoints;
        pm::decode_detection_events_to_edges_with_edge_correlations(
            matcher, events, matched_endpoints);
        IdList selected_edges;
        selected_edges.reserve(matched_endpoints.size() / 2);
        for (size_t i = 0; i < matched_endpoints.size(); i += 2)
        {
            uint32_t first = matched_endpoints[i];
            uint32_t second = matched_endpoints[i + 1] == -1
                ? UINT32_MAX : uint32_t(matched_endpoints[i + 1]);
            if (folded)
            {
                if (first == detector_count)
                    first = UINT32_MAX;
                if (second == detector_count)
                    second = UINT32_MAX;
            }
            auto position = edge_by_endpoints.find(edge_key(first, second));
            if (position == edge_by_endpoints.end())
                throw std::runtime_error("Matched edge missing from original hypergraph");
            selected_edges.push_back(position->second);
        }
        selected_edges = cancel_duplicate_ids(std::move(selected_edges));
        return assign_error_channels(selected_edges, perturbed_weights);
    }

    Candidate
    assign_error_channels(const IdList& selected_edges,
                          const std::vector<double>* perturbed_weights)
    {
        auto assignment_weight = [&] (uint32_t channel_id)
        {
            return perturbed_weights ? (*perturbed_weights)[channel_id] : hyperedges[channel_id].weight;
        };
        for (size_t i = 0; i < selected_edges.size(); ++i)
            selected_edge_position[selected_edges[i]] = i;

        struct Pairing
        {
            uint32_t first_position, second_position, channel_id;
        };
        std::vector<Pairing> pairings;
        DisjointSet components(selected_edges.size());
        IdList degree(selected_edges.size());
        for (size_t i = 0; i < selected_edges.size(); ++i)
        {
            for (uint32_t channel_id : graph_edges[selected_edges[i]].correlated_channel_ids)
            {
                const auto& channel_edges = hyperedges[channel_id].graph_edge_ids;
                const uint32_t other_edge = channel_edges[0] ^ channel_edges[1] ^ selected_edges[i];
                const int32_t other_position = selected_edge_position[other_edge];
                // A pair more expensive than its independent singles cannot
                // occur in an optimal assignment. Ignore those pairings.
                if (other_position > int32_t(i) && assignment_weight(channel_id) <
                        assignment_weight(graph_edges[selected_edges[i]].single_channel_id) +
                        assignment_weight(graph_edges[other_edge].single_channel_id))
                {
                    pairings.push_back({uint32_t(i), uint32_t(other_position), channel_id});
                    components.unite(i, other_position);
                    ++degree[i];
                    ++degree[other_position];
                }
            }
        }

        Candidate result;
        std::map<uint32_t, IdList> blocks;
        for (size_t i = 0; i < selected_edges.size(); ++i)
        {
            if (degree[i])
                blocks[components.find(i)].push_back(i);
            else
                result.error_ids.push_back(graph_edges[selected_edges[i]].single_channel_id);
        }
        std::vector<int32_t> local_position(selected_edges.size(), -1);
        for (auto& [root, positions] : blocks)
        {
            const size_t count = positions.size();
            if (count == 1)
            {
                result.error_ids.push_back(graph_edges[selected_edges[positions[0]]].single_channel_id);
                continue;
            }
            if (count == 2)
            {
                // Preserve the unquantized shortcut and its strict tie rule.
                const uint32_t first_single = graph_edges[selected_edges[positions[0]]].single_channel_id;
                const uint32_t second_single = graph_edges[selected_edges[positions[1]]].single_channel_id;
                double best_cost = assignment_weight(first_single) + assignment_weight(second_single);
                uint32_t best_pair = UINT32_MAX;
                for (auto pairing : pairings)
                {
                    if (components.find(pairing.first_position) == root &&
                            assignment_weight(pairing.channel_id) < best_cost)
                    {
                        best_cost = assignment_weight(pairing.channel_id);
                        best_pair = pairing.channel_id;
                    }
                }
                if (best_pair == UINT32_MAX)
                {
                    result.error_ids.push_back(first_single);
                    result.error_ids.push_back(second_single);
                }
                else
                    result.error_ids.push_back(best_pair);
                continue;
            }

            for (size_t j = 0; j < count; ++j)
                local_position[positions[j]] = j;
            std::vector<double> single_weights;
            single_weights.reserve(count);
            for (size_t j = 0; j < count; ++j)
                single_weights.push_back(
                    assignment_weight(graph_edges[selected_edges[positions[j]]].single_channel_id));

            std::vector<libra_assignment::Pair> assignment_pairs;
            IdList pair_channels;
            for (auto pairing : pairings)
            {
                if (components.find(pairing.first_position) == root)
                {
                    assignment_pairs.push_back({uint32_t(local_position[pairing.first_position]),
                        uint32_t(local_position[pairing.second_position]),
                        assignment_weight(pairing.channel_id)});
                    pair_channels.push_back(pairing.channel_id);
                }
            }
            const auto assignment = libra_assignment::solve(single_weights, assignment_pairs);
            for (size_t j = 0; j < count; ++j)
                if (assignment.mate[j] < 0)
                    result.error_ids.push_back(graph_edges[selected_edges[positions[j]]].single_channel_id);
            for (size_t pair_index : assignment.pair_indices)
                result.error_ids.push_back(pair_channels[pair_index]);
        }
        for (uint32_t edge_id : selected_edges)
            selected_edge_position[edge_id] = -1;

        // Perturbed priors choose assignments only. All candidates are scored
        // with the original channel weights before comparing logical classes.
        result.error_ids = cancel_duplicate_ids(std::move(result.error_ids));
        for (uint32_t channel_id : result.error_ids)
        {
            result.weight += hyperedges[channel_id].weight;
            result.flips_observable ^= hyperedges[channel_id].flips_observable;
        }
        return result;
    }

    double
    relative_cycle_weight(const IdList& cycle_errors, const Candidate& base) const
    {
        double delta = 0;
        for (uint32_t channel_id : cycle_errors)
            delta += (std::binary_search(base.error_ids.begin(), base.error_ids.end(), channel_id)
                ? -1 : 1) * hyperedges[channel_id].weight;
        return delta;
    }

    std::vector<Cycle>
    decompose_difference(const Candidate& base, const Candidate& other)
    {
        IdList difference = symmetric_difference(base.error_ids, other.error_ids);
        DisjointSet components(difference.size());
        touched_detectors.clear();
        for (size_t i = 0; i < difference.size(); ++i)
        {
            for (uint32_t detector : hyperedges[difference[i]].detectors)
            {
                if (detector_owner[detector] < 0)
                {
                    detector_owner[detector] = i;
                    touched_detectors.push_back(detector);
                }
                else
                    components.unite(i, detector_owner[detector]);
            }
        }
        for (uint32_t detector : touched_detectors)
            detector_owner[detector] = -1;

        std::map<uint32_t, Cycle> cycles_by_component;
        for (size_t i = 0; i < difference.size(); ++i)
        {
            auto& cycle = cycles_by_component[components.find(i)];
            const uint32_t channel_id = difference[i];
            cycle.error_ids.push_back(channel_id);
            cycle.flips_observable ^= hyperedges[channel_id].flips_observable;
        }
        std::vector<Cycle> result;
        for (auto& [root, cycle] : cycles_by_component)
        {
            cycle.weight_delta = relative_cycle_weight(cycle.error_ids, base);
            result.push_back(std::move(cycle));
        }
        return result;
    }

    void
    retain_positive_cycle(std::vector<Cycle>& heap, Cycle cycle)
    {
        if (!config.degeneracy || !config.positive_cycle_capacity ||
                cycle.error_ids.empty() || cycle.weight_delta < -1e-9)
            return;
        cycle.weight_delta = std::max(0., cycle.weight_delta);
        for (const auto& retained : heap)
            if (retained.error_ids == cycle.error_ids)
                return;

        // A max heap discards the most expensive retained cycle when full.
        auto compare = [] (const Cycle& first, const Cycle& second)
        {
            return first.weight_delta < second.weight_delta;
        };
        if (heap.size() == config.positive_cycle_capacity)
        {
            if (cycle.weight_delta >= heap.front().weight_delta)
                return;
            std::pop_heap(heap.begin(), heap.end(), compare);
            heap.pop_back();
        }
        heap.push_back(std::move(cycle));
        std::push_heap(heap.begin(), heap.end(), compare);
    }

    void
    synthesize(Candidate& base, const Candidate& other, std::vector<Cycle>& positive_cycles)
    {
        auto cycles = decompose_difference(base, other);
        IdList errors_to_flip;
        std::vector<size_t> logical_cycles, negative_logical_cycles;
        for (size_t i = 0; i < cycles.size(); ++i)
        {
            if (cycles[i].flips_observable)
            {
                logical_cycles.push_back(i);
                if (cycles[i].weight_delta < -1e-9)
                    negative_logical_cycles.push_back(i);
            }
            else if (cycles[i].weight_delta < -1e-9)
            {
                errors_to_flip = symmetric_difference(errors_to_flip, cycles[i].error_ids);
                ++stats.improving_cycles;
            }
        }

        // An even number of logical operators preserves the representative's
        // class. If the negative subset is odd, remove its least beneficial
        // member or add the cheapest nonnegative member, whichever costs less.
        if (negative_logical_cycles.size() & 1)
        {
            size_t remove = negative_logical_cycles[0];
            for (size_t i : negative_logical_cycles)
                if (cycles[i].weight_delta > cycles[remove].weight_delta)
                    remove = i;
            size_t add = SIZE_MAX;
            for (size_t i : logical_cycles)
                if (cycles[i].weight_delta >= -1e-9 &&
                        (add == SIZE_MAX || cycles[i].weight_delta < cycles[add].weight_delta))
                    add = i;
            if (add != SIZE_MAX && cycles[add].weight_delta < -cycles[remove].weight_delta)
                negative_logical_cycles.push_back(add);
            else
                negative_logical_cycles.erase(std::find(
                    negative_logical_cycles.begin(), negative_logical_cycles.end(), remove));
        }
        double logical_delta = 0;
        for (size_t i : negative_logical_cycles)
            logical_delta += cycles[i].weight_delta;
        if (logical_delta < -1e-9)
        {
            for (size_t i : negative_logical_cycles)
                errors_to_flip = symmetric_difference(errors_to_flip, cycles[i].error_ids);
            stats.improving_cycles += negative_logical_cycles.size() / 2;
        }
        if (!errors_to_flip.empty())
        {
            base.error_ids = symmetric_difference(base.error_ids, errors_to_flip);
            base.weight = 0;
            for (uint32_t channel_id : base.error_ids)
                base.weight += hyperedges[channel_id].weight;
            positive_cycles.clear();
        }

        // The base may have changed. Recompute costs before retaining neutral
        // cycles and every pair of logical cycles for the degeneracy estimate.
        for (auto cycle : cycles)
        {
            if (!cycle.flips_observable)
            {
                cycle.weight_delta = relative_cycle_weight(cycle.error_ids, base);
                retain_positive_cycle(positive_cycles, std::move(cycle));
            }
        }
        for (size_t first = 0; first < logical_cycles.size(); ++first)
        {
            for (size_t second = first + 1; second < logical_cycles.size(); ++second)
            {
                Cycle cycle{symmetric_difference(cycles[logical_cycles[first]].error_ids,
                                                cycles[logical_cycles[second]].error_ids), 0, false};
                cycle.weight_delta = relative_cycle_weight(cycle.error_ids, base);
                retain_positive_cycle(positive_cycles, std::move(cycle));
            }
        }
    }

    double
    estimate_log_mass(const Candidate& base, const std::vector<Cycle>& cycles) const
    {
        // Only cycles sharing error channels interact in their relative weight.
        // Independent groups therefore contribute additive log partition sums.
        DisjointSet components(cycles.size());
        std::unordered_map<uint32_t, uint32_t> channel_owner;
        for (size_t i = 0; i < cycles.size(); ++i)
        {
            for (uint32_t channel_id : cycles[i].error_ids)
            {
                auto [position, inserted] = channel_owner.emplace(channel_id, i);
                if (!inserted)
                    components.unite(i, position->second);
            }
        }
        std::map<uint32_t, IdList> groups;
        for (size_t i = 0; i < cycles.size(); ++i)
            groups[components.find(i)].push_back(i);

        double log_mass = -base.weight;
        for (auto& [root, group] : groups)
        {
            std::set<IdList> span{{}};
            if (group.size() < 10)
            {
                // Enumerate the unique binary span, avoiding double counting
                // when the retained cycles are linearly dependent.
                for (uint32_t cycle_index : group)
                {
                    std::vector<IdList> additions;
                    additions.reserve(span.size());
                    for (const auto& errors : span)
                        additions.push_back(symmetric_difference(errors, cycles[cycle_index].error_ids));
                    for (auto& errors : additions)
                        span.insert(std::move(errors));
                }
            }
            else
            {
                // Bound the cost of large overlap groups by including only
                // the empty configuration, individual cycles, and cycle pairs.
                for (uint32_t cycle_index : group)
                    span.insert(cycles[cycle_index].error_ids);
                for (size_t first = 0; first < group.size(); ++first)
                    for (size_t second = first + 1; second < group.size(); ++second)
                        span.insert(symmetric_difference(cycles[group[first]].error_ids,
                                                         cycles[group[second]].error_ids));
            }
            double group_log_mass = -std::numeric_limits<double>::infinity();
            for (const auto& errors : span)
                group_log_mass = log_sum_exp(group_log_mass, -relative_cycle_weight(errors, base));
            log_mass += group_log_mass;
        }
        return log_mass;
    }

    std::pair<SyndromeType, bool>
    reconstruct_syndrome(const IdList& error_ids) const
    {
        SyndromeType syndrome(detector_count);
        bool flips_observable = false;
        for (uint32_t channel_id : error_ids)
        {
            if (channel_id >= hyperedges.size())
                throw std::out_of_range("Libra hyperedge id");
            for (uint32_t detector : hyperedges[channel_id].detectors)
                syndrome[detector] ^= 1;
            flips_observable ^= hyperedges[channel_id].flips_observable;
        }
        return {std::move(syndrome), flips_observable};
    }

    void
    validate_candidate(const Candidate& candidate, SyndromeRef given) const
    {
        auto [syndrome, flips_observable] = reconstruct_syndrome(candidate.error_ids);
        if (flips_observable != candidate.flips_observable)
            throw std::runtime_error("Libra observable invariant failed");
        for (size_t detector = 0; detector < detector_count; ++detector)
            if (syndrome[detector] != given[detector])
                throw std::runtime_error("Libra syndrome invariant failed");
    }

    DetailedResult
    decode(SyndromeRef syndrome)
    {
        ++stats.shots;
        std::vector<uint64_t> events;
        for (size_t word = 0; word < (detector_count + 63) / 64; ++word)
        {
            uint64_t bits = syndrome.u64[word];
            while (bits)
            {
                const size_t detector = word * 64 + std::countr_zero(bits);
                if (detector < detector_count)
                    events.push_back(detector);
                bits &= bits - 1;
            }
        }
        DetailedResult result;
        if (events.empty())
            return result;

        // Obtain one representative of each logical class and their gap.
        auto start = Clock::now();
        Candidate representatives[2];
        representatives[0] = decode_candidate(*complementary_matcher, events, true);
        events.push_back(detector_count);
        representatives[1] = decode_candidate(*complementary_matcher, events, true);
        events.pop_back();
        if (representatives[0].flips_observable || !representatives[1].flips_observable)
            throw std::runtime_error("Libra complementary class invariant failed");
        result.gap_db = std::abs(representatives[0].weight - representatives[1].weight) * 10 / std::log(10.);
        if (config.profile)
            stats.complementary_ns += elapsed_nanoseconds(start);

        // Section V: high-gap syndromes use the original unconstrained
        // correlated matcher. Constrained correlated matching is approximate;
        // its cheaper representative can predict a different logical class.
        if (result.gap_db >= config.gap_threshold_db)
        {
            start = Clock::now();
            auto fallback = decode_candidate(*fallback_matcher, events, false);
            result.fallback_used = true;
            result.prediction = fallback.flips_observable;
            ++stats.fallback_shots;
            representatives[fallback.flips_observable] = std::move(fallback);
            if (config.profile)
                stats.fallback_ns += elapsed_nanoseconds(start);
        }

        std::vector<Cycle> positive_cycles[2];
        if (result.gap_db < config.gap_threshold_db && config.ensemble_size)
        {
            result.ensemble_used = true;
            ++stats.ensemble_shots;
            start = Clock::now();
            std::vector<Candidate> ensemble;
            ensemble.reserve(ensemble_members.size() + 2);
            ensemble.push_back(representatives[0]);
            ensemble.push_back(representatives[1]);
            for (size_t i = 0; i < ensemble_members.size(); ++i)
            {
                prepare_member(i);
                auto candidate = decode_candidate(*ensemble_members[i].matcher, events, false,
                                                  &ensemble_members[i].perturbed_weights);
                // Consecutive identical assignments are redundant. Preserve
                // repetitions after intervening candidates: the base can change.
                if (ensemble.back().error_ids != candidate.error_ids)
                    ensemble.push_back(std::move(candidate));
            }
            if (config.profile)
                stats.ensemble_ns += elapsed_nanoseconds(start);

            start = Clock::now();
            for (size_t pass = 0; pass < (config.degeneracy ? 2 : 1); ++pass)
                for (const auto& candidate : ensemble)
                    for (size_t logical_class = 0; logical_class < 2; ++logical_class)
                        synthesize(representatives[logical_class], candidate, positive_cycles[logical_class]);
            if (config.profile)
                stats.synthesis_ns += elapsed_nanoseconds(start);
        }

        start = Clock::now();
        for (size_t logical_class = 0; logical_class < 2; ++logical_class)
        {
            result.weights[logical_class] = representatives[logical_class].weight;
            result.log_masses[logical_class] = config.degeneracy
                ? estimate_log_mass(representatives[logical_class], positive_cycles[logical_class])
                : -representatives[logical_class].weight;
            if (config.validate)
                validate_candidate(representatives[logical_class], syndrome);
        }
        if (config.profile)
            stats.degeneracy_ns += elapsed_nanoseconds(start);
        if (!result.fallback_used)
            result.prediction = result.log_masses[1] > result.log_masses[0];
        result.errors = std::move(representatives[result.prediction].error_ids);
        return result;
    }
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Libra::Libra(const stim::DetectorErrorModel& dem)
    :Libra(dem, Config{})
{}

Libra::Libra(const stim::DetectorErrorModel& dem, Config config)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    impl_(std::make_unique<Impl>(dem, config))
{}

Libra::~Libra() = default;
Libra::Libra(Libra&&) noexcept = default;

Libra::DetailedResult
Libra::decode_detailed(SyndromeRef syndrome)
{
    return impl_->decode(syndrome);
}

Result
Libra::decode(SyndromeRef syndrome, ObsRef)
{
    const auto start = Clock::now();
    auto result = decode_detailed(syndrome);
    Result out;
    out.flipped_obs = ObsType(1);
    out.flipped_obs[0] = result.prediction;
    out.matching_data.frame_flips = out.flipped_obs;
    out.matching_data.gap = result.gap_db;
    out.matching_data.probability = std::exp(result.log_masses[result.prediction]);
    const auto elapsed_ns = elapsed_nanoseconds(start);
    out.decode_latency_ns = elapsed_ns;
    auto& stats = impl_->stats;
    ++stats.timed_shots;
    stats.decode_total_ns += elapsed_ns;
    stats.decode_min_ns = std::min(stats.decode_min_ns, elapsed_ns);
    stats.decode_max_ns = std::max(stats.decode_max_ns, elapsed_ns);
    return out;
}

void
Libra::prepare_ensemble()
{
    for (size_t i = 0; i < impl_->ensemble_members.size(); ++i)
        impl_->prepare_member(i);
}

Libra::Dimensions
Libra::dimensions() const
{
    return {num_detectors, impl_->hyperedges.size(), impl_->graph_edges.size(),
            impl_->correlation_rule_count};
}

Libra::Statistics
Libra::statistics() const
{
    auto stats = impl_->stats;
    if (!stats.timed_shots)
        stats.decode_min_ns = 0;
    return stats;
}

std::pair<SyndromeType, bool>
Libra::syndrome_of(const IdList& errors) const
{
    return impl_->reconstruct_syndrome(errors);
}

void
Libra::mpi_accumulate()
{
#if defined(ENABLE_MPI)
    auto& stats = impl_->stats;
    uint64_t counters[] = {stats.shots, stats.ensemble_shots, stats.fallback_shots, stats.improving_cycles,
        stats.complementary_ns, stats.ensemble_ns, stats.synthesis_ns, stats.degeneracy_ns, stats.fallback_ns,
        stats.timed_shots, stats.decode_total_ns};
    MPI_Allreduce(MPI_IN_PLACE, counters, 11, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &stats.decode_min_ns, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &stats.decode_max_ns, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    stats = {counters[0], counters[1], counters[2], counters[3], counters[4],
             counters[5], counters[6], counters[7], counters[8], counters[9], counters[10],
             stats.decode_min_ns, stats.decode_max_ns};
#endif
}

void
Libra::print_stats(std::ostream& out) const
{
    const auto stats = statistics();
    out << "Libra shots=" << stats.shots << " ensemble_shots=" << stats.ensemble_shots
        << " fallback_shots=" << stats.fallback_shots << " improving_cycles=" << stats.improving_cycles << '\n';
    const double average_us = stats.timed_shots
        ? double(stats.decode_total_ns) / double(stats.timed_shots) / 1000.0 : 0.0;
    out << "Libra decode latency (us): timed_shots=" << stats.timed_shots
        << " average=" << average_us
        << " min=" << double(stats.decode_min_ns) / 1000.0
        << " max=" << double(stats.decode_max_ns) / 1000.0
        << " includes_empty_syndromes=true scope=full_history\n";
    if (impl_->config.profile)
        out << "Libra phase times (ns): complementary=" << stats.complementary_ns
            << " fallback=" << stats.fallback_ns << " ensemble=" << stats.ensemble_ns
            << " synthesis=" << stats.synthesis_ns << " degeneracy=" << stats.degeneracy_ns << '\n';
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec
