/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Independent implementation of Lee--Li--Bartlett concatenated MWPM.
 *  Optional first-stage randomized ensemble, coalesced by virtual syndrome.
 *  Adapted from the project implementation by OpenAI Codex, 24 October 2026.
 * */
#include "decoder/concat_mwpm/matching.h"
#include <random>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <ostream>

namespace dec::concat
{
namespace
{
using Edge = std::pair<int64_t, int64_t>;

double
xor_probability(double p, double q)
{
    return p + q - 2*p*q;
}

double
error_weight(double p)
{
    if (!(p > 0 && p <= 1)) throw std::invalid_argument("Need 0 < event probability <= 1");
    return -std::log(p);
}

Edge
edge_key(const std::vector<size_t>& support)
{
    if (support.empty() || support.size() > 2) throw std::invalid_argument("Not graphlike");
    return std::minmax(static_cast<int64_t>(support[0]),
                      support.size() == 1 ? int64_t{-1} : static_cast<int64_t>(support[1]));
}

Model
compress(Model model, bool preserve_order)
{
    if (model.colors.empty()) throw std::invalid_argument("Concatenated MWPM needs check detectors");
    std::map<std::pair<std::vector<size_t>, uint64_t>, double> merged;
    std::map<std::pair<std::vector<size_t>, uint64_t>, BoundaryHint> boundary_kind;
    std::map<std::pair<std::vector<size_t>, uint64_t>, bool> supported;
    std::vector<std::pair<std::vector<size_t>, uint64_t>> insertion_order;
    for (auto& e : model.events)
    {
        std::sort(e.detectors.begin(), e.detectors.end());
        if (std::adjacent_find(e.detectors.begin(), e.detectors.end()) != e.detectors.end())
            throw std::invalid_argument("Repeated detector within event");
        for (auto v : e.detectors)
            if (v >= model.colors.size()) throw std::out_of_range("Detector outside model");
        if (e.detectors.empty()) throw std::invalid_argument("Unobservable event not supported");
        error_weight(e.probability);
        const auto key = std::pair{e.detectors,e.observable};
        if (!merged.contains(key)) { boundary_kind[key] = e.boundary_hint; insertion_order.push_back(key); }
        else if (boundary_kind[key] != e.boundary_hint)
            boundary_kind[key] = BoundaryHint::Unknown;
        auto& p = merged[key];
        // Author: OpenAI Codex. C_f probabilities contain supported single
        // faults, rather than also adding suppressed unsupported channels.
        if (model.flag_renormalization && e.flag_supported && !supported[key]) p = 0;
        supported[key] = supported[key] || e.flag_supported;
        if (!model.flag_renormalization || e.flag_supported || !supported[key])
            p = xor_probability(p, e.probability);
    }
    model.events.clear();
    // Contributor: OpenAI Codex. In the opt-in boundary ensemble retain input
    // event order, including PyMatching tie choices, as in the validated model.
    if (preserve_order)
        for (const auto& key : insertion_order)
            model.events.push_back({merged.at(key),key.first,key.second,boundary_kind.at(key),supported.at(key)});
    else
        for (const auto& [key,p] : merged)
            model.events.push_back({p,key.first,key.second,boundary_kind.at(key),supported.at(key)});
    for (int c : model.colors)
        if (c < 0 || c > 2) throw std::invalid_argument("Invalid detector color");
    return model;
}
} // namespace

struct Decoder::Graph
{
    pm::Mwpm matcher;
    std::map<Edge, size_t> records;
    std::vector<double> probabilities;
    std::vector<Edge> endpoints;
    std::vector<size_t> component;
    std::vector<bool> component_boundary;

    Graph(size_t num_detectors, const std::vector<Event>& events, bool uniform)
    {
        stim::DetectorErrorModel dem;
        for (size_t i = 0; i < events.size(); i++)
        {
            const auto& e = events[i];
            probabilities.push_back(e.probability > 0 ? (uniform ? .1 : e.probability) : 0);
            endpoints.push_back(edge_key(e.detectors));
            if (e.probability == 0) continue;
            if (!records.emplace(edge_key(e.detectors), i).second)
                throw std::invalid_argument("Parallel graph edges need explicit treatment");
            std::vector<stim::DemTarget> targets;
            for (auto v : e.detectors) targets.push_back(stim::DemTarget::relative_detector_id(v));
            // Contributor: OpenAI Codex. PyMatching converts DEM probabilities
            // to log odds. Encode q so that log((1-r)/r) == -log(q).
            const double q = uniform ? .1 : e.probability;
            dem.append_error_instruction(q/(1+q),targets, "");
        }
        dem += stim::DetectorErrorModel("detector D" + std::to_string(num_detectors-1));
        matcher = pm::detector_error_model_to_mwpm(dem, 1 << 24, true, false);
        // Contributor: OpenAI Codex. Detect infeasible boundary-parity sectors without
        // invoking/resetting PyMatching on an unmatched odd component.
        component.resize(num_detectors);
        for (size_t v = 0; v < num_detectors; v++) component[v] = v;
        auto root = [&](size_t v)
        {
            while (component[v] != v) v = component[v];
            return v;
        };
        for (const auto& [edge,id] : records)
            if (edge.first >= 0)
            {
                const size_t a = root(edge.first), b = root(edge.second);
                component[b] = a;
            }
        for (size_t v = 0; v < num_detectors; v++) component[v] = root(v);
        component_boundary.resize(num_detectors, false);
        for (const auto& [edge,id] : records)
            if (edge.first < 0) component_boundary[component[edge.second]] = true;
    }

    bool can_solve(const std::vector<uint64_t>& fired) const
    {
        std::vector<bool> parity(component.size(), false);
        for (auto v : fired) parity.at(component.at(v)) = !parity[component[v]];
        for (size_t v = 0; v < parity.size(); v++)
            if (parity[v] && !component_boundary[v]) return false;
        return true;
    }

    size_t audit_weights() const
    {
        const auto& graph = matcher.flooder.graph;
        const auto& search = matcher.search_flooder.graph;
        const double scale = graph.normalising_constant;
        for (const auto& [edge, id] : records)
        {
            const int64_t x = edge.first < 0 ? edge.second : edge.first;
            const int64_t y = edge.first < 0 ? -1 : edge.second;
            const auto& node = graph.nodes[x];
            const auto& search_node = search.nodes[x];
            const auto matching_it = std::find(node.neighbors.begin(),node.neighbors.end(),
                y < 0 ? nullptr : &graph.nodes[y]);
            const auto search_it = std::find(search_node.neighbors.begin(),search_node.neighbors.end(),
                y < 0 ? nullptr : &search.nodes[y]);
            if (matching_it == node.neighbors.end() || search_it == search_node.neighbors.end())
                throw std::runtime_error("Reweighted graph edge was lost");
            const double expected = error_weight(probabilities[id]);
            const double actual = node.neighbor_weights[matching_it-node.neighbors.begin()]/scale;
            const double actual_search = search_node.neighbor_weights[search_it-search_node.neighbors.begin()]/scale;
            if (std::abs(actual-expected) > 2/scale || std::abs(actual_search-expected) > 2/scale)
                throw std::runtime_error("Actual matching/search weights disagree with negative-log probabilities");
        }
        return records.size();
    }

    std::vector<size_t> solve(const std::vector<uint64_t>& fired,
                             std::vector<std::vector<size_t>>* paths = nullptr,
                             std::vector<std::vector<int64_t>>* node_paths = nullptr)
    {
        std::vector<int64_t> edges;
        pm::decode_detection_events_to_edges(matcher, fired, edges);
        std::map<size_t,bool> traced_parity;
        if (paths || node_paths)
        {
            const auto* nodes = matcher.flooder.graph.nodes.data();
            const auto* search_nodes = matcher.search_flooder.graph.nodes.data();
            for (const auto& match : matcher.flooder.match_edges)
            {
                std::vector<size_t> path, unordered;
                size_t from = match.loc_from-nodes;
                size_t to = match.loc_to ? match.loc_to-nodes : SIZE_MAX;
                matcher.search_flooder.iter_edges_on_shortest_path_from_middle(from, to,
                    [&](const pm::SearchGraphEdge& e)
                    {
                        auto* u = e.detector_node;
                        auto* v = u->neighbors[e.neighbor_index];
                        const auto id = records.at(std::minmax(int64_t(u-search_nodes),
                            v ? int64_t(v-search_nodes) : int64_t{-1}));
                        unordered.push_back(id);
                        traced_parity[id] = !traced_parity[id];
                    });
                int64_t current = from;
                std::vector<int64_t> nodes_path{current};
                while (!unordered.empty())
                {
                    const auto it = std::find_if(unordered.begin(),unordered.end(),[&](auto id)
                    { return endpoints[id].first == current || endpoints[id].second == current; });
                    if (it == unordered.end()) throw std::runtime_error("Disconnected concatenated matching trace");
                    const auto id = *it;
                    unordered.erase(it);
                    path.push_back(id);
                    current = endpoints[id].first == current ? endpoints[id].second : endpoints[id].first;
                    nodes_path.push_back(current);
                    if (current < 0 && !unordered.empty()) throw std::runtime_error("Matching trace exits a boundary");
                }
                if (current != (to == SIZE_MAX ? int64_t{-1} : int64_t(to)))
                    throw std::runtime_error("Wrong concatenated matching trace endpoint");
                if (paths) paths->push_back(std::move(path));
                if (node_paths) node_paths->push_back(std::move(nodes_path));
            }
        }
        std::vector<size_t> ids;
        for (size_t i = 0; i < edges.size(); i += 2)
            ids.push_back(records.at(std::minmax(edges[i], edges[i+1])));
        if (paths || node_paths)
        {
            for (auto id : ids) traced_parity[id] = !traced_parity[id];
            for (const auto& [id,odd] : traced_parity)
                if (odd) throw std::runtime_error("Concatenated matching trace differs from selected events");
        }
        return ids;
    }
};

struct Decoder::SubDecoder
{
    std::shared_ptr<Graph> restricted;
    std::shared_ptr<Graph> boundary_graph;
    std::unique_ptr<Graph> monochromatic;
    std::vector<Event> second_events;
    std::vector<Event> first_events;
    std::vector<std::shared_ptr<Graph>> ensemble;
    std::vector<size_t> original_events;
};

size_t
Decoder::ambiguous_boundary_events() const
{
    size_t count = 0;
    for (const auto& e : model.events)
        if (e.boundary_hint == BoundaryHint::Unknown)
            for (int c = 0; c < 3; c++)
            {
                size_t other = 0;
                for (auto v : e.detectors) other += model.colors[v] != c;
                if (other == 1) { count++; break; }
            }
    return count;
}

size_t
Decoder::audit_weights() const
{
    size_t count = 0;
    for (const auto& sub : subdecoders_)
    {
        for (const auto& graph : sub->ensemble) count += graph->audit_weights();
        count += sub->monochromatic->audit_weights();
    }
    return count;
}

Decoder::Decoder(Model input, bool uniform, bool boundary_gap, EnsembleConfig config)
    :model(compress(std::move(input),boundary_gap)), boundary_gap_(boundary_gap), config_(config)
{
    if (config.size == 0 || !(config.spread_scale >= 0) || !std::isfinite(config.spread_scale))
        throw std::invalid_argument("Invalid ensemble configuration");
    const size_t n = model.colors.size();
    for (int c = 0; c < 3; c++)
    {
        auto sub = std::make_shared<SubDecoder>();
        // Contributor: OpenAI Codex. Complete known spatial hyperedges with
        // their missing-color boundaries before restriction. Boundary vertices
        // are n+color; provenance, never the logical label, determines this.
        auto projection = [&](const Event& e)
        {
            std::vector<size_t> other;
            std::array<bool,3> present{};
            for (auto v : e.detectors)
            {
                if (boundary_gap_ && e.boundary_hint == BoundaryHint::Spatial && present[model.colors[v]])
                    throw std::invalid_argument("Spatial boundary provenance requires at most one detector per color");
                present[model.colors[v]] = true;
                if (model.colors[v] != c) other.push_back(v);
            }
            if (boundary_gap_ && e.boundary_hint == BoundaryHint::Spatial)
                for (int color = 0; color < 3; color++)
                    if (color != c && !present[color]) other.push_back(n+color);
            std::sort(other.begin(),other.end());
            if (other.size() > 2) throw std::invalid_argument("Non-graphlike restricted projection");
            return other;
        };
        // Compress first-stage projections, including physical positive-weight
        // boundary--boundary edges from spatial corner data faults.
        std::map<std::vector<size_t>,double> projected;
        std::map<std::vector<size_t>,double> flagged_projected;
        for (const auto& e : model.events)
        {
            const auto other = projection(e);
            if (!other.empty())
            {
                auto& p = projected[other];
                p = xor_probability(p,e.probability);
                if (model.flag_renormalization && e.flag_supported)
                {
                    auto& flagged = flagged_projected[other];
                    flagged = xor_probability(flagged,e.probability);
                }
            }
        }
        std::map<std::vector<size_t>,size_t> virtual_id;
        std::vector<Event> first_events;
        for (const auto& [support,p] : projected)
        {
            virtual_id[support] = n+first_events.size();
            first_events.push_back({flagged_projected.contains(support) ? flagged_projected.at(support) : p,support,0});
        }
        sub->first_events = first_events;
        if (boundary_gap_)
            sub->boundary_graph = std::make_shared<Graph>(n+3,first_events,uniform);
        else
            sub->restricted = std::make_shared<Graph>(n,first_events,uniform);
        sub->ensemble.push_back(boundary_gap_ ? sub->boundary_graph : sub->restricted);
        // Fixed seeded graphs. Alternating spreads make ensemble prefixes nested.
        // Perturb probabilities, not final correction scores; only stage one changes.
        std::mt19937_64 rng(config.seed+0x9e3779b97f4a7c15ULL*(c+1));
        for (size_t member = 1; member < config.size; member++)
        {
            const double sigma = std::log(member%2 ? 2.0 : 4.0)*config.spread_scale;
            std::normal_distribution<double> normal(0,1);
            auto perturbed = first_events;
            for (auto& event : perturbed)
                event.probability = std::clamp((uniform ? .1 : event.probability)*std::exp(sigma*normal(rng)),
                                               1e-12, .499999999);
            sub->ensemble.push_back(std::make_shared<Graph>(boundary_gap_ ? n+3 : n,perturbed,false));
        }
        // Algorithm 1, lines 19--28: keep the c-colored detectors and
        // replace the entire other-color support by its virtual detector.
        // Pure c-colored measurement faults keep both temporal endpoints.
        for (size_t id = 0; id < model.events.size(); id++)
        {
            const auto& e = model.events[id];
            std::vector<size_t> own, other;
            for (auto v : e.detectors) (model.colors[v] == c ? own : other).push_back(v);
            // Preserve the original second stage: only real other-color
            // support creates a virtual check. Boundary-only projections are
            // unobserved, even though their parity was fixed in stage one.
            if (!other.empty()) own.push_back(virtual_id.at(projection(e)));
            if (own.size() > 2) throw std::invalid_argument("Non-graphlike monochromatic projection");
            sub->second_events.push_back({e.probability, own, e.observable});
            sub->original_events.push_back(id);
        }
        sub->monochromatic = std::make_unique<Graph>(n+first_events.size(), sub->second_events, false);
        subdecoders_[c] = std::move(sub);
    }
}

Result
Decoder::decode(const std::vector<uint8_t>& syndrome)
{
    auto candidates = decode_candidates(syndrome, false);
    Result best{0, std::numeric_limits<double>::infinity(), -1, {}, {}, {}};
    for (auto& candidate : candidates)
        if (candidate.weight < best.weight - 1e-9) best = std::move(candidate);
    return best;
}

std::vector<Result>
Decoder::decode_candidates(const std::vector<uint8_t>& syndrome, bool trace, DebugTrace* debug)
{
    if (syndrome.size() != model.colors.size()) throw std::invalid_argument("Wrong syndrome size");
    statistics = {};
    std::vector<Result> out;
    std::array<std::map<std::vector<uint64_t>,bool>,3> seen;
    for (size_t member = 0; member < config_.size; member++)
      for (int c = 0; c < 3; c++)
        for (int mask = boundary_gap_ ? 0 : 3; mask < 4; mask++)
        {
            auto candidate = decode_color(syndrome,c,mask,member,trace,seen[c],debug);
            if (candidate) out.push_back(std::move(*candidate));
        }
    if (out.empty()) throw std::runtime_error("No feasible concatenated matching candidate");
    return out;
}

std::array<Result, 3>
Decoder::decode_colors(const std::vector<uint8_t>& syndrome, bool trace)
{
    auto candidates = decode_candidates(syndrome,trace);
    std::array<Result,3> out;
    for (int c = 0; c < 3; c++)
        out[c] = {0,std::numeric_limits<double>::infinity(),c,{},{},{}};
    for (auto& candidate : candidates)
        if (candidate.weight < out[candidate.color].weight-1e-9)
            out[candidate.color] = std::move(candidate);
    return out;
}

std::optional<Result>
Decoder::decode_color(const std::vector<uint8_t>& syndrome, int c, int mask, size_t member, bool trace,
                      std::map<std::vector<uint64_t>,bool>& seen, DebugTrace* debug)
{
    auto& sub = *subdecoders_[c];
    std::vector<uint64_t> first_fired, second_fired;
    for (size_t v = 0; v < syndrome.size(); v++)
        if (syndrome[v]) (model.colors[v] == c ? second_fired : first_fired).push_back(v);
    Result candidate{0,0,c,{},{},{}};
    candidate.boundary_mask = mask;
    candidate.member = member;
    auto& graph = sub.ensemble[member];
    if (boundary_gap_)
        for (int bit = 0; bit < 2; bit++)
            if (mask & (1 << bit)) first_fired.push_back(model.colors.size()+(c+1+bit)%3);
    // Contributor: OpenAI Codex. Trace the actual matching calls without replay.
    auto node_name = [&](int64_t id,bool first_stage)
    {
        if (id < 0) return std::string("BoundarySink");
        if (static_cast<size_t>(id) < model.colors.size())
            return debug && static_cast<size_t>(id) < debug->detector_labels.size()
                ? debug->detector_labels[id] : "v"+std::to_string(id);
        if (first_stage) return std::string(1,"RGB"[id-model.colors.size()])+"Boundary";
        return "V"+std::to_string(id-model.colors.size());
    };
    auto print_matches = [&](const char* stage,const auto& node_paths,const auto& edge_paths,
                             const Graph& matched_graph,bool first_stage)
    {
        if (!debug) return;
        debug->out << "  " << stage << " matching:\n";
        if (node_paths.empty()) debug->out << "    empty matching\n";
        for (size_t j = 0; j < node_paths.size(); j++)
        {
            const auto& nodes = node_paths[j];
            debug->out << "    match " << node_name(nodes.front(),first_stage) << " <-> "
                       << node_name(nodes.back(),first_stage) << " path: " << node_name(nodes.front(),first_stage);
            for (size_t k = 0; k < edge_paths[j].size(); k++)
            {
                const auto id = edge_paths[j][k];
                const double p = matched_graph.probabilities[id];
                debug->out << " --{" << (first_stage ? "P" : "E")
                    << (first_stage ? id : sub.original_events[id]) << ", p=" << p << ", W=" << error_weight(p);
                if (first_stage) debug->out << ", base_p=" << sub.first_events[id].probability;
                debug->out << "}--> " << node_name(nodes[k+1],first_stage);
            }
            debug->out << '\n';
        }
    };
    if (debug)
    {
        debug->out << "Candidate !" << "RGB"[c] << " (" << "RGB"[(c+1)%3] << "RGB"[(c+2)%3]
                   << ") member=" << member << " boundary_mask=" << mask
                   << (boundary_gap_ ? " (parity sector)" : " (ordinary boundaries)") << '\n';
        debug->out << "  Stage 1 input:";
        for (auto v : first_fired) debug->out << ' ' << node_name(v,true);
        debug->out << '\n';
    }
    if (!graph->can_solve(first_fired))
    {
        if (debug) debug->out << "  Stage 1 infeasible; no matching call\n";
        statistics.infeasible++; return std::nullopt;
    }
    statistics.first_calls++;
    std::vector<std::vector<size_t>> first_paths;
    std::vector<std::vector<int64_t>> first_nodes;
    const auto first_selected = graph->solve(first_fired,debug ? &first_paths : nullptr,debug ? &first_nodes : nullptr);
    print_matches("Stage 1",first_nodes,first_paths,*graph,true);
    for (auto id : first_selected)
    {
        const auto& support = sub.first_events[id].detectors;
        if (std::any_of(support.begin(),support.end(),[&](auto v) { return v < model.colors.size(); }))
            second_fired.push_back(model.colors.size()+id);
        if (trace) candidate.restricted_edges.push_back(support);
    }
    auto canonical = second_fired;
    std::sort(canonical.begin(),canonical.end());
    // Preserve historical K=1 boundary-sector candidates and tie choices.
    if (config_.size > 1 && !seen.emplace(canonical,true).second)
    {
        if (debug) debug->out << "  Coalesced: identical second-stage syndrome already decoded\n";
        statistics.duplicates++; return std::nullopt;
    }
    for (auto id : canonical) if (id >= model.colors.size()) candidate.virtual_detectors.push_back(id);
    if (debug)
    {
        debug->out << "  Virtual detector definitions (selected first-stage edges):\n";
        for (auto id : first_selected)
        {
            debug->out << "    V" << id << " = {";
            for (auto v : sub.first_events[id].detectors) debug->out << node_name(v,true) << ' ';
            debug->out << "}\n";
        }
        debug->out << "  Stage 2 input:";
        for (auto v : second_fired) debug->out << ' ' << node_name(v,false);
        debug->out << '\n';
    }
    if (!sub.monochromatic->can_solve(second_fired))
    {
        if (debug) debug->out << "  Stage 2 infeasible; candidate discarded\n";
        return std::nullopt;
    }
    statistics.second_calls++;
    auto residual = syndrome;
    std::vector<std::vector<size_t>> paths;
    std::vector<std::vector<int64_t>> second_nodes;
    const auto second_selected = sub.monochromatic->solve(second_fired,(trace || debug) ? &paths : nullptr,
                                                        debug ? &second_nodes : nullptr);
    print_matches("Stage 2",second_nodes,paths,*sub.monochromatic,false);
    for (auto id : second_selected)
    {
        const auto original_id = sub.original_events[id];
        const auto& e = model.events[original_id];
        candidate.events.push_back(original_id);
        candidate.observable ^= e.observable;
        candidate.weight += error_weight(e.probability);
        for (auto v : e.detectors) residual[v] ^= 1;
    }
    for (const auto& path : paths)
    {
        std::vector<size_t> original_path;
        for (auto id : path) original_path.push_back(sub.original_events[id]);
        candidate.second_stage_paths.push_back(std::move(original_path));
    }
    if (std::any_of(residual.begin(),residual.end(),[](auto x) { return x != 0; }))
        throw std::runtime_error("Concatenated matching did not reproduce full syndrome");
    if (debug)
    {
        debug->out << "  Correction: weight=" << candidate.weight << " observable=" << candidate.observable << " events:";
        for (auto id : candidate.events) debug->out << " E" << id;
        debug->out << '\n';
    }
    return candidate;
}
} // namespace dec::concat
