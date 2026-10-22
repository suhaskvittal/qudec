/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Independent implementation of Lee--Li--Bartlett concatenated MWPM.
 *  No restriction lifting, exact-decoder fallback, or randomized ensemble.
 * */
#include "decoder/concat_mwpm/matching.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

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
log_odds(double p)
{
    if (!(p > 0 && p < 0.5)) throw std::invalid_argument("Need 0 < event probability < 0.5");
    return std::log((1-p)/p);
}

Edge
edge_key(const std::vector<size_t>& support)
{
    if (support.empty() || support.size() > 2) throw std::invalid_argument("Not graphlike");
    return std::minmax(static_cast<int64_t>(support[0]),
                      support.size() == 1 ? int64_t{-1} : static_cast<int64_t>(support[1]));
}

Model
compress(Model model)
{
    if (model.colors.empty()) throw std::invalid_argument("Concatenated MWPM needs check detectors");
    std::map<std::pair<std::vector<size_t>, uint64_t>, double> merged;
    for (auto& e : model.events)
    {
        std::sort(e.detectors.begin(), e.detectors.end());
        if (std::adjacent_find(e.detectors.begin(), e.detectors.end()) != e.detectors.end())
            throw std::invalid_argument("Repeated detector within event");
        for (auto v : e.detectors)
            if (v >= model.colors.size()) throw std::out_of_range("Detector outside model");
        if (e.detectors.empty()) throw std::invalid_argument("Unobservable event not supported");
        log_odds(e.probability);
        auto& p = merged[{e.detectors, e.observable}];
        p = xor_probability(p, e.probability);
    }
    model.events.clear();
    for (const auto& [key, p] : merged) model.events.push_back({p, key.first, key.second});
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

    Graph(size_t num_detectors, const std::vector<Event>& events, bool uniform)
    {
        stim::DetectorErrorModel dem;
        for (size_t i = 0; i < events.size(); i++)
        {
            const auto& e = events[i];
            probabilities.push_back(uniform ? .1 : e.probability);
            if (!records.emplace(edge_key(e.detectors), i).second)
                throw std::invalid_argument("Parallel graph edges need explicit treatment");
            std::vector<stim::DemTarget> targets;
            for (auto v : e.detectors) targets.push_back(stim::DemTarget::relative_detector_id(v));
            dem.append_error_instruction(uniform ? 0.1 : e.probability, targets, "");
        }
        dem += stim::DetectorErrorModel("detector D" + std::to_string(num_detectors-1));
        matcher = pm::detector_error_model_to_mwpm(dem, 1 << 24, true, false);
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
            const double expected = log_odds(probabilities[id]);
            const double actual = node.neighbor_weights[matching_it-node.neighbors.begin()]/scale;
            const double actual_search = search_node.neighbor_weights[search_it-search_node.neighbors.begin()]/scale;
            if (std::abs(actual-expected) > 2/scale || std::abs(actual_search-expected) > 2/scale)
                throw std::runtime_error("Actual matching/search weights disagree with reweighted log odds");
        }
        return records.size();
    }

    std::vector<size_t> solve(const std::vector<uint64_t>& fired,
                             std::vector<std::vector<size_t>>* paths = nullptr)
    {
        std::vector<int64_t> edges;
        pm::decode_detection_events_to_edges(matcher, fired, edges);
        if (paths)
        {
            const auto* nodes = matcher.flooder.graph.nodes.data();
            const auto* search_nodes = matcher.search_flooder.graph.nodes.data();
            for (const auto& match : matcher.flooder.match_edges)
            {
                std::vector<size_t> path;
                size_t from = match.loc_from-nodes;
                size_t to = match.loc_to ? match.loc_to-nodes : SIZE_MAX;
                matcher.search_flooder.iter_edges_on_shortest_path_from_source(from, to,
                    [&](const pm::SearchGraphEdge& e)
                    {
                        auto* u = e.detector_node;
                        auto* v = u->neighbors[e.neighbor_index];
                        path.push_back(records.at(std::minmax(int64_t(u-search_nodes),
                            v ? int64_t(v-search_nodes) : int64_t{-1})));
                    });
                paths->push_back(std::move(path));
            }
        }
        std::vector<size_t> ids;
        for (size_t i = 0; i < edges.size(); i += 2)
            ids.push_back(records.at(std::minmax(edges[i], edges[i+1])));
        return ids;
    }
};

struct Decoder::SubDecoder
{
    std::unique_ptr<Graph> restricted, monochromatic;
    std::vector<Event> second_events;
    std::vector<Event> first_events;
    std::vector<size_t> original_events;
};

size_t
Decoder::audit_weights() const
{
    size_t count = 0;
    for (const auto& sub : subdecoders_)
        count += sub->restricted->audit_weights()+sub->monochromatic->audit_weights();
    return count;
}

Decoder::Decoder(Model input, bool uniform)
    :model(compress(std::move(input)))
{
    const size_t n = model.colors.size();
    for (int c = 0; c < 3; c++)
    {
        auto sub = std::make_shared<SubDecoder>();
        // Algorithm 1, lines 11--18: compress projections, create one
        // virtual detector per distinct nonempty restricted edge support.
        std::map<std::vector<size_t>, double> projected;
        for (const auto& e : model.events)
        {
            std::vector<size_t> other;
            for (auto v : e.detectors) if (model.colors[v] != c) other.push_back(v);
            if (other.size() > 2) throw std::invalid_argument("Non-graphlike restricted projection");
            if (!other.empty())
            {
                auto& p = projected[other];
                p = xor_probability(p, e.probability);
            }
        }
        std::map<std::vector<size_t>, size_t> virtual_id;
        std::vector<Event> first_events;
        for (const auto& [support, p] : projected)
        {
            virtual_id[support] = n + first_events.size();
            first_events.push_back({p, support, 0});
        }
        sub->restricted = std::make_unique<Graph>(n, first_events, uniform);
        sub->first_events = first_events;
        // Algorithm 1, lines 19--28: keep the c-colored detectors and
        // replace the entire other-color support by its virtual detector.
        // Pure c-colored measurement faults keep both temporal endpoints.
        for (size_t id = 0; id < model.events.size(); id++)
        {
            const auto& e = model.events[id];
            std::vector<size_t> own, other;
            for (auto v : e.detectors) (model.colors[v] == c ? own : other).push_back(v);
            if (!other.empty()) own.push_back(virtual_id.at(other));
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
    auto candidates = decode_colors(syndrome, false);
    Result best{0, std::numeric_limits<double>::infinity(), -1, {}, {}, {}};
    for (auto& candidate : candidates)
        if (candidate.weight < best.weight - 1e-9) best = std::move(candidate);
    return best;
}

std::array<Result, 3>
Decoder::decode_colors(const std::vector<uint8_t>& syndrome, bool trace)
{
    if (syndrome.size() != model.colors.size()) throw std::invalid_argument("Wrong syndrome size");
    std::array<Result, 3> out;
    for (int c = 0; c < 3; c++)
    {
        auto& sub = *subdecoders_[c];
        std::vector<uint64_t> first_fired, second_fired;
        for (size_t v = 0; v < syndrome.size(); v++)
            if (syndrome[v]) (model.colors[v] == c ? second_fired : first_fired).push_back(v);
        Result candidate{0, 0, c, {}, {}, {}};
        for (auto id : sub.restricted->solve(first_fired))
        {
            second_fired.push_back(model.colors.size()+id);
            if (trace) candidate.restricted_edges.push_back(sub.first_events[id].detectors);
        }
        auto residual = syndrome;
        std::vector<std::vector<size_t>> paths;
        for (auto id : sub.monochromatic->solve(second_fired, trace ? &paths : nullptr))
        {
            const auto original_id = sub.original_events[id];
            const auto& e = model.events[original_id];
            candidate.events.push_back(original_id);
            candidate.observable ^= e.observable;
            candidate.weight += log_odds(e.probability);
            for (auto v : e.detectors) residual[v] ^= 1;
        }
        for (const auto& path : paths)
        {
            std::vector<size_t> original_path;
            for (auto id : path) original_path.push_back(sub.original_events[id]);
            candidate.second_stage_paths.push_back(std::move(original_path));
        }
        if (std::any_of(residual.begin(), residual.end(), [](auto x) { return x != 0; }))
            throw std::runtime_error("Concatenated matching did not reproduce full syndrome");
        out[c] = std::move(candidate);
    }
    return out;
}
} // namespace dec::concat
