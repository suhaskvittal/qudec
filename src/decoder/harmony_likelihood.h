/* Contributor: OpenAI Codex. Original-model error decomposition for Harmony. */
#ifndef QUDEC_HARMONY_LIKELIHOOD_H
#define QUDEC_HARMONY_LIKELIHOOD_H

#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <pymatching/sparse_blossom/driver/user_graph.h>
#include <stim.h>
#include <algorithm>
#include <compare>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace dec::harmony_detail {

constexpr uint32_t boundary = UINT32_MAX;
struct Edge {
    uint32_t a;
    uint32_t b;
    bool obs;
    auto operator<=>(const Edge&) const = default;
};
struct Hypothesis {
    double weight = 0;  // Sum of unperturbed log odds; common likelihood constant omitted.
    bool obs = false;
    std::vector<uint32_t> channel_ids;
};

class ErrorModel {
    struct DisjointSet {
        std::vector<uint32_t> parent;
        explicit DisjointSet(size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
        uint32_t find(uint32_t v) { return parent[v] == v ? v : parent[v] = find(parent[v]); }
        void unite(uint32_t a, uint32_t b) { parent[find(a)] = find(b); }
    };
    struct Channel {
        double weight;
        bool obs;
        std::vector<uint32_t> edge_ids;
    };
    std::vector<Edge> edges_;
    std::vector<Channel> channels_;
    std::map<Edge, uint32_t> edge_id_;
    std::vector<uint32_t> best_single_;
    std::map<std::pair<uint32_t, uint32_t>, uint32_t> best_pair_;
    std::unique_ptr<pm::Mwpm> decomposition_;

    static std::vector<uint32_t> odd_ids(std::vector<uint32_t> ids) {
        std::sort(ids.begin(), ids.end());
        std::vector<uint32_t> out;
        for (size_t i = 0; i < ids.size();) {
            size_t j = i + 1;
            while (j < ids.size() && ids[j] == ids[i]) ++j;
            if ((j - i) & 1) out.push_back(ids[i]);
            i = j;
        }
        return out;
    }
    uint32_t insert_edge(Edge edge) {
        if (edge.a > edge.b) std::swap(edge.a, edge.b);
        auto [it, inserted] = edge_id_.emplace(edge, edges_.size());
        if (inserted) edges_.push_back(edge);
        return it->second;
    }
public:
    explicit ErrorModel(const stim::DetectorErrorModel& dem) {
        if (dem.count_observables() != 1 || dem.count_detectors() >= boundary)
            throw std::invalid_argument("Harmony likelihood requires one observable and 32-bit detector IDs");
        DisjointSet basis(dem.count_detectors());
        dem.iter_flatten_error_instructions([&](const auto& inst) {
            inst.for_separated_targets([&](const auto& group) {
                std::vector<uint32_t> detectors;
                for (auto target : group)
                    if (target.is_relative_detector_id()) detectors.push_back(target.val());
                detectors = odd_ids(std::move(detectors));
                if (detectors.size() > 2)
                    throw std::invalid_argument("Harmony requires graphlike decomposed components");
                if (detectors.size() == 2) basis.unite(detectors[0], detectors[1]);
            });
        });
        uint32_t logical_basis = boundary;
        dem.iter_flatten_error_instructions([&](const auto& inst) {
            inst.for_separated_targets([&](const auto& group) {
                std::vector<uint32_t> detectors;
                bool obs = false;
                for (auto target : group) {
                    if (target.is_relative_detector_id()) detectors.push_back(target.val());
                    else if (target.is_observable_id()) obs = !obs;
                }
                if (!obs) return;
                detectors = odd_ids(std::move(detectors));
                if (detectors.size() != 1)
                    throw std::invalid_argument("Harmony observable must label a boundary component");
                const uint32_t component = basis.find(detectors[0]);
                if (logical_basis != boundary && logical_basis != component)
                    throw std::invalid_argument("Harmony observable spans multiple CSS components");
                logical_basis = component;
            });
        });
        if (logical_basis == boundary)
            throw std::invalid_argument("Harmony requires an observable-crossing error");

        dem.iter_flatten_error_instructions([&](const auto& inst) {
            const double p = inst.arg_data[0];
            if (p == 0) return;
            if (!(p > 0 && p < 0.5))
                throw std::invalid_argument("Harmony requires original channel probabilities below 1/2");
            std::vector<uint32_t> detectors;
            bool obs = false;
            for (auto target : inst.target_data) {
                if (target.is_relative_detector_id()) detectors.push_back(target.val());
                else if (target.is_observable_id()) obs = !obs;
            }
            detectors = odd_ids(std::move(detectors));
            if (detectors.empty()) {
                if (obs) throw std::invalid_argument("Harmony does not support undetectable logical channels");
                return;
            }
            std::map<uint32_t, std::vector<uint32_t>> components;
            for (uint32_t detector : detectors) components[basis.find(detector)].push_back(detector);
            if (components.size() > 2)
                throw std::invalid_argument("Harmony channel spans more than two CSS components");
            Channel channel{std::log((1 - p) / p), obs, {}};
            for (const auto& [component, ids] : components) {
                if (ids.empty() || ids.size() > 2)
                    throw std::invalid_argument("Harmony channel has more than two detectors in one basis");
                const bool edge_obs = obs && component == logical_basis;
                channel.edge_ids.push_back(insert_edge({ids[0], ids.size() == 2 ? ids[1] : boundary,
                                                        edge_obs}));
            }
            if (obs && !components.count(logical_basis))
                throw std::invalid_argument("Logical channel has no logical-basis component");
            const uint32_t id = channels_.size();
            channels_.push_back(channel);
            if (channel.edge_ids.size() == 1) {
                if (best_single_.size() < edges_.size()) best_single_.resize(edges_.size(), boundary);
                auto& best = best_single_[channel.edge_ids[0]];
                if (best == boundary || channels_[best].weight > channel.weight) best = id;
            } else {
                const auto endpoints = std::minmax(channel.edge_ids[0], channel.edge_ids[1]);
                auto& best = best_pair_[{endpoints.first, endpoints.second}];
                if (best == 0 || channels_[best - 1].weight > channel.weight) best = id + 1;
            }
        });
        best_single_.resize(edges_.size(), boundary);
        pm::UserGraph graph(edges_.size());
        for (uint32_t edge = 0; edge < edges_.size(); ++edge) {
            const auto channel = best_single_[edge];
            if (channel != boundary)
                graph.add_or_merge_boundary_edge(edge, {}, channels_[channel].weight,
                    1 / (1 + std::exp(channels_[channel].weight)), pm::DISALLOW);
        }
        for (const auto& [pair, encoded_channel] : best_pair_) {
            const auto& channel = channels_[encoded_channel - 1];
            graph.add_or_merge_edge(pair.first, pair.second, {}, channel.weight,
                1 / (1 + std::exp(channel.weight)), pm::DISALLOW);
        }
        decomposition_ = std::make_unique<pm::Mwpm>(graph.to_mwpm(pm::NUM_DISTINCT_WEIGHTS, true));
    }
    size_t graph_edge_count() const { return edges_.size(); }
    Hypothesis reconstruct(const std::vector<Edge>& matched_edges) {
        std::vector<uint32_t> selected;
        selected.reserve(matched_edges.size());
        bool expected_obs = false;
        for (auto edge : matched_edges) {
            if (edge.a > edge.b) std::swap(edge.a, edge.b);
            const auto found = edge_id_.find(edge);
            if (found == edge_id_.end())
                throw std::runtime_error("Selected matching edge has no original channel signature");
            selected.push_back(found->second);
            expected_obs ^= edge.obs;
        }
        selected = odd_ids(std::move(selected));
        std::vector<uint64_t> events(selected.begin(), selected.end());
        std::vector<int64_t> chosen;
        pm::decode_detection_events_to_edges(*decomposition_, events, chosen);
        Hypothesis result;
        std::vector<uint32_t> parity;
        for (size_t i = 0; i < chosen.size(); i += 2) {
            const auto a = static_cast<uint32_t>(chosen[i]);
            const auto b = chosen[i + 1] == -1 ? boundary : static_cast<uint32_t>(chosen[i + 1]);
            uint32_t channel;
            if (b == boundary) channel = best_single_.at(a);
            else {
                const auto pair = std::minmax(a, b);
                channel = best_pair_.at({pair.first, pair.second}) - 1;
            }
            if (channel == boundary) throw std::runtime_error("Unmapped decomposition channel");
            result.channel_ids.push_back(channel);
            result.weight += channels_[channel].weight;
            result.obs ^= channels_[channel].obs;
            for (auto edge_id : channels_[channel].edge_ids) parity.push_back(edge_id);
        }
        if (odd_ids(std::move(parity)) != selected || result.obs != expected_obs)
            throw std::runtime_error("Error decomposition changed syndrome or observable");
        return result;
    }
};

} // namespace dec::harmony_detail
#endif
