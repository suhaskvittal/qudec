// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PYMATCHING2_MWPM_H
#define PYMATCHING2_MWPM_H

#include <array>
#include <map>

#include "pymatching/sparse_blossom/flooder/graph_flooder.h"
#include "pymatching/sparse_blossom/matcher/alternating_tree.h"
#include "pymatching/sparse_blossom/search/search_flooder.h"

namespace pm {

class AltTreeNode;

struct MatchingResult {
    pm::obs_int obs_mask;
    total_weight_int weight;
    MatchingResult();

    bool operator==(const MatchingResult& rhs) const;

    bool operator!=(const MatchingResult& rhs) const;

    MatchingResult(obs_int obs_mask, total_weight_int weight);

    MatchingResult& operator+=(const MatchingResult& rhs);
    MatchingResult operator+(const MatchingResult& rhs) const;
};

struct Mwpm {
    GraphFlooder flooder;
    Arena<AltTreeNode> node_arena;
    SearchFlooder search_flooder;

    Mwpm();
    explicit Mwpm(GraphFlooder flooder);
    Mwpm(GraphFlooder flooder, SearchFlooder search_flooder);
    Mwpm(Mwpm&& other) noexcept;
    Mwpm& operator=(Mwpm&& other) noexcept;

    AltTreeNode* make_child(
        AltTreeNode& parent,
        GraphFillRegion* child_inner_region,
        GraphFillRegion* child_outer_region,
        const CompressedEdge& child_inner_to_outer_edge,
        const CompressedEdge& child_compressed_edge);
    void process_event(const pm::MwpmEvent& event);
    void handle_blossom_shattering(const BlossomShatterEventData& event);
    void shatter_descendants_into_matches_and_freeze(AltTreeNode& alt_tree_node);
    void handle_tree_hitting_boundary(const RegionHitBoundaryEventData& event);
    void handle_region_hit_region(const MwpmEvent event);
    void handle_tree_hitting_match(
        GraphFillRegion* unmatched_region,
        GraphFillRegion* matched_region,
        const CompressedEdge& unmatched_to_matched_edge);
    void handle_tree_hitting_boundary_match(
        GraphFillRegion* unmatched_region,
        GraphFillRegion* matched_region,
        const CompressedEdge& unmatched_to_matched_edge);
    void handle_tree_hitting_self(const RegionHitRegionEventData& event, AltTreeNode* common_ancestor);
    void handle_tree_hitting_other_tree(const RegionHitRegionEventData& event);
    GraphFillRegion* pair_and_shatter_subblossoms_and_extract_matches(GraphFillRegion* region, MatchingResult& res);
    MatchingResult shatter_blossom_and_extract_matches(GraphFillRegion* region);

    GraphFillRegion* pair_and_shatter_subblossoms_and_extract_match_edges(
        GraphFillRegion* region, std::vector<CompressedEdge>& match_edges);
    void shatter_blossom_and_extract_match_edges(GraphFillRegion* region, std::vector<CompressedEdge>& match_edges);
    void extract_paths_from_match_edges(
        const std::vector<CompressedEdge>& match_edges, uint8_t* obs_begin_ptr, pm::total_weight_int& weight);

    void verify_invariants() const;

    void create_detection_event(DetectorNode* node);
    void reset();
};
// Contributor: OpenAI Codex. Sparse, scoped absolute undirected weight overrides.
// Handles remain valid only while the owning Mwpm and graph topology stay fixed.
// Units are the graph's even, nonnegative integer weights, not probabilities.
// Separate nested scopes restore the enclosing BP/ensemble weights on clear.
class EdgeWeightModifications {
   public:
    explicit EdgeWeightModifications(Mwpm& mwpm);
    ~EdgeWeightModifications();
    EdgeWeightModifications(const EdgeWeightModifications&) = delete;
    EdgeWeightModifications& operator=(const EdgeWeightModifications&) = delete;

    // v == -1 denotes the boundary. Reversed endpoints return the same ID.
    size_t edge_id(int64_t u, int64_t v);
    weight_int edge_weight(size_t id) const;
    void modify_edge_weight(size_t id, int64_t weight);
    void modify_edge_weight(int64_t u, int64_t v, int64_t weight);
    void clear_edge_weight_modifications() noexcept;
    // Preserve existing correlated matching's lower-only implied-weight rule.
    void apply_correlations(const std::vector<int64_t>& edges);
    size_t modified_edge_count() const { return touched_.size(); }

   private:
    struct Edge {
        std::array<weight_int*, 4> pointers;
        std::array<weight_int, 4> previous{};
        bool modified{false};
    };
    Mwpm& mwpm_;
    std::map<std::pair<int64_t, int64_t>, size_t> ids_;
    std::vector<Edge> edges_;
    std::vector<size_t> touched_;
    bool legacy_negative_correlations_{false};
};

}  // namespace pm

#endif  // PYMATCHING2_MWPM_H
