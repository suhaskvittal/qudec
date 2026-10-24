/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 * */
#ifndef DECODER_CORRELATED_RESTRICTION_CIRCUIT_h
#define DECODER_CORRELATED_RESTRICTION_CIRCUIT_h

#include "decoder/correlated_restriction/matching.h"

#include <string>
#include <unordered_map>

namespace dec::correlated_restriction
{
using namespace dec::restriction;


struct Variant
{
    uint64_t frame{0};
    double score{0};
    size_t iterations{0}, period{0};
    std::string termination;
    bool chose_initial{false};
};


struct ErrorEvent
{
    double probability;
    std::vector<det_id_type> checks;
    std::vector<size_t> flags;
    uint64_t observable_mask;
    // Spatial data representative with the same folded syndrome and frame.
    std::vector<size_t> spatial_faces;
};

struct CircuitResult
{
    uint64_t observable_mask;
    restriction::Result lifted;
    size_t enabled_flag_errors;
    size_t ensemble_member{0};
    std::array<uint64_t,3> reference_frames{};
    std::array<size_t,3> reference_weights{};
};

struct RepresentativeWeight
{
    size_t event;
    double unflagged_probability;
    double flagged_probability;
    double effective_probability;
    size_t active_flags;
};

struct GraphWeight
{
    int ignored_color;
    edge_type endpoints;
    double probability;
    double log_odds_weight;
    double matching_weight;
    double search_weight;
    double rounding_tolerance;
};

struct ReweightingDiagnostics
{
    // Mean inferred prior in corrected mode; scalar 5p in legacy mode.
    double measurement_probability;
    std::vector<RepresentativeWeight> representatives;
    std::vector<GraphWeight> graph_edges;
};

struct MatchingDiagnostics
{
    // Original DEM detector IDs; -1 is PyMatching's collapsed boundary.
    std::array<std::vector<edge_type>, 3> projected_edges;
    std::array<double, 3> projected_weights{};
    std::vector<size_t> canonical_global_detectors;
};

class CircuitRestrictionDecoder
{
    friend struct CircuitTestAccess;
public:
    // Only the chosen memory CSS sector is matched; flags from BOTH sectors
    // condition its error classes. Other-sector check correlations are ignored.
    CircuitRestrictionDecoder(const stim::DetectorErrorModel&, Lattice,
                              bool memory_x, double physical_p, bool use_flags = true,
                              size_t ensemble_size = 1, bool corrected_flag_weights = true);
    CircuitResult decode(const syndrome_type& full_syndrome, std::string* debug_info = nullptr);
    // Original and single-update candidates, for paired diagnostics.
    std::array<Variant,2> decode_single_pass(const syndrome_type&);
    CircuitResult decode_correlated(const syndrome_type&, std::string* debug_info = nullptr);
    // Every candidate still follows matching -> components -> lifting.
    // Allows a model-specific, time-resolved likelihood score after lifting.
    std::vector<CircuitResult> decode_candidates(const syndrome_type&);
    const std::vector<ErrorEvent>& error_events() const { return errors_; }
    syndrome_type event_syndrome(size_t event) const;
    size_t num_detectors() const { return num_detectors_; }
    ReweightingDiagnostics inspect_reweighting(const syndrome_type&) const;
    MatchingDiagnostics inspect_matching(const syndrome_type&);
    MatchingDiagnostics inspect_one_pass_matching(const syndrome_type&);

private:
    struct GraphEdge
    {
        double probability{0};
        double representative_probability{0};
        size_t event{SIZE_MAX};
        edge_type endpoints;
        std::vector<det_id_type> spatial_path;
    };
    struct ProjectionData
    {
        edge_type key{-2,-2};
        edge_type endpoints;
        std::vector<det_id_type> spatial_path;
    };
    std::vector<std::array<ProjectionData,3>> projections_;
    struct Context
    {
        std::array<pm::Mwpm, 3> matchers;
        std::array<std::map<edge_type, GraphEdge>, 3> edges;
        std::map<size_t, double> representatives;
        size_t enabled_flag_errors{0};
        std::vector<RepresentativeWeight> weights;
        std::vector<std::unique_ptr<Context>> alternatives;
    };
    struct MatchingOutput
    {
        std::vector<Assignment> assignments;
        std::array<std::vector<edge_type>, 3> graph_edges;
    };

    const size_t num_detectors_;
    const double physical_p_;
    const bool use_flags_;
    const size_t ensemble_size_;
    const bool corrected_flag_weights_;
    std::vector<double> flag_measurement_p_;
    restriction::RestrictionDecoder lifter_;
    std::vector<Vertex> vertices_;
    std::array<det_id_type, 3> boundaries_;
    std::vector<det_id_type> global_to_check_, check_to_global_;
    std::vector<size_t> flag_globals_;
    std::vector<ErrorEvent> errors_;
    std::map<std::vector<det_id_type>, std::vector<size_t>> classes_;
    std::vector<det_id_type> canonical_to_check_;
    std::unordered_map<std::string, std::unique_ptr<Context>> contexts_;

    std::unique_ptr<Context> make_context(const syndrome_type&, size_t trial = 0) const;
    CircuitResult decode_context(Context&, const syndrome_type&, std::vector<CircuitResult>* candidates = nullptr, std::string* debug_info = nullptr);
    std::unique_ptr<Context> update_context(const Context&, const State&) const;
    CircuitResult lift_matching(const Context&, const syndrome_type&, const MatchingOutput&,
                               std::vector<CircuitResult>* = nullptr);
    void append_debug(const Context&, const MatchingOutput&, const syndrome_type&,
                      const char* stage, std::string&) const;
    Context& context_for(const syndrome_type&);
    MatchingOutput matching(Context&, const syndrome_type&);
    MatchingOutput matching_positive(Context&, const syndrome_type&);
    std::vector<det_id_type> event_spatial_path(const ErrorEvent&, int ignored_color,
                                             det_id_type x, det_id_type y) const;
};

} // namespace dec::correlated_restriction
#endif
