/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Colored-DEM adapter for concatenated MWPM; Bell flags are used when present.
 * */
#ifndef DECODER_CONCAT_MWPM_CIRCUIT_h
#define DECODER_CONCAT_MWPM_CIRCUIT_h

#include "decoder/concat_mwpm/matching.h"
#include "decoder/concat_mwpm/flag_bp.h"
#include <string>
#include <tuple>

namespace dec::concat
{

// Author: OpenAI Codex. Chamberland et al., NJP 22, 023019, Sec. III D.
struct FlagRenormalization
{
    bool enabled{false};
    double physical_error_rate{0};
};

struct CircuitFault
{
    double probability;
    std::vector<size_t> detectors;
    Frame observable;
    // Contributor: OpenAI Codex. Reserved native-generator noise tags.
    BoundaryHint boundary_hint{BoundaryHint::Unknown};
};

struct FlagWeight
{
    std::vector<size_t> checks; // Local memory-sector detector IDs.
    Frame observable;
    double unflagged_probability;
    double effective_probability;
    size_t relevant_flags;
    size_t fired_flags;
    // Contributor: OpenAI Codex. Incidence is diagnostic only; unflagged-event penalties are disabled.
    std::vector<size_t> incident_observed_flags;
    double outside_flag_penalty{1}; // Paper mode: p^m for unsupported ordinary events.
    size_t selected_variant{SIZE_MAX}; // Index in the equivalence class; SIZE_MAX = disabled.
    size_t flag_distance{SIZE_MAX}; // |selected flags XOR full observed flags|.
    std::vector<size_t> selected_flags;
    double representative_probability{0}; // Before conditional normalization.
    double conditional_denominator{1};
};

struct CircuitResult
{
    Result correction;
    std::array<Result, 3> colors;
    std::vector<FlagWeight> weights;
    Model model;
    size_t audited_graph_edges;
    // Contributor: OpenAI Codex. All feasible ensemble/color candidates.
    std::vector<Result> candidates;
    size_t ambiguous_boundary_events{0};
};

class CircuitDecoder
{
public:
    // Uses the full original DEM and infers each flag's isolated-flip prior.
    CircuitDecoder(const stim::DetectorErrorModel&, bool memory_x,
                   bool use_flags = true, size_t cache_capacity = 16,
                   EnsembleConfig ensemble = {}, FlagRenormalization renormalization = {},
                   bool belief_propagation = false, bool boundary_gap = false);
    Frame decode(const std::vector<uint8_t>& full_syndrome, std::string* debug_info = nullptr);
    // Author: OpenAI Codex. Natural-log gap from the most recent serial decode.
    std::optional<double> candidate_gap() const { return candidate_gap_; }

    CircuitResult inspect(const std::vector<uint8_t>& full_syndrome, DebugTrace* debug = nullptr);
    std::vector<uint8_t> fault_syndrome(const std::vector<size_t>& ids) const;
    Frame fault_observable(const std::vector<size_t>& ids) const;

    const size_t num_detectors;
    const std::vector<CircuitFault> faults;
    const std::map<uint64_t, std::vector<double>> coordinates;
    std::vector<size_t> check_globals, flag_globals;
    size_t class_count() const { return classes_.size(); }
    const std::vector<double>& flag_measurement_probabilities() const { return measurement_p_; }
    // Contributor: OpenAI Codex. Full per-flag parity marginals, not flag-only priors.
    const std::vector<double>& flag_flip_probabilities() const { return flag_p_; }

    // Contributor: OpenAI Codex. Class membership depends only on check support.
    // Logical effects are per-variant metadata, never part of the class key.
    struct Variant
    {
        std::vector<size_t> flags;
        double probability;
        Frame observable;
    };
    struct Class
    {
        std::vector<size_t> checks; // Local memory-sector detector IDs, including time.
        std::vector<size_t> flags;
        std::vector<Variant> variants;
        BoundaryHint boundary_hint{BoundaryHint::Unknown};
        bool explicit_boundary_provenance{false};
    };
    // Inspect detector-only class construction and each member's logical effect.
    const std::vector<Class>& equivalence_classes() const { return classes_; }

private:
    struct Context
    {
        std::unique_ptr<Decoder> decoder;
        std::vector<FlagWeight> weights;
        Model input_model;
    };
    std::vector<double> flag_p_; // Marginal P(flag=1), indexed by original detector ID.
    std::vector<double> measurement_p_; // Indexed by original global flag ID.
    const bool use_flags_;
    const bool boundary_gap_;
    // Contributor: OpenAI Codex. Reused for each independently flag-conditioned context.
    const EnsembleConfig ensemble_;
    const FlagRenormalization renormalization_;
    // Author: OpenAI Codex. Flags and opposite-basis checks condition BP;
    // only memory-basis checks enter either matching stage.
    std::vector<size_t> bp_aux_globals_;
    // Author: OpenAI Codex. Posterior auxiliary inference replaces representative selection.
    std::unique_ptr<FlagBP> flag_bp_;
    std::vector<int> colors_;
    std::vector<Class> classes_;
    // Author: OpenAI Codex. One persistent decoder replaces whole-pattern caching.
    Context context_;
    std::string current_flags_;
    bool context_ready_{false};
    std::optional<double> candidate_gap_;

    std::map<size_t,std::vector<size_t>> flag_neighborhoods(const std::vector<uint8_t>&) const;
    Context make_context(const std::vector<uint8_t>&) const;
    // Author: OpenAI Codex. Fixed class/frame slots and reusable per-shot numeric buffers.
    void update_bp_context(const std::vector<uint8_t>&);
    std::vector<size_t> bp_frame_offsets_,bp_variable_slots_;
    std::vector<Frame> bp_frame_labels_,bp_frames_;
    std::vector<double> bp_frame_mass_,bp_probabilities_;
    Context& context(const std::vector<uint8_t>&);
    std::vector<uint8_t> memory_syndrome(const std::vector<uint8_t>&) const;
};
} // namespace dec::concat
#endif
