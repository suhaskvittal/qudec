/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Lee, Li, Bartlett, Quantum 9, 1609 (2025), Secs. 3 and 4.2.
 *  C++ concatenated MWPM decoder; PyMatching performs MWPM.
 * */
#ifndef DECODER_CONCAT_MWPM_MATCHING_h
#define DECODER_CONCAT_MWPM_MATCHING_h

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include "decoder/concat_mwpm/frame.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <iosfwd>
#include <string>
#include <vector>

namespace dec::concat
{

// Contributor: OpenAI Codex. Spatial/temporal provenance retained for diagnostics.
enum class BoundaryHint { Unknown, Spatial, Temporal };

struct Event
{
    double probability;
    std::vector<size_t> detectors;
    Frame observable;
    BoundaryHint boundary_hint{BoundaryHint::Unknown};
    // Author: OpenAI Codex. Raw single-fault contribution to the paper's C_f.
    bool flag_supported{false};
};

struct Model
{
    std::vector<int> colors;
    std::vector<Event> events;
    bool flag_renormalization{false}; // Prefer supported contributions after projection.
    bool boundary_gap{false}; // Add spatial boundary parity vertices to stage one.
};

struct Result
{
    Frame observable;
    double weight;
    int color;
    std::vector<size_t> events;
    std::vector<std::vector<size_t>> restricted_edges;
    // Ordered original-model event IDs on each second-stage matched path.
    std::vector<std::vector<size_t>> second_stage_paths;
    // Contributor: OpenAI Codex. Ensemble provenance and coalesced second-stage input.
    size_t member{0};
    std::vector<size_t> virtual_detectors;
    // Author: OpenAI Codex. Best sampled opposing-frame margin, natural-log units.
    // Missing means no opposing frame was found; this is not infinite confidence.
    std::optional<double> candidate_gap;
    // Author: OpenAI Codex. Parity sector of the two retained boundaries in stage one.
    int boundary_mask{3};

};

// Contributor: OpenAI Codex. Fixed first-stage graph ensemble; member zero is unperturbed.
struct EnsembleConfig
{
    size_t size{1};
    uint64_t seed{20261024};
    double spread_scale{1};
};

struct Statistics
{
    size_t first_calls{0}, second_calls{0}, duplicates{0}, infeasible{0};
    size_t ensemble_updates{0}; // OpenAI Codex: graphs prepared on demand.
};

// Contributor: OpenAI Codex. Names map local checks to original DEM detectors.
struct DebugTrace
{
    std::ostream& out;
    std::vector<std::string> detector_labels;
};

class Decoder
{
public:
    // Section 3 uses uniform restricted-lattice weights for code capacity.
    // Algorithm 1 uses compressed projected probabilities for noisy readout.
    explicit Decoder(Model, bool uniform_restricted_weights = false,
                     EnsembleConfig config = {},
                     bool edge_weight_updates = false);
    Result decode(const std::vector<uint8_t>& syndrome);
    std::array<Result, 3> decode_colors(const std::vector<uint8_t>& syndrome, bool trace = true);
    // Three color candidates per ensemble member, coalesced by second-stage syndrome.
    std::vector<Result> decode_candidates(const std::vector<uint8_t>& syndrome, bool trace = true, DebugTrace* debug = nullptr);
    size_t ambiguous_boundary_events() const;
    // Verify the actual integer matching/search weights against -log(probability).
    size_t audit_weights();
    // Author: OpenAI Codex. Reuse union-topology matching workspaces; zero
    // probability events are disabled. Updates may select a subset of the
    // constructor's event supports and change their frames and probabilities.
    void update(Model);
    // Author: OpenAI Codex. Same ordered, positive event supports; no model copies.
    void update_dense(const std::vector<double>&,const std::vector<Frame>&);
    void update_dense(const std::vector<double>& probabilities,const std::vector<uint64_t>& low_frames) {
        update_dense(probabilities, std::vector<Frame>(low_frames.begin(), low_frames.end()));
    }
    Model model;
    Statistics statistics;

private:
    void refresh_weights();
    void prepare_ensemble(int color,size_t member);
    struct Graph;
    struct SubDecoder;
    std::array<std::shared_ptr<SubDecoder>, 3> subdecoders_;
    bool boundary_gap_;
    EnsembleConfig config_;
    bool uniform_;
    bool edge_weight_updates_; // OpenAI Codex: sparse API for BP/renormalization.
    std::map<std::vector<size_t>,size_t> event_slots_;
    std::vector<size_t> active_events_;
    std::optional<Result> decode_color(const std::vector<uint8_t>&, int color, size_t member, int boundary_mask, bool trace,
        std::map<std::vector<uint64_t>,bool>& seen, DebugTrace* debug);
};

} // namespace dec::concat
#endif
