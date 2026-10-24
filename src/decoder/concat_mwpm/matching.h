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

// Contributor: OpenAI Codex. Provenance of projected singleton endpoints.
// Unknown endpoints retain colorless sink access in every boundary-parity sector.
enum class BoundaryHint { Unknown, Spatial, Temporal };

struct Event
{
    double probability;
    std::vector<size_t> detectors;
    uint64_t observable;
    BoundaryHint boundary_hint{BoundaryHint::Unknown};
    // Author: OpenAI Codex. Raw single-fault contribution to the paper's C_f.
    bool flag_supported{false};
};

struct Model
{
    std::vector<int> colors;
    std::vector<Event> events;
    bool flag_renormalization{false}; // Prefer supported contributions after projection.
};

struct Result
{
    uint64_t observable;
    double weight;
    int color;
    std::vector<size_t> events;
    std::vector<std::vector<size_t>> restricted_edges;
    // Ordered original-model event IDs on each second-stage matched path.
    std::vector<std::vector<size_t>> second_stage_paths;
    // Contributor: OpenAI Codex. Bits fire boundary parity vertices of colors (color+1)%3 and (color+2)%3.
    // In boundary-gap mode zero means even and one means odd. Spatial boundary
    // vertex IDs in restricted_edges are model.colors.size()+color.
    int boundary_mask{3};
    // Contributor: OpenAI Codex. Ensemble provenance and coalesced second-stage input.
    size_t member{0};
    std::vector<size_t> virtual_detectors;
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
                     bool boundary_gap = false, EnsembleConfig config = {});
    Result decode(const std::vector<uint8_t>& syndrome);
    std::array<Result, 3> decode_colors(const std::vector<uint8_t>& syndrome, bool trace = true);
    // Contributor: OpenAI Codex. Three original or every feasible boundary-parity candidate (normally six).
    std::vector<Result> decode_candidates(const std::vector<uint8_t>& syndrome, bool trace = true, DebugTrace* debug = nullptr);
    size_t ambiguous_boundary_events() const;
    // Verify the actual integer matching/search weights against -log(probability).
    size_t audit_weights() const;
    // Author: OpenAI Codex. Reuse union-topology matching workspaces; zero
    // probability events are disabled. Updates may select a subset of the
    // constructor's event supports and change their frames and probabilities.
    void update(Model);
    Model model;
    Statistics statistics;

private:
    struct Graph;
    struct SubDecoder;
    std::array<std::shared_ptr<SubDecoder>, 3> subdecoders_;
    bool boundary_gap_;
    EnsembleConfig config_;
    bool uniform_;
    std::map<std::vector<size_t>,size_t> event_slots_;
    std::vector<size_t> active_events_;
    std::optional<Result> decode_color(const std::vector<uint8_t>&, int color, int mask, size_t member, bool trace,
        std::map<std::vector<uint64_t>,bool>& seen, DebugTrace* debug);
};

} // namespace dec::concat
#endif
