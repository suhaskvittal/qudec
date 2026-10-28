/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Concatenated MWPM, Lee--Li--Bartlett, Quantum 9, 1609 (2025).
 * */
#ifndef DECODER_CONCAT_MWPM_h
#define DECODER_CONCAT_MWPM_h

#include "decoder/types.h"
#include "stats.h"

#include <iosfwd>
#include <limits>
#include <memory>

namespace dec
{
namespace concat { class CircuitDecoder; }

// Color-code decoder with two PyMatching stages per color. Each instance
// owns mutable matching workspaces; use one instance per sampler worker.
class ConcatMWPM
{
public:
    struct config_type
    {
        bool memory_x{false};
        bool use_flags{true};
        size_t flag_cache_capacity{16}; // Compatibility; graphs now update in place without a pattern cache.
        bool boundary_gap{false}; // Include complementary spatial boundary parity sectors.
        // Contributor: OpenAI Codex. Includes one unperturbed graph per color.
        size_t ensemble_size{1};
        uint64_t ensemble_seed{20261024};
        bool debug{false}; // Populate Result::debug_info on each decode.
        // Author: OpenAI Codex. Opt-in Sec. III D probability renormalization.
        bool flag_renormalization{false};
        double physical_error_rate{0};
        // Author: OpenAI Codex. Auxiliary-detector BP posterior preprocessing for both stages.
        bool belief_propagation{false};
    };
    struct Statistics
    {
        uint64_t shots{0};
        uint64_t decode_total_ns{0};
        uint64_t decode_min_ns{std::numeric_limits<uint64_t>::max()};
        uint64_t decode_max_ns{0};
        uint64_t candidate_gap_available{0},candidate_gap_unavailable{0};
    };

    const size_t num_detectors;
    const size_t num_observables;

    // Pass the UNDECOMPOSED DEM, with detector coordinates
    // (x,y,time,color,basis,is_flag). Check colors are 0/1/2; flags are
    // colorless; X-check extraction flags condition Z memory and Z-check
    // extraction flags condition X memory.
    explicit ConcatMWPM(const stim::DetectorErrorModel&);
    ConcatMWPM(const stim::DetectorErrorModel&, config_type);
    ~ConcatMWPM();
    ConcatMWPM(const ConcatMWPM&) = delete;
    ConcatMWPM& operator=(const ConcatMWPM&) = delete;

    // Truth is used only for signed confidence statistics, never decoding/selection.
    Result decode(SyndromeRef, ObsRef);
    Statistics statistics() const;
    // Author: OpenAI Codex. dB; signed positive for correct, negative for incorrect.
    StatsHistogram<double> s_signed_gap{"CONCAT_BP_SIGNED_GAP_DB",-256,256,256};
    StatsHistogram<double> s_unsigned_gap{"CONCAT_BP_UNSIGNED_GAP_DB",0,256,256};
    StatsHistogram<double> s_unsigned_gap_errors{"CONCAT_BP_UNSIGNED_GAP_DB_ERRORS",0,256,256};
    void print_stats(std::ostream&) const;
    void mpi_accumulate();

private:
    std::unique_ptr<concat::CircuitDecoder> decoder_;
    Statistics stats_;
    bool debug_{false};
    bool belief_propagation_{false};
};
} // namespace dec
#endif
