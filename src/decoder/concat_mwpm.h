/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Concatenated MWPM, Lee--Li--Bartlett, Quantum 9, 1609 (2025).
 * */
#ifndef DECODER_CONCAT_MWPM_h
#define DECODER_CONCAT_MWPM_h

#include "decoder/types.h"

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
        size_t flag_cache_capacity{16};
    };
    struct Statistics
    {
        uint64_t shots{0};
        uint64_t decode_total_ns{0};
        uint64_t decode_min_ns{std::numeric_limits<uint64_t>::max()};
        uint64_t decode_max_ns{0};
    };

    const size_t num_detectors;
    const size_t num_observables;

    // Pass the UNDECOMPOSED DEM, with detector coordinates
    // (x,y,time,color,basis,is_flag). Check colors are 0/1/2; flags are
    // colorless and flags from both bases condition the matching graphs.
    explicit ConcatMWPM(const stim::DetectorErrorModel&);
    ConcatMWPM(const stim::DetectorErrorModel&, config_type);
    ~ConcatMWPM();
    ConcatMWPM(const ConcatMWPM&) = delete;
    ConcatMWPM& operator=(const ConcatMWPM&) = delete;

    // The true observable argument is unused; only the syndrome is decoded.
    Result decode(SyndromeRef, ObsRef);
    Statistics statistics() const;
    void print_stats(std::ostream&) const;
    void mpi_accumulate();

private:
    std::unique_ptr<concat::CircuitDecoder> decoder_;
    Statistics stats_;
};
} // namespace dec
#endif
