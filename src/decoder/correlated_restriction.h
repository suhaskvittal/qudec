/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Single-pass correlated restriction decoder, authored by OpenAI Codex.
 * */
#ifndef DECODER_CORRELATED_RESTRICTION_h
#define DECODER_CORRELATED_RESTRICTION_h

#include "decoder/types.h"

#include <iosfwd>
#include <limits>
#include <memory>

namespace dec
{
namespace correlated_restriction { class CircuitRestrictionDecoder; class Decoder; }

// cc666 matching -> correlation update -> matching -> components -> lifting. Each instance
// owns mutable matching workspaces; use one instance per sampler worker.
class CorrelatedRestrictionDecoder
{
public:
    struct config_type
    {
        uint32_t distance{3};
        double physical_p{0.01};
        bool memory_x{false};
        bool use_flags{true};
        bool correlated{true};
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
    CorrelatedRestrictionDecoder(const stim::DetectorErrorModel&, config_type);
    // One round selects capacity; d rounds select phenomenological noise.
    CorrelatedRestrictionDecoder(uint32_t distance, uint32_t rounds, double p, bool correlated = true);
    ~CorrelatedRestrictionDecoder();
    CorrelatedRestrictionDecoder(const CorrelatedRestrictionDecoder&) = delete;
    CorrelatedRestrictionDecoder& operator=(const CorrelatedRestrictionDecoder&) = delete;

    // The true observable argument is unused; only the syndrome is decoded.
    Result decode(SyndromeRef, ObsRef);
    Statistics statistics() const;
    void print_stats(std::ostream&) const;
    void mpi_accumulate();

private:
    std::unique_ptr<correlated_restriction::CircuitRestrictionDecoder> circuit_;
    std::unique_ptr<correlated_restriction::Decoder> model_;
    bool correlated_;
    Statistics stats_;
};
} // namespace dec
#endif
