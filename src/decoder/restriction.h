/*
 *  author: OpenAI Codex
 *  date:   21 October 2026
 *
 *  Replaces the unfinished restriction-decoder skeleton by Suhas Vittal.
 * */
#ifndef DECODER_RESTRICTION_h
#define DECODER_RESTRICTION_h

#include "decoder/types.h"
#include "decoder/restriction/spacetime.h"

#include <memory>

namespace dec
{

// Triangular cc666, one logical observable. Uses vendored PyMatching directly.
// Every candidate follows restricted matching -> connected components -> lifting.
class RestrictionDecoder
{
public:
    struct config_type
    {
        uint32_t distance{3};
        bool memory_x{false};
        double physical_p{0.01};
        // Graph reweighting is implemented; circuit-level lifting is experimental.
        // The validated spacetime baseline is phenomenological noise.
        bool use_flags{false};
        size_t ensemble_size{1};
    };

    const size_t num_detectors;
    const size_t num_observables;

    // Undecomposed DEM, at least two extraction rounds, detector coordinates
    // (x,y,time,color,basis,is_flag), as emitted by qudec's color-code circuits.
    RestrictionDecoder(const stim::DetectorErrorModel&, config_type);

    // Code capacity: one syndrome bit per lattice check, perfect measurements.
    // The validated code-capacity ensemble uses 16 trials with a fixed seed.
    explicit RestrictionDecoder(restriction::Lattice, size_t ensemble_size = 16,
                                uint64_t seed = 0);

    Result decode(SyndromeRef, ObsRef);

private:
    std::unique_ptr<restriction::RestrictionDecoder> spatial_;
    std::unique_ptr<restriction::CircuitRestrictionDecoder> spacetime_;
};

} // namespace dec
#endif
