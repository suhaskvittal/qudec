/*
 *  contributor: OpenAI Codex
 *  purpose: Compose uncorrelated PyMatching, PDZ and Libra.
 * */

#ifndef DECODER_PDZ_LIBRA_h
#define DECODER_PDZ_LIBRA_h

#include "decoder/libra.h"
#include "decoder/surface_code.h"
#include "preprocessor/pdz.h"

#include <array>
#include <cstdint>
#include <iosfwd>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// One instance per worker. The DEM must outlive the instance, as required by PDZ.
class PDZLibra
{
public:
    struct Statistics
    {
        uint64_t sampled_shots;
        uint64_t sampled_failures;
        uint64_t input_events;
        uint64_t residual_events;
    };

    PDZLibra(const stim::DetectorErrorModel&, pp::PDZ::config_type, Libra::Config);

    Result decode(SyndromeRef, ObsRef);
    Statistics statistics() const;
    void print_stats(std::ostream&) const;
    void mpi_accumulate();
private:
    PyMatching initial_;
    pp::PDZ predecoder_;
    Libra residual_;
    std::array<uint64_t, 4> counts_{};
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec

#endif
