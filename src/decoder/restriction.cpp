/*
 *  author: OpenAI Codex
 *  date:   21 October 2026
 * */
#include "decoder/restriction.h"

#include <chrono>
#include <stdexcept>
#include <utility>
#include <sstream>

namespace dec
{

RestrictionDecoder::RestrictionDecoder(const stim::DetectorErrorModel& dem, config_type conf)
    :num_detectors(dem.count_detectors()), num_observables(dem.count_observables()), debug_(conf.debug),
     spacetime_(std::make_unique<restriction::CircuitRestrictionDecoder>(dem,
         restriction::make_code_capacity_lattice(false, conf.distance, conf.physical_p),
         conf.memory_x, conf.physical_p, conf.use_flags, conf.ensemble_size))
{
}

RestrictionDecoder::RestrictionDecoder(restriction::Lattice lattice, size_t ensemble_size,
                                       uint64_t seed, bool debug)
    :num_detectors(lattice.num_checks), num_observables(1), debug_(debug),
     spatial_(std::make_unique<restriction::RestrictionDecoder>(std::move(lattice),
                                                               ensemble_size, seed))
{
}

Result
RestrictionDecoder::decode(SyndromeRef syndrome, ObsRef)
{
    if (syndrome.num_bits_padded() < num_detectors)
        throw std::invalid_argument("Restriction decoder syndrome is too short");
    const auto start = std::chrono::steady_clock::now();
    restriction::syndrome_type bits(num_detectors);
    for (size_t i = 0; i < num_detectors; i++) bits[i] = syndrome[i];
    uint64_t mask;
    Result result;
    restriction::Result lifted;
    if (spatial_)
    {
        lifted = spatial_->decode(bits);
        mask = lifted.observable_mask;
    }
    else
    {
        const auto correction = spacetime_->decode(bits);
        mask = correction.observable_mask;
        if (debug_) lifted = correction.lifted;
    }
    if (debug_)
    {
        std::ostringstream out;
        out << "Restriction debug: " << (spatial_ ? "spatial" : "spacetime")
            << " observable=" << mask
            << "\nLifted correction:";
        for (auto q : lifted.correction) out << " q" << q;
        out << "\nAssignments=" << lifted.diagnostics.assignments << " components=" << lifted.diagnostics.components
            << " local_lifts=" << lifted.diagnostics.local_lifts << '\n';
        result.debug_info = out.str();
    }
    result.flipped_obs = ObsType(num_observables);
    result.flipped_obs[0] = mask & 1;
    result.matching_data.frame_flips = result.flipped_obs;
    result.decode_latency_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start).count();
    return result;
}

} // namespace dec
