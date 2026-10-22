/*
 *  author: OpenAI Codex
 *  date:   21 October 2026
 * */
#include "decoder/restriction.h"

#include <chrono>
#include <stdexcept>
#include <utility>

namespace dec
{

RestrictionDecoder::RestrictionDecoder(const stim::DetectorErrorModel& dem, config_type conf)
    :num_detectors(dem.count_detectors()), num_observables(dem.count_observables()),
     spacetime_(std::make_unique<restriction::CircuitRestrictionDecoder>(dem,
         restriction::make_code_capacity_lattice(false, conf.distance, conf.physical_p),
         conf.memory_x, conf.physical_p, conf.use_flags, conf.ensemble_size))
{
}

RestrictionDecoder::RestrictionDecoder(restriction::Lattice lattice, size_t ensemble_size,
                                       uint64_t seed)
    :num_detectors(lattice.num_checks), num_observables(1),
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
    const auto mask = spatial_ ? spatial_->decode(bits).observable_mask
                              : spacetime_->decode(bits).observable_mask;
    Result result;
    result.flipped_obs = ObsType(num_observables);
    result.flipped_obs[0] = mask & 1;
    result.matching_data.frame_flips = result.flipped_obs;
    result.decode_latency_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start).count();
    return result;
}

} // namespace dec
