/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 * */
#include "decoder/correlated_restriction.h"
#include "decoder/correlated_restriction/circuit.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ostream>
#include <stdexcept>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec
{

CorrelatedRestrictionDecoder::CorrelatedRestrictionDecoder(const stim::DetectorErrorModel& dem, config_type conf)
    :num_detectors(dem.count_detectors()), num_observables(dem.count_observables()), correlated_(conf.correlated)
{
    if (!std::isfinite(conf.physical_p) || conf.physical_p < 0 || conf.physical_p >= 0.1)
        throw std::invalid_argument("Correlated restriction requires 0 <= p < 0.1");
    // At p=0 the DEM has no error terms; the positive scale is used only by
    // the spatial geometry constructor and the unused flag fallback prior.
    const double scale = conf.physical_p == 0 ? 1e-12 : conf.physical_p;
    circuit_ = std::make_unique<correlated_restriction::CircuitRestrictionDecoder>(dem,
        restriction::make_code_capacity_lattice(false,conf.distance,scale),
        conf.memory_x,scale,conf.use_flags,1,true);
}

CorrelatedRestrictionDecoder::CorrelatedRestrictionDecoder(uint32_t distance, uint32_t rounds, double p, bool correlated)
    :num_detectors(restriction::make_code_capacity_lattice(false,distance,p).num_checks*
                   (rounds == 1 ? 1 : rounds+1)), num_observables(1), correlated_(correlated)
{
    if (!std::isfinite(p) || p < 0 || p >= 0.5)
        throw std::invalid_argument("Capacity/phenomenological correlation requires 0 < p < 0.5");
    model_ = std::make_unique<correlated_restriction::Decoder>(distance,rounds,p,correlated);
}

CorrelatedRestrictionDecoder::~CorrelatedRestrictionDecoder() = default;

Result
CorrelatedRestrictionDecoder::decode(SyndromeRef syndrome, ObsRef)
{
    if (syndrome.num_bits_padded() < num_detectors)
        throw std::invalid_argument("Correlated restriction syndrome is too short");
    const auto start = std::chrono::steady_clock::now();
    std::vector<uint8_t> bits(num_detectors);
    for (size_t v = 0; v < num_detectors; v++) bits[v] = syndrome[v];
    uint64_t frame;
    if (circuit_)
            frame = correlated_ ? circuit_->decode_correlated(bits).observable_mask
                                : circuit_->decode(bits).observable_mask;
    else
    {
        const auto out = model_->decode(bits);
        if (!std::isfinite(out.weight)) throw std::runtime_error("Correlated restriction lifting failed");
        frame = out.observable;
    }
    Result result;
    result.flipped_obs = ObsType(num_observables);
    result.flipped_obs[0] = frame & 1;
    result.matching_data.frame_flips = result.flipped_obs;
    result.decode_latency_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-start).count();
    stats_.shots++;
    stats_.decode_total_ns += result.decode_latency_ns;
    stats_.decode_min_ns = std::min(stats_.decode_min_ns,result.decode_latency_ns);
    stats_.decode_max_ns = std::max(stats_.decode_max_ns,result.decode_latency_ns);
    return result;
}

CorrelatedRestrictionDecoder::Statistics
CorrelatedRestrictionDecoder::statistics() const
{
    auto out = stats_;
    if (!out.shots) out.decode_min_ns = 0;
    return out;
}

void
CorrelatedRestrictionDecoder::print_stats(std::ostream& out) const
{
    const auto stats = statistics();
    out << (correlated_ ? "CorrelatedRestrictionDecoder shots=" : "RestrictionDecoder shots=") << stats.shots
        << " average_decode_us=" << (stats.shots ? double(stats.decode_total_ns)/stats.shots/1000 : 0)
        << " min_decode_us=" << double(stats.decode_min_ns)/1000
        << " max_decode_us=" << double(stats.decode_max_ns)/1000 << '\n';
}

void
CorrelatedRestrictionDecoder::mpi_accumulate()
{
#if defined(ENABLE_MPI)
    uint64_t sums[]{stats_.shots,stats_.decode_total_ns};
    MPI_Allreduce(MPI_IN_PLACE,sums,2,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
    stats_.shots = sums[0];
    stats_.decode_total_ns = sums[1];
    MPI_Allreduce(MPI_IN_PLACE,&stats_.decode_min_ns,1,MPI_UINT64_T,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&stats_.decode_max_ns,1,MPI_UINT64_T,MPI_MAX,MPI_COMM_WORLD);
#endif
}
} // namespace dec
