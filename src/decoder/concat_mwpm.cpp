/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 * */
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"

#include <algorithm>
#include <chrono>
#include <ostream>
#include <stdexcept>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec
{

ConcatMWPM::ConcatMWPM(const stim::DetectorErrorModel& dem)
    :ConcatMWPM(dem,config_type{})
{
}

ConcatMWPM::ConcatMWPM(const stim::DetectorErrorModel& dem, config_type conf)
    :num_detectors(dem.count_detectors()), num_observables(dem.count_observables()),
     decoder_(std::make_unique<concat::CircuitDecoder>(dem,conf.memory_x,
                                                       conf.use_flags,conf.flag_cache_capacity,
                                                       conf.boundary_gap,concat::EnsembleConfig{conf.ensemble_size,conf.ensemble_seed})), debug_(conf.debug)
{
}

ConcatMWPM::~ConcatMWPM() = default;

Result
ConcatMWPM::decode(SyndromeRef syndrome, ObsRef)
{
    if (syndrome.num_bits_padded() < num_detectors)
        throw std::invalid_argument("Concatenated MWPM syndrome is too short");
    const auto start = std::chrono::steady_clock::now();
    std::vector<uint8_t> bits(num_detectors);
    for (size_t v = 0; v < num_detectors; v++) bits[v] = syndrome[v];
    Result result;
    const auto frame = decoder_->decode(bits,debug_ ? &result.debug_info : nullptr);
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

ConcatMWPM::Statistics
ConcatMWPM::statistics() const
{
    auto out = stats_;
    if (!out.shots) out.decode_min_ns = 0;
    return out;
}

void
ConcatMWPM::print_stats(std::ostream& out) const
{
    const auto stats = statistics();
    out << "ConcatMWPM shots=" << stats.shots
        << " average_decode_us=" << (stats.shots ? double(stats.decode_total_ns)/stats.shots/1000 : 0)
        << " min_decode_us=" << double(stats.decode_min_ns)/1000
        << " max_decode_us=" << double(stats.decode_max_ns)/1000 << '\n';
}

void
ConcatMWPM::mpi_accumulate()
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
