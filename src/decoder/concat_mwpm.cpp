/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Concatenated MWPM, Lee--Li--Bartlett, Quantum 9, 1609 (2025).
 * */
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"

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

ConcatMWPM::ConcatMWPM(const stim::DetectorErrorModel& dem)
    :ConcatMWPM(dem,config_type{})
{
}

ConcatMWPM::ConcatMWPM(const stim::DetectorErrorModel& dem, config_type conf)
    :num_detectors(dem.count_detectors()), num_observables(dem.count_observables()),
     decoder_(std::make_unique<concat::CircuitDecoder>(dem,conf.memory_x,
                                                       conf.use_flags,conf.flag_cache_capacity,
                                                       concat::EnsembleConfig{conf.ensemble_size,conf.ensemble_seed},
                                                       concat::FlagRenormalization{conf.flag_renormalization,conf.physical_error_rate},
                                                       conf.belief_propagation,conf.boundary_gap)),
     debug_(conf.debug),belief_propagation_(conf.belief_propagation)
{
}

ConcatMWPM::~ConcatMWPM() = default;

Result
ConcatMWPM::decode(SyndromeRef syndrome, ObsRef truth)
{
    if (syndrome.num_bits_padded() < num_detectors)
        throw std::invalid_argument("Concatenated MWPM syndrome is too short");
    const auto start = std::chrono::steady_clock::now();
    std::vector<uint8_t> bits(num_detectors);
    for (size_t v = 0; v < num_detectors; v++) bits[v] = syndrome[v];

    Result result;
    const auto frame=decoder_->decode(bits,debug_ ? &result.debug_info:nullptr);
    const auto candidate_gap=decoder_->candidate_gap();
    if (belief_propagation_ && candidate_gap)
    {
        result.has_candidate_gap=true;
        result.candidate_gap_db=*candidate_gap*(10/std::log(10.0));
    }
    result.flipped_obs = ObsType(num_observables);
    for (auto id : frame.ids()) result.flipped_obs[id] = 1;
    result.matching_data.frame_flips = result.flipped_obs;
    result.decode_latency_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-start).count();
    if (belief_propagation_)
    {
        if (result.has_candidate_gap)
        {
            stats_.candidate_gap_available++;
            const double gap=result.candidate_gap_db;
            s_unsigned_gap.add(gap);
            bool error = false;
            if (truth.num_bits_padded() >= num_observables)
                for (size_t id = 0; id < num_observables; ++id)
                    error |= truth[id] != result.flipped_obs[id];
            if (error) s_unsigned_gap_errors.add(gap);
            if (truth.num_bits_padded()) s_signed_gap.add(error ? -gap:gap);
            result.matching_data.gap=error ? -gap:gap;
        }
        else stats_.candidate_gap_unavailable++;
        if (debug_) result.debug_info+=result.has_candidate_gap ?
            "Candidate frame gap (dB): "+std::to_string(result.candidate_gap_db)+"\n":
            "Candidate frame gap unavailable: no opposing-frame candidate found.\n";
    }
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
    uint64_t confidence[]{stats_.candidate_gap_available,stats_.candidate_gap_unavailable};
    MPI_Allreduce(MPI_IN_PLACE,confidence,2,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
    stats_.candidate_gap_available=confidence[0];stats_.candidate_gap_unavailable=confidence[1];
    s_signed_gap.mpi_accumulate();s_unsigned_gap.mpi_accumulate();s_unsigned_gap_errors.mpi_accumulate();
#endif
}
} // namespace dec
