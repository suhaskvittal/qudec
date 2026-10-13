/*
 *  contributor: OpenAI Codex
 * */

#include "decoder/pdz_libra.h"

#include <ostream>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PDZLibra::PDZLibra(const stim::DetectorErrorModel& dem, pp::PDZ::config_type pdz_config,
                  Libra::Config libra_config)
    :initial_(dem, {.enable_gap_estimation=false, .report_matching_data=true}),
    predecoder_(dem, pdz_config), residual_(dem, libra_config)
{
    residual_.prepare_ensemble();
}

Result
PDZLibra::decode(SyndromeRef syndrome, ObsRef truth)
{
    // Preserve the original sample for the sampler and its error callbacks.
    const auto initial = initial_.decode(syndrome, truth);
    SyndromeType remaining(predecoder_.num_detectors);
    remaining = syndrome;
    const auto partial = predecoder_.apply(remaining, initial);

    ObsType residual_truth(predecoder_.num_observables);
    residual_truth = truth;
    residual_truth ^= partial.flipped_obs;
    auto result = residual_.decode(remaining, residual_truth);
    result.flipped_obs ^= partial.flipped_obs;
    // Gap and probability retain Libra's residual-syndrome scores.
    result.matching_data.frame_flips = result.flipped_obs;

    bool failed = false;
    for (size_t i = 0; i < predecoder_.num_observables; i++)
        failed |= result.flipped_obs[i] != truth[i];
    counts_[0]++;
    counts_[1] += failed;
    counts_[2] += syndrome.popcnt();
    counts_[3] += remaining.popcnt();
    return result;
}

PDZLibra::Statistics
PDZLibra::statistics() const
{
    return {counts_[0], counts_[1], counts_[2], counts_[3]};
}

void
PDZLibra::print_stats(std::ostream& out) const
{
    const auto stats = statistics();
    out << "PDZLibra sampled_shots=" << stats.sampled_shots
        << " sampled_failures=" << stats.sampled_failures
        << " input_detection_events=" << stats.input_events
        << " residual_detection_events=" << stats.residual_events << '\n';
    initial_.print_stats(out);
    residual_.print_stats(out);
}

void
PDZLibra::mpi_accumulate()
{
#if defined(ENABLE_MPI)
    MPI_Allreduce(MPI_IN_PLACE, counts_.data(), counts_.size(), MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
    initial_.mpi_accumulate();
    residual_.mpi_accumulate();
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec
