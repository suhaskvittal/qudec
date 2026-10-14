/* Contributor: OpenAI Codex. Configurable two-decoder PDZ pipeline. */
#ifndef DECODER_PDZ_h
#define DECODER_PDZ_h

#include "preprocessor/pdz.h"

#include <array>
#include <chrono>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec
{

// Non-owning decoder pointers: both decoders and the DEM must outlive this
// pipeline. Use one instance per worker. Initial must report matched endpoints.
template <typename Initial, typename Residual>
class PDZDecoder
{
public:
    struct Statistics
    {
        uint64_t sampled_shots;
        uint64_t sampled_failures;
        uint64_t input_events;
        uint64_t residual_events;
        uint64_t bypassed_shots;
        uint64_t residual_called_shots;
        uint64_t residual_total_ns;
        uint64_t frame_changes;
        uint64_t pipeline_total_ns;
    };

    struct DetailedResult
    {
        Result initial;
        Result partial;
        Result result;
        SyndromeType remaining{0};
        bool bypassed{false};
        bool ensembled{false};
    };

    PDZDecoder(const stim::DetectorErrorModel& dem, Initial* initial, Residual* residual,
               pp::PDZ::config_type config = {}, double removal_threshold = 0.5,
               std::string initial_name = "initial", std::string residual_name = "residual")
        : initial_(initial), residual_(residual), predecoder_(dem, config),
          removal_threshold_(removal_threshold), initial_name_(std::move(initial_name)),
          residual_name_(std::move(residual_name))
    {
        if (!initial_ || !residual_)
            throw std::invalid_argument("PDZ requires two nonnull decoder pointers");
        if (!std::isfinite(removal_threshold_) || removal_threshold_ < 0 || removal_threshold_ > 1)
            throw std::invalid_argument("PDZ removal threshold must be in [0, 1]");
        if (initial_->num_detectors != dem.count_detectors() ||
            residual_->num_detectors != dem.count_detectors() ||
            initial_->num_observables != dem.count_observables() ||
            residual_->num_observables != dem.count_observables())
            throw std::invalid_argument("PDZ decoder geometry does not match the DEM");
        if constexpr (requires { initial_->conf.report_matching_data; })
            if (!initial_->conf.report_matching_data)
                throw std::invalid_argument("PDZ initial decoder must enable report_matching_data");
    }

    pp::PDZ& predecoder() { return predecoder_; }

    DetailedResult decode_detailed(SyndromeRef syndrome, ObsRef truth)
    {
        const auto start = std::chrono::steady_clock::now();
        DetailedResult out;
        out.initial = initial_->decode(syndrome, truth);
        out.remaining = SyndromeType(predecoder_.num_detectors);
        out.remaining = syndrome;
        out.partial = predecoder_.apply(out.remaining, out.initial);
        const size_t input_hw = syndrome.popcnt();
        const size_t post_hw = out.remaining.popcnt();
        out.bypassed = post_hw == 0 || (input_hw > 0 &&
            double(input_hw - post_hw) >= removal_threshold_ * double(input_hw));
        out.result = out.initial;
        if (!out.bypassed)
        {
            ObsType residual_truth(predecoder_.num_observables);
            residual_truth = truth;
            residual_truth ^= out.partial.flipped_obs;
            uint64_t ensemble_before = 0;
            if constexpr (requires { residual_->statistics().ensemble_shots; })
                ensemble_before = residual_->statistics().ensemble_shots;
            out.result = residual_->decode(out.remaining, residual_truth);
            if constexpr (requires { residual_->statistics().ensemble_shots; })
                out.ensembled = residual_->statistics().ensemble_shots > ensemble_before;
            counts_[5] += out.result.decode_latency_ns;
            out.result.flipped_obs ^= out.partial.flipped_obs;
            out.result.matching_data.frame_flips = out.result.flipped_obs;
        }
        bool failed = false;
        for (size_t k = 0; k < predecoder_.num_observables; ++k)
            failed |= out.result.flipped_obs[k] != truth[k];
        ++counts_[0];
        counts_[1] += failed;
        counts_[2] += input_hw;
        counts_[3] += post_hw;
        counts_[4] += out.bypassed;
        counts_[6] += out.result.flipped_obs != out.initial.flipped_obs;
        out.result.decode_latency_ns = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - start).count());
        counts_[7] += out.result.decode_latency_ns;
        return out;
    }

    Result decode(SyndromeRef syndrome, ObsRef truth)
    {
        return decode_detailed(syndrome, truth).result;
    }

    Statistics statistics() const
    {
        return {counts_[0], counts_[1], counts_[2], counts_[3], counts_[4],
                counts_[0] - counts_[4], counts_[5], counts_[6], counts_[7]};
    }

    void print_stats(std::ostream& out) const
    {
        const uint64_t calls = counts_[0] - counts_[4];
        out << "PDZ initial=" << initial_name_ << " residual=" << residual_name_
            << " sampled_shots=" << counts_[0] << " sampled_failures=" << counts_[1]
            << " input_detection_events=" << counts_[2]
            << " residual_detection_events=" << counts_[3]
            << " bypassed_shots=" << counts_[4] << " residual_called_shots=" << calls
            << " frame_changes=" << counts_[6]
            << " pdz_removal_threshold=" << removal_threshold_ << '\n';
        out << "PDZ residual weighted decode latency (us): called_fraction="
            << (counts_[0] ? double(calls) / counts_[0] : 0)
            << " mean_when_called=" << (calls ? double(counts_[5]) / calls / 1000 : 0)
            << " weighted=" << (counts_[0] ? double(counts_[5]) / counts_[0] / 1000 : 0)
            << " scope=full_history excludes_initial_decoder=true\n";
        out << "PDZ full pipeline mean latency (us): average="
            << (counts_[0] ? double(counts_[7]) / counts_[0] / 1000 : 0) << '\n';
        initial_->print_stats(out);
        residual_->print_stats(out);
    }

    void mpi_accumulate()
    {
#if defined(ENABLE_MPI)
        MPI_Allreduce(MPI_IN_PLACE, counts_.data(), static_cast<int>(counts_.size()),
                      MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
        initial_->mpi_accumulate();
        residual_->mpi_accumulate();
    }

private:
    Initial* initial_;
    Residual* residual_;
    pp::PDZ predecoder_;
    double removal_threshold_;
    std::string initial_name_;
    std::string residual_name_;
    // Shots, failures, input/residual bits, bypasses, residual ns, frame changes, total ns.
    std::array<uint64_t, 8> counts_{};
};

} // namespace dec
#endif
