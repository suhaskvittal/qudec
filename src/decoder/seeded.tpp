/* Contributor: OpenAI Codex. Template implementation of SeededMatching. */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <set>
#include <stdexcept>
#include <type_traits>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec
{

template <typename D1, typename D2>
SeededMatching<D1, D2>::SeededMatching(const stim::DetectorErrorModel& dem, D1& seed,
                                      D2& second, SeededMatchingConfig config)
    : num_detectors(dem.count_detectors()), num_observables(dem.count_observables()),
      seed_(seed), second_(second), config_(config)
{
    static_assert(std::is_same_v<D1, PyMatching> || std::is_same_v<D1, UnionFind>);
    static_assert(std::is_same_v<D2, Libra> || std::is_same_v<D2, CorrelatedPyMatching>);
    if (seed.num_detectors != num_detectors || seed.num_observables != num_observables ||
        second.num_detectors != num_detectors || second.num_observables != num_observables)
        throw std::invalid_argument("Seeded decoder dimensions disagree");
    if (!seed.conf.report_matching_data)
        throw std::invalid_argument("Seed decoder must report matching data");
    if constexpr (std::is_same_v<D2, Libra>)
    {
        if (!std::isfinite(config_.libra_gap_threshold_db) || config_.libra_gap_threshold_db < 0)
            throw std::invalid_argument("Libra seed gap threshold must be finite and nonnegative");
        if (config_.libra_gap_mode == SeedGapMode::Cluster && !seed.conf.enable_cluster_gap)
            throw std::invalid_argument("Cluster-gap gating requires cluster gap from the seed");
        if (config_.libra_gap_mode == SeedGapMode::Complementary)
        {
            if constexpr (std::is_same_v<D1, PyMatching>)
                complementary_gap_ = std::make_unique<PyMatching>(dem,
                    PyMatching::config_type{.enable_gap_estimation=true});
            else
                throw std::invalid_argument("Complementary gap requires a PyMatching seed");
        }
        if (config_.postprocess)
        {
            libra_conditioning_ = std::make_unique<CorrelatedPyMatching>(dem);
            postprocessor_ = std::make_unique<SeededPostprocessor>(
                libra_conditioning_->matcher(), num_detectors, num_observables);
        }
    }
    else if (config_.postprocess)
        postprocessor_ = std::make_unique<SeededPostprocessor>(
            second_.matcher(), num_detectors, num_observables);
}

template <typename D1, typename D2>
SeededMatching<D1, D2>::~SeededMatching() = default;

template <typename D1, typename D2>
Result SeededMatching<D1, D2>::decode(SyndromeRef syndrome, ObsRef truth)
{
    const auto started = std::chrono::steady_clock::now();
    ++shots_;
    input_events_ += syndrome.popcnt();
    Result seed_result = seed_.decode(syndrome, truth);
    const auto finish = [&](Result result)
    {
        result.matching_data.frame_flips = result.flipped_obs;
        result.decode_latency_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - started).count();
        pipeline_total_ns_ += result.decode_latency_ns;
        bool failure = false;
        for (size_t obs = 0; obs < num_observables; ++obs)
            failure |= result.flipped_obs[obs] != truth[obs];
        failures_ += failure;
        return result;
    };

    if constexpr (std::is_same_v<D2, Libra>)
    {
        if (syndrome.popcnt() == 0)
        {
            ++gap_bypasses_;
            return finish(std::move(seed_result));
        }
        double gap = 0;
        if (config_.libra_gap_mode == SeedGapMode::Cluster)
        {
            if (!seed_result.has_cluster_gap)
                throw std::logic_error("Seed result did not provide cluster gap");
            gap = seed_result.cluster_gap_db;
        }
        else
            gap = std::abs(complementary_gap_->decode(syndrome, truth).matching_data.gap);
        if (gap >= config_.libra_gap_threshold_db)
        {
            ++gap_bypasses_;
            return finish(std::move(seed_result));
        }
    }

    // Physical edges are a set: the same edge appearing in two pair paths is
    // selected once when conditioning the correlated decoding graph.
    std::vector<int64_t> selected_edges;
    std::set<std::pair<int64_t, int64_t>> seen;
    for (const auto& pair : seed_result.matching_data.assignments)
        for (const auto& edge : pair.correction_edges)
        {
            int64_t u = edge.d1, v = edge.d2;
            if (v >= 0 && u > v) std::swap(u, v);
            if (seen.emplace(u, v).second)
            {
                selected_edges.push_back(u);
                selected_edges.push_back(v);
            }
        }

    auto* conditioning = [&]() -> CorrelatedPyMatching*
    {
        if constexpr (std::is_same_v<D2, CorrelatedPyMatching>) return &second_;
        else return libra_conditioning_.get();
    }();
    uint64_t reweight_ns = 0;
    if (conditioning)
    {
        if constexpr (std::is_same_v<D2, CorrelatedPyMatching>)
        {
            const auto reweight_started = std::chrono::steady_clock::now();
            conditioning->prepare(selected_edges);
            reweight_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - reweight_started).count();
        }
        else
            conditioning->prepare(selected_edges);
    }
    // Keep the conditioned graph alive through the optional pair checks and
    // final matching pass, including on exception paths.
    struct RestoreGuard
    {
        CorrelatedPyMatching* decoder;
        ~RestoreGuard() { if (decoder) decoder->restore(); }
    } guard{conditioning};

    SyndromeType residual(num_detectors);
    residual = syndrome;
    ObsType fixed_frame(num_observables);
    if (postprocessor_)
    {
        auto filtered = postprocessor_->apply(syndrome, seed_result.matching_data);
        residual = std::move(filtered.remaining);
        fixed_frame = std::move(filtered.fixed_frame);
        frozen_pairs_ += filtered.frozen_pairs;
    }
    residual_events_ += residual.popcnt();

    Result final_result;
    if (residual.popcnt() == 0)
        final_result.flipped_obs = std::move(fixed_frame);
    else
    {
        ++second_calls_;
        if constexpr (std::is_same_v<D2, Libra>)
            final_result = second_.decode_with_ensemble(residual, truth);
        else
            final_result = second_.decode_prepared(residual, truth);
        // The correlated second call includes graph reweighting, which happens
        // before optional postprocessing. Libra's separate conditioning graph
        // serves postprocessing and is not part of its second decoder call.
        const uint64_t second_call_ns = final_result.decode_latency_ns + reweight_ns;
        second_total_ns_ += second_call_ns;
        if (config_.latency_histogram)
            second_decode_latency_us_.add(double(second_call_ns) / 1000.0);
        final_result.flipped_obs ^= fixed_frame;
    }
    return finish(std::move(final_result));
}

template <typename D1, typename D2>
void SeededMatching<D1, D2>::print_stats(std::ostream& out) const
{
    out << "SeededMatching shots=" << shots_ << " failures=" << failures_
        << " second_calls=" << second_calls_ << " gap_bypasses=" << gap_bypasses_
        << " frozen_pairs=" << frozen_pairs_ << " input_events=" << input_events_
        << " residual_events=" << residual_events_ << '\n';
    out << "SeededMatching latency (us): average="
        << (shots_ ? double(pipeline_total_ns_) / shots_ / 1000.0 : 0)
        << " second_weighted="
        << (shots_ ? double(second_total_ns_) / shots_ / 1000.0 : 0) << '\n';
    if (config_.latency_histogram)
    {
        out << "Second decoder latency histogram: 1000 logarithmic bins (1 ns to 1 s; boundaries in us)\n";
        out << second_decode_latency_us_.to_string_full() << '\n';
    }
    seed_.print_stats(out);
    second_.print_stats(out);
    if (complementary_gap_)
    {
        const auto stats = complementary_gap_->statistics();
        out << "SeededMatching complementary gap latency (us): calls=" << stats.timed_shots
            << " average=" << (stats.timed_shots
                ? double(stats.decode_total_ns) / stats.timed_shots / 1000.0 : 0) << '\n';
    }
}

template <typename D1, typename D2>
void SeededMatching<D1, D2>::mpi_accumulate()
{
#if defined(ENABLE_MPI)
    uint64_t counters[]{shots_, failures_, second_calls_, gap_bypasses_, frozen_pairs_,
                        input_events_, residual_events_, pipeline_total_ns_, second_total_ns_};
    MPI_Allreduce(MPI_IN_PLACE, counters, 9, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    shots_ = counters[0]; failures_ = counters[1]; second_calls_ = counters[2];
    gap_bypasses_ = counters[3]; frozen_pairs_ = counters[4]; input_events_ = counters[5];
    residual_events_ = counters[6]; pipeline_total_ns_ = counters[7]; second_total_ns_ = counters[8];
#endif
    if (config_.latency_histogram) second_decode_latency_us_.mpi_accumulate();
    seed_.mpi_accumulate();
    second_.mpi_accumulate();
    if (complementary_gap_) complementary_gap_->mpi_accumulate();
}

} // namespace dec
