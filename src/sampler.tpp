/*
 *  contributor: OpenAI Codex (corrected rare event sampling).
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#include "rare_event.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <memory>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

// Author: OpenAI Codex. A decoder group shares one shot stream. Helpers enter
// every collective decode, but only its first rank samples and counts outcomes.
namespace sampler_detail
{
struct Group
{
    int world_rank{0},world_size{1};
    size_t cores{1},index{0},count{1};
    bool root{true};
};
template<class D> Group group(const D& decoder)
{
    Group g;
#if defined(ENABLE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD,&g.world_rank);MPI_Comm_size(MPI_COMM_WORLD,&g.world_size);
#endif
    g.cores=dec::decoding_core_count(decoder);
    if (!g.cores || g.cores>size_t(g.world_size) || size_t(g.world_size)%g.cores)
        throw std::invalid_argument("MPI world must contain complete decoding_core_count() groups");
    g.index=g.world_rank/g.cores;g.count=g.world_size/g.cores;g.root=g.world_rank%g.cores==0;
    return g;
}
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class D, class ErrorCallback> double
estimate_logical_error_rate(const stim::DetectorErrorModel& dem, 
                            D& decoder, 
                            ExperimentConfig conf, 
                            const ErrorCallback& error_callback)
{
    if (conf.method == "monte_carlo")
        return monte_carlo_sampler(dem, decoder, conf, error_callback);
    else if (conf.method == "rare_event")
        return rare_event_sampler(dem, decoder, conf, error_callback);

    std::cerr << "unknown sampler method: " << conf.method << _die{};
    return 0.0;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class D, class ErrorCallback> double
monte_carlo_sampler(const stim::DetectorErrorModel& dem, 
                    D& decoder, 
                    ExperimentConfig conf, 
                    const ErrorCallback& error_callback)
{
    const auto group=sampler_detail::group(decoder);
    const int world_rank=group.world_rank;

    const bool can_print_progress = (conf.print_progress && conf.verbosity == 0 && world_rank == 0);

    constexpr size_t SHOTS_PER_BATCH{1024};

    // initialize sampler:
    std::mt19937_64 rng;
    rng.seed(conf.seed + group.index);
    std::unique_ptr<stim::DemSampler<64>> sampler;
    if (group.root) sampler=std::make_unique<stim::DemSampler<64>>(dem,std::move(rng),SHOTS_PER_BATCH);
    dec::SyndromeType helper_syndrome(dem.count_detectors());
    dec::ObsType helper_truth(dem.count_observables());

    const size_t num_observables = dem.count_observables();

    int64_t t{0},
            errors{0};
    while (t < conf.monte_carlo.max_samples && errors < conf.monte_carlo.stop_at_error_count)
    {
        if (group.root) sampler->resample(false);
        // detectors and observables are in `sampler.det_buffer` and `sampler.obs_buffer`
        // (each stored as [event x shot]); transpose so each row is one shot.
        SyndromeTable detector_table(0,0),obs_table(0,0);
        if (group.root)
        {
            detector_table=sampler->det_buffer.transposed();
            obs_table=sampler->obs_buffer.transposed();
        }

        const int64_t batch_samples = std::min<int64_t>(conf.monte_carlo.max_samples - t,
                                                       group.count * SHOTS_PER_BATCH);
        const size_t local_samples = batch_samples / group.count + (group.index < batch_samples % group.count);
        uint64_t new_errors{0};
        for (size_t i = 0; i < local_samples; i++)
        {
            auto syndrome=group.root ? detector_table[i]:dec::SyndromeRef(helper_syndrome);
            auto truth=group.root ? obs_table[i]:dec::ObsRef(helper_truth);
            auto result = decoder.decode(syndrome,truth);
            if (!group.root) continue;
            // Contributor: OpenAI Codex. Opt-in diagnostic text is decoder-owned.
            if constexpr (requires { result.debug_info; })
                if (!result.debug_info.empty()) std::cout << result.debug_info;
            if (is_decoding_error(result, truth, num_observables))
            {
                new_errors++;
                if (conf.verbosity)
                    error_callback(syndrome,truth,result);
            }
        }

#if defined(ENABLE_MPI)
        MPI_Allreduce(MPI_IN_PLACE, &new_errors, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
        errors += new_errors;
        t += batch_samples;

        if (can_print_progress)
        {
            double progress_percent = 100.0 * static_cast<double>(t) / static_cast<double>(conf.monte_carlo.max_samples);
            double logical_error_rate = static_cast<double>(errors) / static_cast<double>(t);
            std::cout << "[" << t << " / " << conf.monte_carlo.max_samples << "] (" << progress_percent << "%)"
                        << " -- errors = " << errors
                        << " , logical error rate = " << logical_error_rate
                        << "\n";
        }
    }

    double logical_error_rate = static_cast<double>(errors) / static_cast<double>(t);
    return logical_error_rate;
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class D, class ErrorCallback> double
rare_event_sampler(const stim::DetectorErrorModel& dem, 
                    D& decoder, 
                    ExperimentConfig conf, 
                    const ErrorCallback& error_callback)
{
    if (conf.rare_event.start_level != 0 || conf.rare_event.skip_levels_after_no_errors_found != 0)
        throw std::invalid_argument("rare_event must include levels starting at K=0; skipping levels biases the estimate");
    if (conf.rare_event.max_level <= 0 || conf.rare_event.samples_per_level <= 0 ||
        conf.rare_event.max_errors_per_level <= 0 ||
        !std::isfinite(conf.rare_event.tail_tolerance) || conf.rare_event.tail_tolerance < 0 ||
        conf.rare_event.tail_tolerance >= 1 ||
        !std::isfinite(conf.rare_event.relative_tail_tolerance) ||
        conf.rare_event.relative_tail_tolerance < 0 || conf.rare_event.relative_tail_tolerance > 1)
        throw std::invalid_argument("Invalid rare_event sample counts, count limit, or tail tolerance");

    const auto group=sampler_detail::group(decoder);
    const int world_rank=group.world_rank;
    const size_t limit = std::min<size_t>(conf.rare_event.max_level, dem.count_errors() + 1);
    ConditionalErrorSampler conditional(dem, limit);
    const double tail_bound = conditional.tail_probability_bound(limit);
    if (tail_bound > conf.rare_event.tail_tolerance)
        throw std::invalid_argument("rare_event omitted-count bound exceeds --rare-tail-tolerance; increase --rare-max-level");
    const auto& probabilities = conditional.count_probabilities();
    // Sum backwards instead of subtracting from 1: small tail masses would
    // otherwise lose precision. Include the bound beyond the configured limit.
    std::vector<double> remaining_bounds(limit + 1, tail_bound);
    long double remaining = tail_bound;
    for (size_t k = limit; k-- > 0;)
    {
        remaining += probabilities[k];
        remaining_bounds[k] = std::nextafter(double(std::min(remaining, 1.0L)), 1.0);
    }
    size_t first_omitted = limit;
    double omitted_bound = tail_bound;
    const size_t observables = dem.count_observables();
    Rng rng(static_cast<uint64_t>(conf.seed) + group.index);
    dec::SyndromeType helper_syndrome(dem.count_detectors());
    dec::ObsType helper_truth(dem.count_observables());
    double ler = 0;
    constexpr uint64_t batch_size = 256;
    const bool progress = conf.print_progress && world_rank == 0;
    for (size_t k = 0; k < limit; ++k) {
        if (conf.rare_event.relative_tail_tolerance > 0 && ler > 0 &&
                remaining_bounds[k] < conf.rare_event.relative_tail_tolerance * ler)
        {
            first_omitted = k;
            omitted_bound = remaining_bounds[k];
            if (world_rank == 0)
                std::cout << "Rare-event early stop: first_omitted_level=" << k
                          << " accumulated_ler=" << ler
                          << " remaining_contribution_bound=" << omitted_bound
                          << " relative_tail_tolerance=" << conf.rare_event.relative_tail_tolerance << '\n';
            break;
        }
        if (!(probabilities[k] > 0)) continue;
        uint64_t samples = 0, errors = 0;
        const uint64_t maximum = conf.rare_event.samples_per_level;
        // Synchronize batch stopping across ranks, using global errors/counts.
        while (samples < maximum && errors < uint64_t(conf.rare_event.max_errors_per_level)) {
            const uint64_t count = std::min<uint64_t>(maximum - samples, batch_size * group.count);
            const size_t local = count / group.count + (group.index < count % group.count);
            auto [syndromes, truths] = conditional.sample_syndromes(k,group.root ? local:0,rng);
            uint64_t new_errors = 0;
            for (size_t i = 0; i < local; ++i) {
                auto syndrome=group.root ? syndromes[i]:dec::SyndromeRef(helper_syndrome);
                auto truth=group.root ? truths[i]:dec::ObsRef(helper_truth);
                auto result = decoder.decode(syndrome,truth);
                if (!group.root) continue;
                if constexpr (requires { result.debug_info; })
                    if (!result.debug_info.empty()) std::cout << result.debug_info;
                if (is_decoding_error(result, truth, observables)) {
                    ++new_errors;
                    if (conf.verbosity) error_callback(syndrome,truth,result);
                }
            }
#if defined(ENABLE_MPI)
            MPI_Allreduce(MPI_IN_PLACE, &new_errors, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
            samples += count;
            errors += new_errors;
        }
        const double conditional_rate = double(errors) / samples;
        ler += probabilities[k] * conditional_rate;
        if (progress)
            std::cout << "[K=" << k << "] samples=" << samples << " errors=" << errors
                      << " probability=" << probabilities[k] << " contribution="
                      << probabilities[k] * conditional_rate << '\n';
    }
    if (world_rank == 0)
        std::cout << std::setprecision(12) << "Rare-event omitted-count contribution bound: "
                  << omitted_bound << " (K >= " << first_omitted << ")\n";

    return ler;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class ResultType, class ObsRefType> bool
is_decoding_error(const ResultType& r, ObsRefType obs, size_t obs_count)
{
    bool any_mismatch{false};
    for (size_t i = 0; i < obs_count; i++)
        any_mismatch |= (r.flipped_obs[i] != obs[i]);
    return any_mismatch;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
