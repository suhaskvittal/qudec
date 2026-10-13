/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#include "rare_event.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

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
    int world_rank{0};
    int world_size{1};
#if defined(ENABLE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif

    const bool can_print_progress = (conf.print_progress && conf.verbosity == 0 && world_rank == 0);

    constexpr size_t SHOTS_PER_BATCH{1024};

    // initialize sampler:
    std::mt19937_64 rng;
    rng.seed(conf.seed + world_rank);
    stim::DemSampler<64> sampler(dem, std::move(rng), SHOTS_PER_BATCH);

    const size_t num_observables = dem.count_observables();

    int64_t t{0},
            errors{0};
    while (t < conf.monte_carlo.max_samples && errors < conf.monte_carlo.stop_at_error_count)
    {
        sampler.resample(false);
        // detectors and observables are in `sampler.det_buffer` and `sampler.obs_buffer`
        // (each stored as [event x shot]); transpose so each row is one shot.
        auto detector_table = sampler.det_buffer.transposed();
        auto obs_table = sampler.obs_buffer.transposed();

        const int64_t batch_samples = std::min<int64_t>(conf.monte_carlo.max_samples - t,
                                                       world_size * SHOTS_PER_BATCH);
        const size_t local_samples = batch_samples / world_size + (world_rank < batch_samples % world_size);
        uint64_t new_errors{0};
        for (size_t i = 0; i < local_samples; i++)
        {
            auto result = decoder.decode(detector_table[i], obs_table[i]);
            if (is_decoding_error(result, obs_table[i], num_observables))
            {
                new_errors++;
                if (conf.verbosity)
                    error_callback(detector_table[i], obs_table[i], result);
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
        conf.rare_event.tail_tolerance >= 1)
        throw std::invalid_argument("Invalid rare_event sample counts, count limit, or tail tolerance");

    int world_rank = 0, world_size = 1;
#if defined(ENABLE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
    const size_t limit = std::min<size_t>(conf.rare_event.max_level, dem.count_errors() + 1);
    ConditionalErrorSampler conditional(dem, limit);
    const double tail_bound = conditional.tail_probability_bound(limit);
    if (tail_bound > conf.rare_event.tail_tolerance)
        throw std::invalid_argument("rare_event omitted-count bound exceeds --rare-tail-tolerance; increase --rare-max-level");
    const auto& probabilities = conditional.count_probabilities();
    const size_t observables = dem.count_observables();
    Rng rng(static_cast<uint64_t>(conf.seed) + world_rank);
    double ler = 0;
    constexpr uint64_t batch_size = 256;
    const bool progress = conf.print_progress && world_rank == 0;
    for (size_t k = 0; k < limit; ++k) {
        if (!(probabilities[k] > 0)) continue;
        uint64_t samples = 0, errors = 0;
        const uint64_t maximum = conf.rare_event.samples_per_level;
        // Synchronize batch stopping across ranks, using global errors/counts.
        while (samples < maximum && errors < uint64_t(conf.rare_event.max_errors_per_level)) {
            const uint64_t count = std::min<uint64_t>(maximum - samples, batch_size * world_size);
            const size_t local = count / world_size + (uint64_t(world_rank) < count % world_size);
            auto [syndromes, truths] = conditional.sample_syndromes(k, local, rng);
            uint64_t new_errors = 0;
            for (size_t i = 0; i < local; ++i) {
                auto result = decoder.decode(syndromes[i], truths[i]);
                if (is_decoding_error(result, truths[i], observables)) {
                    ++new_errors;
                    if (conf.verbosity) error_callback(syndromes[i], truths[i], result);
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
                  << tail_bound << " (K >= " << limit << ")\n";

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
