/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

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

        uint64_t new_errors{0};
        for (size_t i = 0; i < SHOTS_PER_BATCH; i++)
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
        t += world_size*SHOTS_PER_BATCH;

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
    // generate probability polynomial: 
    Poly prob_x = compute_probability_polynomial(dem, conf.rare_event.max_level);

    // now start sampling DEM for errors:
    double ler{0.0};

    int world_rank{0};
    int world_size{1};
#if defined(ENABLE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
    const bool can_print_progress = (conf.print_progress && conf.verbosity == 0 && world_rank == 0);

    // choose start level to be `k` with highest `prob_x[k]`
    for (size_t k = conf.rare_event.start_level; k < conf.rare_event.max_level; k++)
        if (prob_x[k] > prob_x[conf.rare_event.start_level])
            conf.rare_event.start_level = k;

    std::mt19937_64 rng;
    rng.seed(conf.seed + world_rank);
    bool no_errors_found_yet{true};
    for (size_t k = conf.rare_event.start_level; k < conf.rare_event.max_level; k++)
    {
        if (prob_x[k] < ler)
            continue;

        int local_samples = conf.rare_event.samples_per_level;
#if defined(ENABLE_MPI)
        local_samples = (local_samples/world_size) + (world_rank==0 ? (local_samples % world_size) : 0);
#endif
        auto [syndromes, obs_array] = generate_syndromes_with_k_errors(dem, k, local_samples, rng);

        if (world_rank == 0)
        {
            if (can_print_progress)
                (std::cout << "[ sampling " << k << " errors ]: ").flush();
            else
                std::cout << "[ ERROR COUNT = " << k << " ] ==========================\n";
        }

        // decode the given syndromes and estimate the LER for `k` errors
        uint64_t error_count{0};
        uint64_t samples{0};

        uint64_t prev_progress_error_count{0};
        int local_max_errors = conf.rare_event.max_errors_per_level;
#if defined(ENABLE_MPI)
        local_max_errors = std::max(conf.rare_event.max_errors_per_level / world_size, int64_t{4});
#endif
        int progress_print_frequency = local_samples / 20;
        for (int i = 0; i < local_samples && error_count < local_max_errors; i++)
        {
#if !defined(ENABLE_MPI)
            if (can_print_progress && (i % progress_print_frequency == 0))
            {
                int progress_tick_value = error_count - prev_progress_error_count;
                if (progress_tick_value > 0)
                    std::cout << " " << progress_tick_value;
                else
                    std::cout << " .";
                std::cout.flush();
                prev_progress_error_count = error_count;
            }
#endif

            // decode and determine if logical error occurred:
            auto result = decoder.decode(syndromes[i], obs_array[i]);
            if (is_decoding_error(result, obs_array[i], dem.count_observables()))
            {
                error_count++;
                if (conf.verbosity)
                    error_callback(syndromes[i], obs_array[i], result);
            }
            samples++;
        }

#if defined(ENABLE_MPI)
        MPI_Allreduce(MPI_IN_PLACE, &error_count, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &samples,      1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif

        double ler_given_k_errors = static_cast<double>(error_count) / static_cast<double>(samples);
        double contrib = ler_given_k_errors * prob_x[k];

        ler += contrib;

        if (can_print_progress)
        {
            std::cout << "\t" << error_count << " of " << samples << " samples had errors"
                        << ", P(error | " << k << " errors) = " << ler_given_k_errors
                        << ", P(" << k << " errors) = " << prob_x[k]
                        << ", total contribution = " << contrib
                        << ", cumulative LER = " << ler
                        << "\n";
        }

        // skip levels if there are no errors
        no_errors_found_yet &= (error_count == 0);
        if (no_errors_found_yet)
            k += conf.rare_event.skip_levels_after_no_errors_found;
    }

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
