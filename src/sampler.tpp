/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#include "decoder/logger.h"

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class D> double
estimate_logical_error_rate(const stim::DetectorErrorModel& dem, D& decoder, EXPERIMENT_CONFIG conf)
{
    // generate probability polynomial: 
    POLY prob_x = compute_probability_polynomial(dem, conf.max_level);

    // now start sampling DEM for errors:
    double ler{0.0};

    int world_rank{0};
    int world_size{1};
#if defined(ENABLE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
    const bool can_print_progress = (conf.print_progress && conf.verbosity == 0 && world_rank == 0);

    decoder::LOGGER logger;
    std::mt19937_64 rng;
    rng.seed(conf.seed + world_rank);
    bool no_errors_found_yet{true};
    for (size_t k = conf.start_level; k < conf.max_level; k++)
    {
        if (prob_x[k] < 0.005*ler)
            continue;

        int local_samples = conf.samples_per_level;
#if defined(ENABLE_MPI)
        local_samples = (local_samples/world_size) + (world_rank==0 ? (local_samples % world_size) : 0);
#endif
        auto [syndromes, obs_array] = generate_syndromes_with_k_errors(dem, k, local_samples, rng);

        if (can_print_progress)
            (std::cout << "[ sampling " << k << " errors ]: ").flush();


        // decode the given syndromes and estimate the LER for `k` errors
        uint64_t error_count{0};
        uint64_t samples{0};

        uint64_t prev_progress_error_count{0};
        int local_max_errors = conf.max_errors_per_level;
#if defined(ENABLE_MPI)
        local_max_errors = std::max(conf.max_errors_per_level / world_size, int64_t{4});
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
            auto result = decoder.decode(syndromes[i], logger);
            bool any_mismatch{false};
            for (size_t j = 0; j < dem.count_observables(); j++)
                any_mismatch |= (result.flipped_obs[j] != obs_array[i][j]);
            if (any_mismatch)
                error_count++;
            samples++;

            // print out debug/error info:
            logger.dump_info(std::cout, conf.verbosity);
            if (any_mismatch && conf.verbosity > 0)
                logger.dump_error(std::cout);
            logger.reset();
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
                        << "\n";
        }

        // skip levels if there are no errors
        no_errors_found_yet &= (error_count == 0);
        if (no_errors_found_yet)
            k += conf.skip_levels_after_no_errors_found;
    }

    return ler;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
