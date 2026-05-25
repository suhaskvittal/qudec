/*
 *  author: Suhas Vittal
 *  date:   10 May 2026
 * */

#ifndef SAMPLER_h
#define SAMPLER_h

#include <stim.h>

#include <cstddef>
#include <random>
#include <utility>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

struct EXPERIMENT_CONFIG
{
    int64_t verbosity{0};
    int64_t samples_per_level{10'000};
    int64_t max_errors_per_level{25};
    int64_t start_level{1};
    int64_t max_level{128};

    int64_t seed{0};

    bool print_progress{false};

    int64_t skip_levels_after_no_errors_found{0};
};

/*
 * Estimates the logical error rate of the given decoder on the
 * given error model.
 * */

template <class D_TYPE>
double estimate_logical_error_rate(const stim::DetectorErrorModel&, D_TYPE&, EXPERIMENT_CONFIG);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

using POLY = std::vector<double>;
using SYNDROME_TABLE = stim::simd_bit_table<64>;
using PROBLEM = std::pair<SYNDROME_TABLE, SYNDROME_TABLE>;
using RNG = std::mt19937_64;

/*
 * This function returns an array containing the probability of `K`
 * errors occuring -- this probability is located at index `K` of the
 * output array.
 *
 * This output array is computed using a generating function that is
 * the product of terms of the form:
 *      `(1-p) + p*x`
 * for each error in the DEM (`p` is the probability of the error). We
 * truncate the polynomial to a size of `max_errors` since the 
 * probablility of `K` errors is `O(p^K)`.
 *
 * This array is need for accurately estimating the logical error rate.
 * */
POLY compute_probability_polynomial(const stim::DetectorErrorModel& dem, size_t max_errors);

/*
 * Returns `count` randomly generated syndromes, each with `k` errors, from the given
 * DEM.
 * */
PROBLEM generate_syndromes_with_k_errors(const stim::DetectorErrorModel&, size_t k, size_t count, RNG&);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#include "sampler.tpp"

#endif // SAMPLER_h
