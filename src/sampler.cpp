/*
 *  contributor: OpenAI Codex (corrected rare event sampling).
 * Independent DEM count probabilities and exact conditional sampling.
 */
#include "sampler.h"
#include "rare_event.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

Problem generate_syndromes_with_k_errors(const stim::DetectorErrorModel& dem,
                                         size_t k, size_t count, Rng& rng) {
    if (k > dem.count_errors())
        throw std::out_of_range("Requested more active channels than the DEM contains");
    ConditionalErrorSampler sampler(dem, k + 1);
    return sampler.sample_syndromes(k, count, rng);
}

Poly compute_probability_polynomial(const stim::DetectorErrorModel& dem, size_t max_errors) {
    if (!max_errors) return {};
    std::vector<long double> coefficients(max_errors, 0);
    coefficients[0] = 1;
    size_t instructions = 0;
    dem.iter_flatten_error_instructions([&](const auto& inst) {
        const long double p = inst.arg_data[0];
        if (!std::isfinite(p) || p < 0 || p > 1)
            throw std::invalid_argument("Invalid DEM error probability");
        const size_t highest = std::min(++instructions, max_errors - 1);
        for (size_t k = highest; k > 0; --k)
            coefficients[k] = (1 - p) * coefficients[k] + p * coefficients[k - 1];
        coefficients[0] *= 1 - p;
    });
    return Poly(coefficients.begin(), coefficients.end());
}
