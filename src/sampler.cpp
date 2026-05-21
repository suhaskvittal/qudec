/*
 *  author: Suhas Vittal
 *  date:   10 May 2026
 * */

#include "sampler.h"

#include <cassert>
#include <queue>
#include <unordered_set>
#include <vector>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

using dem_type = stim::DetectorErrorModel;

std::vector<size_t> _generate_k_random_numbers_without_replacement(size_t k, size_t max, RNG&);

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PROBLEM
generate_syndromes_with_k_errors(const dem_type& dem, size_t k, size_t count, RNG& rng)
{
    // First determine what errors we will sample for each trial.
    // We store this data in a priority queue so it is O(1) to determine whether
    // an error should be sampled by any trial.

    struct pq_entry { size_t error; size_t trial; };
    struct pq_cmp
    {
        bool
        operator()(const pq_entry& a, const pq_entry& b) const
        {
            // return if `a` is larger than `b` as priority queue default implementation is maxheap
            return (a.error > b.error) || ((a.error == b.error) && (a.trial > b.trial));
        }
    };

    using pq_type = std::priority_queue<pq_entry, std::vector<pq_entry>, pq_cmp>;

    // Initialize the priority queue:

    pq_type sampled_errors;
    for (size_t i = 0; i < count; i++)
        for (size_t e : _generate_k_random_numbers_without_replacement(k, dem.count_errors(), rng))
            sampled_errors.push({e, i});

    // go through the DEM to initialize `dets` and `obs`

    SYNDROME_TABLE dets(count, dem.count_detectors());
    SYNDROME_TABLE obs(count, dem.count_observables());
    size_t error_idx{0};
    dem.iter_flatten_error_instructions(
            [&dets, &obs, &sampled_errors, &error_idx] (const auto& inst)
            {
                while (sampled_errors.size() > 0)
                {
                    auto [e, trial] = sampled_errors.top();
                    if (e != error_idx)
                        break;
                    sampled_errors.pop();
                    // update `dets` and `obs`
                    inst.for_separated_targets(
                            [&dets, &obs, trial] (const auto& grp)
                            {
                                for (const auto& t : grp)
                                {
                                    if (t.is_relative_detector_id())
                                        dets[trial][t.val()] ^= 1;
                                    else if (t.is_observable_id())
                                        obs[trial][t.val()] ^= 1;
                                }
                            });
                }
                error_idx++;
            });
    return PROBLEM{dets, obs};
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
 
POLY
compute_probability_polynomial(const dem_type& dem, size_t max_errors)
{
    POLY out(max_errors, 0.0);
    out[0] = 1.0;
    dem.iter_flatten_error_instructions(
            [max_errors, &out] (const auto& inst)
            {
                assert(inst.type == stim::DemInstructionType::DEM_ERROR);
                double prob = inst.arg_data[0];
                // modify `out`:
                POLY prev(out);
                out[0] *= (1-prob);
                for (size_t i = 1; i < max_errors; i++)
                    out[i] = (1-prob)*prev[i] + prob*prev[i-1];
            });
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Helper functions
 * */

namespace
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::vector<size_t>
_generate_k_random_numbers_without_replacement(size_t k, size_t m, RNG& rng)
{
    std::vector<size_t> out(k);
    std::unordered_set<size_t> out_set;
    out_set.reserve(k);

    for (size_t i = 0; i < k; i++)
    {
        size_t e;
        do { e = rng() % m; } while (out_set.count(e));
        out[i] = e;
        out_set.insert(e);
    }

    return out;
}
 
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // anon
 
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
