/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#include "decoder/surface_code.h"

#include <pymatching/sparse_blossom/driver/user_graph.h>

#include <PerfectMatching.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PYMATCHING::PYMATCHING(const stim::DetectorErrorModel& dem)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    mwpm_(pm::detector_error_model_to_mwpm(dem, pm::NUM_DISTINCT_WEIGHTS))
{}

result_type
PYMATCHING::decode(syndrome_ref syn)
{
    // Collect indices of fired detectors.
    std::vector<uint64_t> det_events;
    for (size_t i = 0; i < num_detectors; i++)
        if (syn[i])
            det_events.push_back(i);
    // Decode.
    result_type res;
    res.flipped_obs = obs_type(num_observables);
    pm::total_weight_int weight = 0;
    pm::decode_detection_events(mwpm_, det_events, res.flipped_obs.u8, weight, false);
    return res;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

using det_id_type = BLOSSOMV::det_id_type;
using adj_entry_type = BLOSSOMV::adj_entry_type;
using adj_list_type = BLOSSOMV::adj_list_type;
using mwpm_edge_type = BLOSSOMV::mwpm_edge_type;
using matching_problem_type = BLOSSOMV::matching_problem_type;

constexpr det_id_type BOUNDARY_ID{-1};

/*
 * Scale used when quantizing edge weights to integers. Blossom-5's `REAL` type is
 * `int` by default, so weights cannot be true 64-bit; `1e6` gives ample resolution
 * for ground-truth accuracy while keeping matched paths well within `INT_MAX/2`.
 * */
constexpr double WEIGHT_SCALE = 1e6;

/*
 * Quantizes a weight to a fixed high-resolution integer.
 * */
uint64_t
_quantize(double w)
{
    return static_cast<uint64_t>(std::round(WEIGHT_SCALE * w));
}

/*
 * Safe adjacency list update that also handles the case where the given detector is
 * already present in the adjacency list (merges the parallel-edge probabilities).
 * */
void
_update_adjacency_list(adj_list_type& adj, det_id_type d, double p, obs_ref frame_flips)
{
    auto adj_it = std::find_if(adj.begin(), adj.end(),
                        [d] (const auto& e) { return e.d == d; });
    if (adj_it != adj.end())
        adj_it->pr = (1-adj_it->pr)*p + (1-p)*adj_it->pr;
    else
        adj.push_back(adj_entry_type{d, p, obs_type{frame_flips}});
}

/*
 * Data structures used when computing pairwise distances via Dijkstra. The weight `w`
 * accumulates the raw (unquantized) `-log(pr)` of each edge along the path; it is only
 * quantized once, after the shortest path is finalized.
 * */
struct distance_type
{
    double   w = std::numeric_limits<double>::max();
    obs_type frame_flips;
};

struct distance_queue_entry
{
    det_id_type d;
    double      w;
};

struct distance_cmp
{
    bool operator()(const distance_queue_entry& a, const distance_queue_entry& b) const { return a.w > b.w; }
};

using distance_queue_type = std::priority_queue<distance_queue_entry,
                                                std::vector<distance_queue_entry>,
                                                distance_cmp>;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

BLOSSOMV::BLOSSOMV(const stim::DetectorErrorModel& dem)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    adj_matrix_(num_detectors)
{
    // setup `adj_matrix_`
    dem.iter_flatten_error_instructions(
            [this] (const auto& inst)
            {
                const double pr = inst.arg_data[0];
                inst.for_separated_targets(
                        [this, pr] (const auto& grp)
                        {
                            std::vector<det_id_type> dets;
                            obs_type frame_flips(num_observables);
                            for (const auto& t : grp)
                            {
                                if (t.is_relative_detector_id())
                                    dets.push_back(static_cast<det_id_type>(t.val()));
                                else if (t.is_observable_id())
                                    frame_flips[t.val()] ^= 1;
                            }

                            det_id_type d1 = dets[0],
                                        d2 = (dets.size() == 1) ? BOUNDARY_ID : dets[1];
                            _update_adjacency_list(adj_matrix_[d1], d2, pr, frame_flips);
                            if (d2 == BOUNDARY_ID)
                                _update_adjacency_list(boundary_adjacency_, d1, pr, frame_flips);
                            else
                                _update_adjacency_list(adj_matrix_[d2], d1, pr, frame_flips);
                        });
            });
}

const adj_list_type&
BLOSSOMV::adj_matrix(det_id_type d) const
{
    return (d == BOUNDARY_ID) ? boundary_adjacency_ : adj_matrix_[d];
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type
BLOSSOMV::decode(syndrome_ref syndrome)
{
    auto detectors = collect_detection_events(syndrome);
    if (detectors.empty())
        return result_type{.flipped_obs=obs_type(num_observables)};

    auto mp = synthesize_matching_problem(std::move(detectors));
    return solve_matching_problem(mp);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::vector<det_id_type>
BLOSSOMV::collect_detection_events(syndrome_ref syndrome) const
{
    std::vector<det_id_type> detectors;
    for (det_id_type i = 0; i < static_cast<det_id_type>(num_detectors); i++)
        if (syndrome[i])
            detectors.push_back(i);

    // Append the boundary only when the count is odd so that a perfect matching
    // exists. The boundary is still a transit node during distance computation
    // regardless, which lets pairs of events route a shortest path through it.
    if (detectors.size() % 2 == 1)
        detectors.push_back(BOUNDARY_ID);

    return detectors;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

matching_problem_type
BLOSSOMV::synthesize_matching_problem(std::vector<det_id_type> detectors) const
{
    const size_t n = detectors.size();

    // map a detector id to its slot in `dist` (boundary lives at index `num_detectors`)
    auto slot = [this] (det_id_type d) -> size_t
                {
                    return (d == BOUNDARY_ID) ? num_detectors : static_cast<size_t>(d);
                };

    // Compute pairwise distances via Dijkstra over the full decoding graph: run from
    // each detector and read off the upper triangle (distances are symmetric).
    const distance_type fill_val{.frame_flips=obs_type(num_observables)};
    std::vector<distance_type> dist(num_detectors+1, fill_val);

    std::vector<mwpm_edge_type> edges;
    edges.reserve(n*(n-1)/2);

    for (size_t ii = 0; ii+1 < n; ii++)
    {
        const det_id_type d1 = detectors[ii];

        std::fill(dist.begin(), dist.end(), fill_val);
        dist[slot(d1)].w = 0.0;
        distance_queue_type pq;
        pq.push({d1, 0.0});
        while (pq.size() > 0)
        {
            auto e = std::move(pq.top());
            pq.pop();
            const det_id_type z1 = e.d;
            const double w1 = e.w;
            const size_t i = slot(z1);
            if (w1 != dist[i].w)
                continue;
            // relax every neighbor (full graph, including the boundary as transit).
            // Accumulate the raw `-log(pr)`; quantization happens once below.
            for (const auto& x : adj_matrix(z1))
            {
                const size_t j = slot(x.d);
                const double w2 = w1 + (-std::log(x.pr));
                if (w2 < dist[j].w)
                {
                    dist[j].w = w2;
                    dist[j].frame_flips = dist[i].frame_flips ^ x.frame_flips;
                    pq.push({x.d, w2});
                }
            }
        }

        // create mwpm edges to all detectors later in `detectors`
        for (size_t jj = ii+1; jj < n; jj++)
        {
            const det_id_type d2 = detectors[jj];
            const size_t j = slot(d2);
            // `dist[j].w` is the accumulated `-log` path weight: the matching's error
            // probability is `exp(-w)`, and the quantized weight is `_quantize(w)`.
            edges.push_back(mwpm_edge_type{ .d1=d1,
                                            .d2=d2,
                                            .pr=std::exp(-dist[j].w),
                                            .w_qu=_quantize(dist[j].w),
                                            .frame_flips=std::move(dist[j].frame_flips) });
        }
    }

    return matching_problem_type{ .detectors=std::move(detectors), .edges=std::move(edges) };
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type
BLOSSOMV::solve_matching_problem(const matching_problem_type& mp) const
{
    using assignment_type = MATCHING_DATA::assignment_type;

    // map a detector id to its matching-problem node index
    std::unordered_map<det_id_type, size_t> idx_map;
    idx_map.reserve(mp.detectors.size());
    for (size_t i = 0; i < mp.detectors.size(); i++)
        idx_map[mp.detectors[i]] = i;

    const size_t n = mp.detectors.size(),
                 m = mp.edges.size();

    b5::PerfectMatching pm(n, m);
    pm.options.verbose = false;
    for (size_t k = 0; k < m; k++)
    {
        const auto& e = mp.edges[k];
        [[ maybe_unused ]] auto _k = pm.AddEdge(idx_map[e.d1], idx_map[e.d2], e.w_qu);
        assert(k == static_cast<size_t>(_k));
    }

    pm.Solve();

    // Apply the frame flips of every matched edge.
    result_type out{.flipped_obs=obs_type(num_observables)};
    for (size_t k = 0; k < m; k++)
    {
        if (pm.GetSolution(k))
        {
            const auto& e = mp.edges[k];
            out.flipped_obs ^= e.frame_flips;
            // update matching data:
            out.matching_data.add(e);
        }
    }

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder
