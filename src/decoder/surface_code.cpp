/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 *  contributor: OpenAI Codex
 * */

#include "decoder/surface_code.h"

#include <pymatching/sparse_blossom/driver/user_graph.h>

#include <PerfectMatching.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

/*
 * Builds a DEM that turns the logical observable into an *explicit boundary node*: each
 * observable target `Lk` is folded into a detector node at index `num_detectors + k`.
 * pymatching represents observables as edge bit-masks, not graph nodes, so there is
 * otherwise nothing to flip; this node is the handle. Because the observable is carried
 * only by boundary (single-detector) edges in a surface-code memory DEM, folding it moves
 * exactly those edges onto the node, and its parity then equals the logical class:
 * decoding with the node unfired vs. fired enumerates the two parity classes (see
 * `PyMatching::decode`).
 *
 * Two structural properties of the memory DEM make this correct and safe:
 *   - Every observable-crossing component is single-detector, so the fold yields clean
 *     two-detector edges (never a three-detector hyperedge that MWPM cannot match). The
 *     assert guards other code families where that may not hold.
 *   - Only the observable-crossing boundary edges are rerouted onto the node; all other
 *     boundary edges (the opposite boundary and the entire complementary-basis subgraph)
 *     are left untouched, so the virtual boundary still exists and odd-parity syndromes
 *     remain matchable. The complementary-basis subgraph never touches the node, so it
 *     contributes identical weight to both decode passes and cancels in the gap.
 * */
stim::DetectorErrorModel
_build_gap_dem(const stim::DetectorErrorModel& dem)
{
    const uint64_t num_detectors = dem.count_detectors();
    stim::DetectorErrorModel out;
    std::vector<stim::DemTarget> targets;
    dem.iter_flatten_error_instructions(
            [&] (const stim::DemInstruction& inst)
            {
                targets.clear();
                bool first_group = true;
                inst.for_separated_targets(
                        [&] (const auto& grp)
                        {
                            if (!first_group)
                                targets.push_back(stim::DemTarget::separator());
                            first_group = false;

                            size_t num_dets = 0;
                            for (const auto& t : grp)
                            {
                                if (t.is_observable_id())
                                {
                                    targets.push_back(stim::DemTarget::relative_detector_id(num_detectors + t.val()));
                                    num_dets++;
                                }
                                else
                                {
                                    if (t.is_relative_detector_id())
                                        num_dets++;
                                    targets.push_back(t);
                                }
                            }
                            [[ maybe_unused ]] const bool graphlike = (num_dets <= 2);
                            assert(graphlike && "folding observable into a detector produced a hyperedge");
                        });
                out.append_error_instruction(inst.arg_data[0], targets, "");
            });
    return out;
}

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PyMatching::PyMatching(const stim::DetectorErrorModel& dem)
    :PyMatching(dem, config_type{})
{
}

PyMatching::PyMatching(const stim::DetectorErrorModel& dem, config_type conf)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    conf(conf),
    obs_det_id_(dem.count_detectors()),
    mwpm_(pm::detector_error_model_to_mwpm(
                conf.enable_gap_estimation ? _build_gap_dem(dem) : dem,
                pm::NUM_DISTINCT_WEIGHTS))
{
    if (conf.enable_gap_estimation)
        assert(num_observables == 1 && "complementary gap estimation supports exactly one observable");
    if (conf.report_matching_data && num_observables > sizeof(pm::obs_int) * 8)
        throw std::invalid_argument("PyMatching matched-edge reporting supports at most 64 observables");

    norm_const_ = mwpm_.flooder.graph.normalising_constant;
    if (norm_const_ <= 0.0)
        norm_const_ = 1.0;
    // An integer edge weight equals `round(-ln(p/(1-p)) * norm_const_)`, so decibels are
    // `10 * log10(...) = (10/ln(10)) * (-ln(p/(1-p)))`, i.e. `10/(ln(10)*norm_const_)` per
    // unit of quantized weight.
    decibels_per_w_ = 10.0 / (std::log(10.0) * norm_const_);
}

Result
PyMatching::decode(SyndromeRef syn, ObsRef obs)
{
    const auto start = std::chrono::steady_clock::now();
    auto record_latency = [&]
    {
        const auto elapsed_ns = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - start).count());
        ++stats_.timed_shots;
        stats_.decode_total_ns += elapsed_ns;
        stats_.decode_min_ns = std::min(stats_.decode_min_ns, elapsed_ns);
        stats_.decode_max_ns = std::max(stats_.decode_max_ns, elapsed_ns);
        return elapsed_ns;
    };
    pm::total_weight_int w_primary = 0;
    auto out = internal_decode(syn, /*fire_obs_det=*/false, w_primary,
                               conf.report_matching_data && !conf.enable_gap_estimation);
    if (!conf.enable_gap_estimation)
    {
        out.decode_latency_ns = record_latency();
        return out;
    }

    // Complementary pass: fire the observable detector to force the opposite parity class.
    pm::total_weight_int w_complement = 0;
    if (conf.report_matching_data)
    {
        // The endpoint API does not return solution weight; retain the established weight
        // passes for gap estimation, then report endpoints for the selected logical class.
        internal_decode(syn, /*fire_obs_det=*/true, w_complement);
        const bool predict_flip = (w_complement < w_primary);
        pm::total_weight_int ignored_weight = 0;
        out = internal_decode(syn, predict_flip, ignored_weight, true);
    }
    else
        internal_decode(syn, /*fire_obs_det=*/true, w_complement);

    // set output observable:
    const bool predict_flip = (w_complement < w_primary);
    out.flipped_obs = ObsType(num_observables);
    out.flipped_obs[0] = predict_flip ? 1 : 0;
    if (conf.report_matching_data)
        out.matching_data.frame_flips = out.flipped_obs;

    // update gap distributions
    if (predict_flip)
        std::swap(w_primary, w_complement);
    double g = (w_complement-w_primary) * decibels_per_w_;

    // negate gap if this is an error:
    s_unsigned_gap.add(g);
    if (obs[0] != out.flipped_obs[0])
    {
        s_unsigned_gap_errors.add(g);
        g = -g;
    }
    out.matching_data.gap = g;
    s_signed_gap.add(g);

    out.decode_latency_ns = record_latency();
    return out;
}

Result
PyMatching::internal_decode(SyndromeRef syn, bool fire_obs_det, pm::total_weight_int& weight_out,
                            bool report_matching_data)
{
    // Collect indices of fired detectors.
    std::vector<uint64_t> det_events;
    for (size_t i = 0; i < num_detectors; i++)
        if (syn[i])
            det_events.push_back(i);
    // The observable detector has the largest index, so appending keeps `det_events` sorted.
    if (fire_obs_det)
        det_events.push_back(obs_det_id_);

    // Decode. `pm::decode_detection_events` writes one *byte* per observable (observable `k`'s
    // flip goes in `obs_begin_ptr[k]`), which is a different layout from `ObsType`'s bit-packed
    // storage (where `flipped_obs[k]` is bit `k`) -- so we decode into a byte buffer first and
    // then translate it into the bit-packed `ObsType` that the rest of the codebase expects.
    Result res;
    res.flipped_obs = ObsType(num_observables);
    if (report_matching_data)
    {
        pm::decode_detection_events_to_match_edges(mwpm_, det_events);
        const auto* node_begin = mwpm_.flooder.graph.nodes.data();
        res.matching_data.frame_flips = ObsType(num_observables);
        for (const auto& edge : mwpm_.flooder.match_edges)
        {
            pm::obs_int edge_obs_mask = edge.obs_mask;
            if (conf.enable_gap_estimation &&
                ((edge.loc_from && static_cast<size_t>(edge.loc_from - node_begin) == obs_det_id_) ||
                 (edge.loc_to && static_cast<size_t>(edge.loc_to - node_begin) == obs_det_id_)))
                edge_obs_mask ^= 1;
            for (size_t k = 0; k < num_observables; k++)
                if ((edge_obs_mask >> k) & 1)
                    res.flipped_obs[k] = !res.flipped_obs[k];
            const auto id = [node_begin, this](const pm::DetectorNode* node) -> int64_t
            {
                if (!node) return -1;
                const auto index = static_cast<size_t>(node - node_begin);
                return index == obs_det_id_ ? -1 : static_cast<int64_t>(index);
            };
            const int64_t d1 = id(edge.loc_from), d2 = id(edge.loc_to);
            if (d1 < 0 && d2 < 0) continue;
            MatchingData::assignment_type assignment{
                .d1 = d1,
                .d2 = d2,
                .pr = 1.0,
                .w_qu = 0,
                .frame_flips = ObsType(num_observables)};
            for (size_t k = 0; k < num_observables; k++)
                assignment.frame_flips[k] = (edge_obs_mask >> k) & 1;
            res.matching_data.add(assignment);
        }
        weight_out = 0;
    }
    else
    {
        std::vector<uint8_t> obs_bytes(num_observables, 0);
        pm::total_weight_int weight = 0;
        pm::decode_detection_events(mwpm_, det_events, obs_bytes.data(), weight, false);
        weight_out = weight;
        for (size_t i = 0; i < num_observables; i++)
            res.flipped_obs[i] = obs_bytes[i] & 1;
        res.matching_data.total_weight = weight;
    }
    return res;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PyMatching::Statistics
PyMatching::statistics() const
{
    auto stats = stats_;
    if (!stats.timed_shots)
        stats.decode_min_ns = 0;
    return stats;
}

void
PyMatching::mpi_accumulate()
{
#if defined(ENABLE_MPI)
    uint64_t sums[]{stats_.timed_shots, stats_.decode_total_ns};
    MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    stats_.timed_shots = sums[0];
    stats_.decode_total_ns = sums[1];
    MPI_Allreduce(MPI_IN_PLACE, &stats_.decode_min_ns, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &stats_.decode_max_ns, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
#endif
    s_signed_gap.mpi_accumulate();
    s_unsigned_gap.mpi_accumulate();
    s_unsigned_gap_errors.mpi_accumulate();
}

void
PyMatching::print_stats(std::ostream& ostrm) const
{
    const auto stats = statistics();
    const double average_us = stats.timed_shots
        ? double(stats.decode_total_ns) / double(stats.timed_shots) / 1000.0 : 0.0;
    ostrm << "PyMatching decode latency (us): timed_shots=" << stats.timed_shots
          << " average=" << average_us
          << " min=" << double(stats.decode_min_ns) / 1000.0
          << " max=" << double(stats.decode_max_ns) / 1000.0
          << " includes_empty_syndromes=true scope=full_history\n";
    // The gap histogram is only populated when complementary-gap estimation is enabled;
    // otherwise it is empty (and its mean would divide by a zero count).
    if (conf.enable_gap_estimation)
    {
        ostrm << s_signed_gap.to_string_full()
                << "\n" << s_unsigned_gap.to_string_full()
                << "\n" << s_unsigned_gap_errors.to_string_full()
                << "\nnormalization constant = " << norm_const_ 
                << "\ndecibels per weight = " << decibels_per_w_
                << "\n";
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

using DetIdType = BlossomV::det_id_type;
using AdjEntryType = BlossomV::adj_entry_type;
using AdjListType = BlossomV::adj_list_type;
using MwpmEdgeType = BlossomV::mwpm_edge_type;
using MatchingProblemType = BlossomV::matching_problem_type;

constexpr DetIdType BOUNDARY_ID{-1};

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
_update_adjacency_list(AdjListType& adj, DetIdType d, double p, ObsRef frame_flips)
{
    auto adj_it = std::find_if(adj.begin(), adj.end(),
                        [d] (const auto& e) { return e.d == d; });
    if (adj_it != adj.end())
        adj_it->pr = (1-adj_it->pr)*p + (1-p)*adj_it->pr;
    else
        adj.push_back(AdjEntryType{d, p, ObsType{frame_flips}});
}

/*
 * Data structures used when computing pairwise distances via Dijkstra. The weight `w`
 * accumulates the raw (unquantized) `-log(pr)` of each edge along the path; it is only
 * quantized once, after the shortest path is finalized.
 * */
struct distance_type
{
    double   w = std::numeric_limits<double>::max();
    ObsType frame_flips;
};

struct distance_queue_entry
{
    DetIdType d;
    double      w;
};

struct distance_cmp
{
    bool operator()(const distance_queue_entry& a, const distance_queue_entry& b) const { return a.w > b.w; }
};

using DistanceQueueType = std::priority_queue<distance_queue_entry,
                                                std::vector<distance_queue_entry>,
                                                distance_cmp>;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

BlossomV::BlossomV(const stim::DetectorErrorModel& dem)
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
                            std::vector<DetIdType> dets;
                            ObsType frame_flips(num_observables);
                            for (const auto& t : grp)
                            {
                                if (t.is_relative_detector_id())
                                    dets.push_back(static_cast<DetIdType>(t.val()));
                                else if (t.is_observable_id())
                                    frame_flips[t.val()] ^= 1;
                            }

                            DetIdType d1 = dets[0],
                                        d2 = (dets.size() == 1) ? BOUNDARY_ID : dets[1];
                            _update_adjacency_list(adj_matrix_[d1], d2, pr, frame_flips);
                            if (d2 == BOUNDARY_ID)
                                _update_adjacency_list(boundary_adjacency_, d1, pr, frame_flips);
                            else
                                _update_adjacency_list(adj_matrix_[d2], d1, pr, frame_flips);
                        });
            });
}

const AdjListType&
BlossomV::adj_matrix(DetIdType d) const
{
    return (d == BOUNDARY_ID) ? boundary_adjacency_ : adj_matrix_[d];
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Result
BlossomV::decode(SyndromeRef syndrome, ObsRef)
{
    auto detectors = collect_detection_events(syndrome);
    if (detectors.empty())
        return Result{.flipped_obs=ObsType(num_observables)};

    auto mp = synthesize_matching_problem(std::move(detectors));
    return solve_matching_problem(mp);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::vector<DetIdType>
BlossomV::collect_detection_events(SyndromeRef syndrome) const
{
    std::vector<DetIdType> detectors;
    for (DetIdType i = 0; i < static_cast<DetIdType>(num_detectors); i++)
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

MatchingProblemType
BlossomV::synthesize_matching_problem(std::vector<DetIdType> detectors) const
{
    const size_t n = detectors.size();

    // map a detector id to its slot in `dist` (boundary lives at index `num_detectors`)
    auto slot = [this] (DetIdType d) -> size_t
                {
                    return (d == BOUNDARY_ID) ? num_detectors : static_cast<size_t>(d);
                };

    // Compute pairwise distances via Dijkstra over the full decoding graph: run from
    // each detector and read off the upper triangle (distances are symmetric).
    const distance_type fill_val{.frame_flips=ObsType(num_observables)};
    std::vector<distance_type> dist(num_detectors+1, fill_val);

    std::vector<MwpmEdgeType> edges;
    edges.reserve(n*(n-1)/2);

    for (size_t ii = 0; ii+1 < n; ii++)
    {
        const DetIdType d1 = detectors[ii];

        std::fill(dist.begin(), dist.end(), fill_val);
        dist[slot(d1)].w = 0.0;
        DistanceQueueType pq;
        pq.push({d1, 0.0});
        while (pq.size() > 0)
        {
            auto e = std::move(pq.top());
            pq.pop();
            const DetIdType z1 = e.d;
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
            const DetIdType d2 = detectors[jj];
            const size_t j = slot(d2);
            // `dist[j].w` is the accumulated `-log` path weight: the matching's error
            // probability is `exp(-w)`, and the quantized weight is `_quantize(w)`.
            edges.push_back(MwpmEdgeType{ .d1=d1,
                                            .d2=d2,
                                            .pr=std::exp(-dist[j].w),
                                            .w_qu=_quantize(dist[j].w),
                                            .frame_flips=std::move(dist[j].frame_flips) });
        }
    }

    return MatchingProblemType{ .detectors=std::move(detectors), .edges=std::move(edges) };
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Result
BlossomV::solve_matching_problem(const MatchingProblemType& mp) const
{
    using AssignmentType = MatchingData::assignment_type;

    // map a detector id to its matching-problem node index
    std::unordered_map<DetIdType, size_t> idx_map;
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
    Result out{.flipped_obs=ObsType(num_observables)};
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

} // namespace dec
