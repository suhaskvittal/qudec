/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 *  contributor: OpenAI Codex
 * */

#include "decoder/pymatching.h"

#include <pymatching/sparse_blossom/driver/user_graph.h>


#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

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
                pm::NUM_DISTINCT_WEIGHTS,
                conf.correlated || conf.report_matching_data, conf.correlated)),
    cluster_gap_graph_(conf.enable_cluster_gap ? num_detectors : 0)
{
    if (conf.enable_cluster_gap)
    {
        if (conf.enable_gap_estimation || conf.correlated)
            throw std::invalid_argument("Cluster gap requires ordinary uncorrelated PyMatching; disable --gap and --correlated");
        validate_cluster_gap_dem(dem);
        const auto* begin = mwpm_.flooder.graph.nodes.data();
        for (size_t u = 0; u < mwpm_.flooder.graph.nodes.size(); u++)
        {
            const auto& node = mwpm_.flooder.graph.nodes[u];
            for (size_t k = 0; k < node.neighbors.size(); k++)
            {
                const auto* other = node.neighbors[k];
                const size_t v = other ? other - begin :
                    ((node.neighbor_observables[k] & 1) ? cluster_gap_graph_.boundary_one
                                                      : cluster_gap_graph_.boundary_zero);
                if (other && u > v)
                    continue;
                cluster_gap_graph_.add_edge(u, v, node.neighbor_weights[k]);
            }
        }
        if (mwpm_.flooder.negative_weight_sum != 0)
            throw std::invalid_argument("Cluster gap does not support negative matching weights");
        cluster_gap_graph_.finish();
    }
    if (conf.correlated && (conf.enable_gap_estimation || conf.report_matching_data))
        throw std::invalid_argument("Correlated PyMatching currently supports observable decoding only; disable gap estimation and matched-endpoint reporting");
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
        if (conf.enable_cluster_gap && obs[0] != out.flipped_obs[0])
            s_cluster_gap_errors.add(out.cluster_gap_db);
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
    std::function<void(const pm::Mwpm&)> before_extraction;
    if (conf.enable_cluster_gap)
        before_extraction = [&](const pm::Mwpm& matching)
        {
            std::vector<std::pair<size_t, double>> sources;
            for (size_t v : det_events)
            {
                if (!cluster_gap_graph_.relevant[v])
                    continue;
                double radius = 0;
                auto* region = matching.flooder.graph.nodes[v].region_that_arrived;
                while (region)
                {
                    radius += region->radius.get_distance_at_time(matching.flooder.queue.cur_time);
                    region = region->blossom_parent;
                }
                if (radius < 0)
                    throw std::logic_error("MWPM cluster radius was negative");
                sources.push_back({v, radius});
            }
            res.cluster_gap_db = cluster_gap_graph_.gap(cluster_gap_graph_.ball_coverage(sources))
                                 * decibels_per_w_;
            res.has_cluster_gap = true;
            s_cluster_gap.add(res.cluster_gap_db);
        };
    if (report_matching_data)
    {
        pm::decode_detection_events_to_match_edges(mwpm_, det_events, before_extraction);
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
            if (!conf.enable_gap_estimation)
            {
                const auto& search_nodes = mwpm_.search_flooder.graph.nodes;
                const size_t from = static_cast<size_t>(edge.loc_from - node_begin);
                const size_t to = edge.loc_to ? static_cast<size_t>(edge.loc_to - node_begin) : SIZE_MAX;
                mwpm_.search_flooder.iter_edges_on_shortest_path_from_middle(
                    from, to, [&](const pm::SearchGraphEdge& path_edge)
                    {
                        const auto* source = path_edge.detector_node;
                        const auto* neighbor = source->neighbors[path_edge.neighbor_index];
                        uint64_t observable_mask = 0;
                        for (size_t observable : source->neighbor_observable_indices[path_edge.neighbor_index])
                            observable_mask ^= uint64_t(1) << observable;
                        assignment.correction_edges.push_back({
                            static_cast<int64_t>(source - search_nodes.data()),
                            neighbor ? static_cast<int64_t>(neighbor - search_nodes.data()) : -1,
                            observable_mask});
                    });
            }
            res.matching_data.add(assignment);
        }
        weight_out = 0;
    }
    else
    {
        std::vector<uint8_t> obs_bytes(num_observables, 0);
        pm::total_weight_int weight = 0;
        pm::decode_detection_events(mwpm_, det_events, obs_bytes.data(), weight, conf.correlated, before_extraction);
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
    s_cluster_gap.mpi_accumulate();
    s_cluster_gap_errors.mpi_accumulate();
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
          << " includes_empty_syndromes=true scope=full_history"
          << " correlated=" << (conf.correlated ? "true" : "false") << '\n';
    if (conf.enable_cluster_gap)
    {
        ostrm << "Cluster gap (dB): count=" << stats.timed_shots;
        if (stats.timed_shots)
            ostrm << " " << s_cluster_gap.to_string_some();
        ostrm << "\n";
        if (conf.print_cluster_gap_histogram)
            ostrm << s_cluster_gap.to_string_full() << "\n"
                  << s_cluster_gap_errors.to_string_full() << "\n";
    }
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

} // namespace dec
