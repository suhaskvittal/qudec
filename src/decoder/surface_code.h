/* 
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#ifndef DECODER_SURFACE_CODE_h
#define DECODER_SURFACE_CODE_h

#include "decoder/common.h"
#include "stats.h"

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>

#include <iosfwd>
#include <vector>

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class PyMatching
{
public:
    const size_t num_detectors;
    const size_t num_observables;

    const bool estimate_complementary_gap;

    /*
     * `s_gap` is the complementary gap distribution.
     * `s_gap_error` are logical errors at a given complementary gap.
     * */
    StatsHistogram<double> s_signed_gap{"SIGNED_GAP", -256, 256, 256},
                           s_unsigned_gap{"UNSIGNED_GAP", 0, 256, 256},
                           s_unsigned_gap_errors{"UNSIGNED_GAP_ERRORS", 0, 256, 256};
private:
    /*
     * When `estimate_complementary_gap` is set, the matching graph is built from an
     * augmented DEM in which the single logical observable is folded into an explicit
     * boundary node at index `obs_det_id_` (= `num_detectors`). The node's parity equals
     * the logical class, so decoding with it unfired vs. fired yields the two parity
     * classes and their weight margin is the gap (see `_build_gap_dem` in the .cpp).
     * */
    const size_t obs_det_id_;

    pm::Mwpm mwpm_;

    /*
     * Weight->decibel conversion factor read from the matching graph's normalising
     * constant `C`: an integer edge weight equals `round(-ln(p/(1-p)) * C)`, so
     * `decibels_per_w_ = 10 / (ln(10) * C)` and a natural-log weight is `w_qu / C`.
     * */
    double norm_const_{1.0};
    double decibels_per_w_{1.0};
public:
    PyMatching(const stim::DetectorErrorModel&, bool enable_gap_estimation);

    result_type decode(SyndromeRef, ObsRef);

    void print_stats(std::ostream&) const;
    
    void 
    mpi_accumulate()
    {
        s_signed_gap.mpi_accumulate();
        s_unsigned_gap.mpi_accumulate();
        s_unsigned_gap_errors.mpi_accumulate();
    }
private:
    result_type internal_decode(SyndromeRef, bool fire_obs_det, pm::total_weight_int& weight_out);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * `BlossomV` is a global Minimum-Weight Perfect Matching (MWPM) decoder built on
 * top of Blossom-5. Unlike `ClusterMatch`, it does not split the syndrome into
 * clusters: it matches all detection events at once. As such, it is an exact
 * software reference / ground-truth decoder whose accuracy should track `PyMatching`.
 *
 * The decoder works in three steps:
 *   (1) The constructor encaches the decoding graph (adjacency matrix + boundary
 *       adjacency) from the detector error model.
 *   (2) On each syndrome, pairwise distances between all detection events are
 *       computed via Dijkstra over the full decoding graph.
 *   (3) Blossom-5 computes the MWPM and the Pauli frame flips of the matched edges
 *       are applied as the correction.
 * */
class BlossomV
{
public:
    using det_id_type = int64_t;

    /*
     * Decoding graph edge: the other detector `d`, its probability `pr`, and the
     * Pauli frame flips `frame_flips` associated with the edge.
     * */
    struct adj_entry_type
    {
        det_id_type d;
        double      pr;
        ObsType    frame_flips;
    };

    using adj_list_type = std::vector<adj_entry_type>;

    /*
     * An edge in the matching problem: both endpoints, a quantized weight, and the
     * accumulated Pauli frame flips along the shortest path between the endpoints.
     * */
    using mwpm_edge_type = MatchingData::assignment_type;

    /*
     * A matching problem: the detectors to be matched and the complete graph of
     * edges between them.
     * */
    struct matching_problem_type
    {
        std::vector<det_id_type>    detectors;
        std::vector<mwpm_edge_type> edges;
    };

    const size_t num_detectors;
    const size_t num_observables;
public:
    BlossomV(const stim::DetectorErrorModel&);

    result_type decode(SyndromeRef, ObsRef);

    void print_stats(std::ostream&) const {}
    void mpi_accumulate() {}
private:
    /*
     * `collect_detection_events()` gathers all flipped detectors, appending the
     * boundary when their count is odd so that a perfect matching exists.
     * */
    std::vector<det_id_type> collect_detection_events(SyndromeRef) const;

    /*
     * `synthesize_matching_problem()` computes the pairwise distances between all
     * detection events (Dijkstra over the full decoding graph).
     * */
    matching_problem_type synthesize_matching_problem(std::vector<det_id_type>) const;

    /*
     * `solve_matching_problem()` runs Blossom-5 and returns the correction implied
     * by the min-weight perfect matching.
     * */
    result_type solve_matching_problem(const matching_problem_type&) const;

    /*
     * Decoding graph adjacency data:
     * */
    std::vector<adj_list_type> adj_matrix_;
    adj_list_type boundary_adjacency_;

    const adj_list_type& adj_matrix(det_id_type) const;
};

} // namespace decoder

#endif // DECODER_SURFACE_CODE_h
