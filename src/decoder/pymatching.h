/* 
 *  author: Suhas Vittal
 *  date:   16 May 2026
 *  contributor: OpenAI Codex
 * */

#ifndef DECODER_PYMATCHING_h
#define DECODER_PYMATCHING_h

#include "decoder/types.h"
#include "decoder/cluster_gap.h"
#include "stats.h"

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>

#include <iosfwd>
#include <limits>
#include <vector>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class PyMatching
{
public:
    struct config_type
    {
        bool enable_gap_estimation{false};
        // Report matched detector endpoints and each pair's observable parity.
        // Boundary endpoints use -1. Weights are not populated in this mode.
        bool report_matching_data{false};
        // Use PyMatching's existing two-pass correlated MWPM backend.
        bool correlated{false};
        bool enable_cluster_gap{false};
        bool print_cluster_gap_histogram{false};
    };

    // Full-value cluster output; histogram bins are printed only on request.
    StatsHistogram<double> s_cluster_gap{"CLUSTER_GAP_DB", 0, 256, 256};
    StatsHistogram<double> s_cluster_gap_errors{"CLUSTER_GAP_DB_ERRORS", 0, 256, 256};

    struct Statistics
    {
        uint64_t timed_shots{0};
        uint64_t decode_total_ns{0};
        uint64_t decode_min_ns{std::numeric_limits<uint64_t>::max()};
        uint64_t decode_max_ns{0};
    };

    const size_t num_detectors;
    const size_t num_observables;

    const config_type conf;

    /*
     * `s_gap` is the complementary gap distribution.
     * `s_gap_error` are logical errors at a given complementary gap.
     * */
    StatsHistogram<double> s_signed_gap{"SIGNED_GAP", -256, 256, 256},
                           s_unsigned_gap{"UNSIGNED_GAP", 0, 256, 256},
                           s_unsigned_gap_errors{"UNSIGNED_GAP_ERRORS", 0, 256, 256};
private:
    /*
     * When `conf.enable_gap_estimation` is set, the matching graph is built from an
     * augmented DEM in which the single logical observable is folded into an explicit
     * boundary node at index `obs_det_id_` (= `num_detectors`). The node's parity equals
     * the logical class, so decoding with it unfired vs. fired yields the two parity
     * classes and their weight margin is the gap (see `_build_gap_dem` in the .cpp).
     * */
    const size_t obs_det_id_;

    pm::Mwpm mwpm_;
    ClusterGapGraph cluster_gap_graph_;
    Statistics stats_;

    /*
     * Weight->decibel conversion factor read from the matching graph's normalising
     * constant `C`: an integer edge weight equals `round(-ln(p/(1-p)) * C)`, so
     * `decibels_per_w_ = 10 / (ln(10) * C)` and a natural-log weight is `w_qu / C`.
     * */
    double norm_const_{1.0};
    double decibels_per_w_{1.0};
public:
    explicit PyMatching(const stim::DetectorErrorModel&);
    PyMatching(const stim::DetectorErrorModel&, config_type);

    Result decode(SyndromeRef, ObsRef);

    Statistics statistics() const;
    void print_stats(std::ostream&) const;
    void mpi_accumulate();
private:
    Result internal_decode(SyndromeRef, bool fire_obs_det, pm::total_weight_int& weight_out,
                                bool report_matching_data=false);
};


} // namespace dec

#endif // DECODER_PYMATCHING_h
