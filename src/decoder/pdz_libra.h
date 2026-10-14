/*
 *  contributor: OpenAI Codex
 *  purpose: Compose uncorrelated PyMatching, PDZ and Libra.
 * */

#ifndef DECODER_PDZ_LIBRA_h
#define DECODER_PDZ_LIBRA_h

#include "decoder/libra.h"
#include "decoder/pdz.h"
#include "decoder/surface_code.h"
#include "preprocessor/pdz.h"
#include "stats.h"

#include <array>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// One instance per worker. The DEM must outlive the instance, as required by PDZ.
class PDZLibra
{
public:
    struct Statistics
    {
        uint64_t sampled_shots;
        uint64_t sampled_failures;
        uint64_t input_events;
        uint64_t residual_events;
        uint64_t bypassed_shots;
        uint64_t libra_called_shots;
        double libra_called_fraction;
        double libra_mean_latency_us;
        double libra_weighted_latency_us;
    };

    // X is the number of detector bits removed by PDZ / 2. Unit-width bin k
    // covers 2k or 2k+1 removed bits: boundary matches can remove one detector.
    StatsHistogram<double> s_pdz_removed_pairs;
    // Restricted to shots that enter Libra's ensemble.
    StatsHistogram<double> s_ensemble_by_pdz_removed_pairs;
    // Same bins, restricted to ensemble shots whose final pipeline frame
    // differs from the initial uncorrelated PyMatching frame.
    StatsHistogram<double> s_ensemble_frame_changes_by_pdz_removed_pairs;

    // The same three populations, using only removed memory-basis detector
    // bits / 2. Opposing-basis pairs are excluded; boundary matches still
    // contribute half a pair because they remove one syndrome bit.
    StatsHistogram<double> s_pdz_same_basis_removed_pairs;
    StatsHistogram<double> s_ensemble_by_pdz_same_basis_removed_pairs;
    StatsHistogram<double> s_ensemble_frame_changes_by_pdz_same_basis_removed_pairs;

    // Unit-width bin k corresponds to exactly k fired detectors remaining
    // after PDZ, including opposing-basis detectors. Count all shots first,
    // then the subset that enters Libra's ensemble.
    StatsHistogram<double> s_post_pdz_hamming_weight;
    StatsHistogram<double> s_ensemble_by_post_pdz_hamming_weight;
    // Ensemble shots whose final pipeline frame differs from initial PyMatching,
    // grouped by the same post-PDZ Hamming-weight bins.
    StatsHistogram<double> s_ensemble_frame_changes_by_post_pdz_hamming_weight;

    // At least removal_threshold of fired bits removed => return initial
    // PyMatching unchanged. Empty residuals always bypass; 1.0 requires complete
    // removal and therefore disables bypass for nonempty residual syndromes.
    PDZLibra(const stim::DetectorErrorModel&, pp::PDZ::config_type, Libra::Config,
             double removal_threshold = 0.5);
    ~PDZLibra();

    // Optional per-worker gzip recording. Start before decoding, and finish
    // explicitly to check the gzip close status. Existing files are rejected.
    void record_syndromes(const std::string& path, uint32_t rank);
    void finish_recording();

    Result decode(SyndromeRef, ObsRef);
    Statistics statistics() const;
    void print_stats(std::ostream&) const;
    void mpi_accumulate();
private:
    struct SyndromeRecorder;
    std::unique_ptr<SyndromeRecorder> recorder_;
    PyMatching initial_;
    Libra residual_;
    PDZDecoder<PyMatching, Libra> pipeline_;
    const double removal_threshold_;
    std::array<uint64_t, 5> counts_{};
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec

#endif
