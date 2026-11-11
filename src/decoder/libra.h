/*
 *  contributor: OpenAI Codex
 *  date:        13 October 2026
 *  purpose:     Expose the Libra matching-synthesis decoder.
 * */

#ifndef QUDEC_DECODER_LIBRA_H
#define QUDEC_DECODER_LIBRA_H

#include "decoder/types.h"

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <limits>
#include <vector>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Matching synthesis, C. Jones, arXiv:2408.12135 (Sections III--V).
 *
 * Libra requires a single logical observable and graphlike components in each
 * CSS basis. Each instance owns mutable matching workspaces: use one instance
 * per worker. Returned error ids refer to the canonicalized error channels,
 * which merge DEM instructions with identical detector and observable action.
 * */
class Libra
{
public:
    struct Config
    {
        size_t ensemble_size = 60;
        double gap_threshold_db = 20;
        size_t positive_cycle_capacity = 30;
        uint64_t seed = 0x4c49425241ULL;
        bool degeneracy = true;
        bool validate = false;
        bool profile = false;
        bool latency_histogram = false;
    };

    struct DetailedResult
    {
        bool prediction = false;
        bool ensemble_used = false;
        bool fallback_used = false;
        double gap_db = 0;
        double weights[2]{};
        double log_masses[2]{};
        std::vector<uint32_t> errors;
    };

    struct Dimensions
    {
        size_t detectors, hyperedges, graph_edges, correlation_rules;
    };

    struct Statistics
    {
        uint64_t shots = 0;
        uint64_t ensemble_shots = 0;
        uint64_t fallback_shots = 0;
        uint64_t improving_cycles = 0;
        uint64_t complementary_ns = 0;
        uint64_t ensemble_ns = 0;
        uint64_t synthesis_ns = 0;
        uint64_t degeneracy_ns = 0;
        uint64_t fallback_ns = 0;
        // Full wall time of public decode() calls, including empty syndromes.
        // decode_detailed() calls are not included in these timing counters.
        uint64_t timed_shots = 0;
        uint64_t decode_total_ns = 0;
        uint64_t decode_min_ns = std::numeric_limits<uint64_t>::max();
        uint64_t decode_max_ns = 0;
    };

    const size_t num_detectors;
    const size_t num_observables;

    explicit Libra(const stim::DetectorErrorModel&);
    Libra(const stim::DetectorErrorModel&, Config);
    ~Libra();
    Libra(Libra&&) noexcept;
    Libra(const Libra&) = delete;
    Libra& operator=(const Libra&) = delete;

    /*
     * ObsRef is ground truth for the generic sampler API; prediction ignores it.
     * matching_data.gap is the unsigned initial complementary gap in dB.
     * matching_data.probability is an unnormalized class score, not confidence.
     * */
    Result decode(SyndromeRef, ObsRef);
    // Contributor: OpenAI Codex. The seed decoder has already selected this
    // syndrome for ensembling, independently of Libra's internal gap.
    Result decode_with_ensemble(SyndromeRef, ObsRef);
    DetailedResult decode_detailed(SyndromeRef);

    // Prepare immutable perturbed priors before timing decoding.
    void prepare_ensemble();

    Dimensions dimensions() const;
    Statistics statistics() const;
    void print_stats(std::ostream&) const;
    void mpi_accumulate();

    // Reconstruct a returned error assignment for correctness checks.
    std::pair<SyndromeType, bool> syndrome_of(const std::vector<uint32_t>&) const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec

#endif // QUDEC_DECODER_LIBRA_H
