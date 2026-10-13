#ifndef QUDEC_DECODER_LIBRA_VARIANT_H
#define QUDEC_DECODER_LIBRA_VARIANT_H

#include "decoder/common.h"
#include <memory>
#include <cstdint>
#include <iosfwd>
#include <vector>

namespace dec {

// Configurable matching and fallback ablations of Libra.
// Matching synthesis, C. Jones, arXiv:2408.12135 (Sections III--V).
// A single logical observable and graphlike components in each CSS basis
// are required. Instances own mutable matching workspaces: one per worker.
class LibraVariant {
public:
    struct Config {
        size_t ensemble_size = 60;
        double gap_threshold_db = 20;
        size_t positive_cycle_capacity = 30;
        uint64_t seed = 0x4c49425241ULL;
        bool degeneracy = true;
        bool validate = false;
        bool profile = false;
        bool correlated_matching = true;
        bool large_gap_fallback = true;
    };
    struct DetailedResult {
        bool prediction = false;
        bool ensemble_used = false;
        bool fallback_used = false;
        double gap_db = 0;
        // Prediction of the cheaper initial constrained representative.
        bool initial_prediction = false;
        double weights[2]{};
        double log_masses[2]{};
        std::vector<uint32_t> errors;
    };
    struct Dimensions {
        size_t detectors, hyperedges, graph_edges, correlation_rules;
    };
    struct Statistics {
        uint64_t shots = 0, ensemble_shots = 0, fallback_shots = 0, improving_cycles = 0;
        uint64_t complementary_ns = 0, ensemble_ns = 0, synthesis_ns = 0, degeneracy_ns = 0;
        uint64_t fallback_ns = 0;
    };
    const size_t num_detectors;
    const size_t num_observables;

    explicit LibraVariant(const stim::DetectorErrorModel&);
    LibraVariant(const stim::DetectorErrorModel&, Config);
    ~LibraVariant();
    LibraVariant(LibraVariant&&) noexcept;
    LibraVariant(const LibraVariant&) = delete;
    LibraVariant& operator=(const LibraVariant&) = delete;

    result_type decode(SyndromeRef, ObsRef);
    DetailedResult decode_detailed(SyndromeRef);
    // Prepare immutable perturbed priors before timing decoding.
    void prepare_ensemble();
    Dimensions dimensions() const;
    Statistics statistics() const;
    void print_stats(std::ostream&) const;
    void mpi_accumulate() {}
    // Reconstruct a returned error assignment for correctness checks.
    std::pair<SyndromeType, bool> syndrome_of(const std::vector<uint32_t>&) const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
#endif
