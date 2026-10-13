#ifndef QUDEC_DECODER_LIBRA_H
#define QUDEC_DECODER_LIBRA_H

#include "decoder/common.h"
#include <memory>
#include <cstdint>
#include <iosfwd>
#include <vector>

namespace dec {

// Matching synthesis, C. Jones, arXiv:2408.12135 (Sections III--V).
// A single logical observable and graphlike components in each CSS basis
// are required. Instances own mutable matching workspaces: one per worker.
class Libra {
public:
    struct Config {
        size_t ensemble_size = 60;
        double gap_threshold_db = 20;
        size_t positive_cycle_capacity = 30;
        uint64_t seed = 0x4c49425241ULL;
        bool degeneracy = true;
        bool validate = false;
        bool profile = false;
    };
    struct DetailedResult {
        bool prediction = false;
        bool ensemble_used = false;
        bool fallback_used = false;
        double gap_db = 0;
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

    explicit Libra(const stim::DetectorErrorModel&);
    Libra(const stim::DetectorErrorModel&, Config);
    ~Libra();
    Libra(Libra&&) noexcept;
    Libra(const Libra&) = delete;
    Libra& operator=(const Libra&) = delete;

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
