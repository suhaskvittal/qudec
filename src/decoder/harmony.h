/* Contributor: OpenAI Codex. Harmony prototype, arXiv:2401.12434. */
#ifndef QUDEC_HARMONY_H
#define QUDEC_HARMONY_H

#include "decoder/types.h"
#include <stim.h>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <memory>

namespace dec {

// Harmonized correlated matching, Shutty, Newman, and Villalonga,
// arXiv:2401.12434. The DEM should have decomposed CSS error components.
// Each instance owns mutable matching workspaces; use one instance per worker.
// Pooling uses the paper's sum-likelihood method by default, with optional
// vote pooling and unanimous-sentinel escalation to the full ensemble.
class Harmony {
public:
    enum class Pooling { sum_likelihood, vote };
    struct Config {
        size_t ensemble_size = 100;
        size_t sentinel_size = 4;
        bool layered = true;
        uint64_t seed = 0x4841524d4f4e59ULL;
        double alpha_first = 1.0;
        double alpha_second = 0.8;
        double alpha_conditional = 0.5;
        Pooling pooling = Pooling::sum_likelihood;
    };
    struct Statistics {
        uint64_t shots = 0;
        uint64_t expanded_shots = 0;
        uint64_t members_run = 0;
        uint64_t decode_total_ns = 0;
        uint64_t decode_min_ns = std::numeric_limits<uint64_t>::max();
        uint64_t decode_max_ns = 0;
    };
    const size_t num_detectors;
    const size_t num_observables;
    explicit Harmony(const stim::DetectorErrorModel&);
    Harmony(const stim::DetectorErrorModel&, Config);
    ~Harmony();
    Harmony(const Harmony&) = delete;
    Harmony& operator=(const Harmony&) = delete;
    // Construct every member before sampling, keeping initialization out of latency.
    void prepare_ensemble();
    Result decode(SyndromeRef, ObsRef);
    Statistics statistics() const;
    void print_stats(std::ostream&) const;
    void mpi_accumulate();
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace dec
#endif
