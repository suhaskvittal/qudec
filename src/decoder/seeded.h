/* Contributor: OpenAI Codex. Seeded correlation and optional pair retention. */
#ifndef QUDEC_DECODER_SEEDED_H
#define QUDEC_DECODER_SEEDED_H

#include "decoder/libra.h"
#include "decoder/pymatching.h"
#include "decoder/uf.h"

#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <vector>

namespace dec
{

// Correlated PyMatching with the first pass supplied by another decoder.
// prepare() and restore() bracket one shot; decode_prepared() may run on either
// the full syndrome or a residual syndrome while its graph is conditioned.
class CorrelatedPyMatching
{
public:
    const size_t num_detectors;
    const size_t num_observables;

    explicit CorrelatedPyMatching(const stim::DetectorErrorModel&);
    ~CorrelatedPyMatching();
    CorrelatedPyMatching(const CorrelatedPyMatching&) = delete;
    CorrelatedPyMatching& operator=(const CorrelatedPyMatching&) = delete;

    void prepare(const std::vector<int64_t>& selected_edges);
    Result decode_prepared(SyndromeRef, ObsRef);
    void restore() noexcept;
    pm::Mwpm& matcher() { return matcher_; }

    void print_stats(std::ostream&) const;
    void mpi_accumulate();

private:
    pm::Mwpm matcher_;
    bool prepared_{false};
    uint64_t calls_{0};
    uint64_t total_ns_{0};
};

// Reuses the unperturbed shortest route for each competing connection. Its
// scored cost uses the current conditioned graph. A pair is frozen only when
// its selected physical route is unchanged and no individual connection is
// strictly cheaper. Complete pairs are always kept or removed together.
class SeededPostprocessor
{
public:
    struct Output
    {
        SyndromeType remaining;
        ObsType fixed_frame;
        size_t frozen_pairs{0};
    };

    SeededPostprocessor(pm::Mwpm&, size_t detectors, size_t observables);
    ~SeededPostprocessor();
    SeededPostprocessor(const SeededPostprocessor&) = delete;
    SeededPostprocessor& operator=(const SeededPostprocessor&) = delete;

    Output apply(SyndromeRef, const MatchingData&);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

enum class SeedGapMode { Cluster, Complementary };

struct SeededMatchingConfig
{
    bool postprocess{false};
    SeedGapMode libra_gap_mode{SeedGapMode::Cluster};
    double libra_gap_threshold_db{40.0};
};

// D1 is PyMatching or UnionFind. D2 is CorrelatedPyMatching or Libra. Decoder
// instances are supplied by the caller and must outlive this pipeline.
template <typename D1, typename D2>
class SeededMatching
{
public:
    const size_t num_detectors;
    const size_t num_observables;

    SeededMatching(const stim::DetectorErrorModel&, D1&, D2&, SeededMatchingConfig = {});
    ~SeededMatching();

    Result decode(SyndromeRef, ObsRef);
    void print_stats(std::ostream&) const;
    void mpi_accumulate();

private:
    D1& seed_;
    D2& second_;
    SeededMatchingConfig config_;
    std::unique_ptr<CorrelatedPyMatching> libra_conditioning_;
    std::unique_ptr<SeededPostprocessor> postprocessor_;
    std::unique_ptr<PyMatching> complementary_gap_;
    uint64_t shots_{0};
    uint64_t failures_{0};
    uint64_t second_calls_{0};
    uint64_t gap_bypasses_{0};
    uint64_t frozen_pairs_{0};
    uint64_t input_events_{0};
    uint64_t residual_events_{0};
    uint64_t pipeline_total_ns_{0};
    uint64_t second_total_ns_{0};
};

} // namespace dec

#include "decoder/seeded.tpp"
#endif
