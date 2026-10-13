/*
 *  contributor: OpenAI Codex
 * */

#ifndef DECODER_UF_h
#define DECODER_UF_h

#include "decoder/types.h"

#include <cstdint>
#include <iosfwd>
#include <limits>
#include <vector>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Weight-aware union-find decoding of graphlike surface-code DEMs.
 * Separated components are independent; no correlation reweighting is used.
 * Odd clusters grow until they become even or touch a boundary. Growth fills
 * an edge from both ends, with length proportional to |log((1-p)/p)|.
 * A spanning forest of the grown edges is then peeled to recover a correction.
 *
 * Weights are rounded to multiples of 1 / weight_scale (natural log units).
 * Unlike MWPM, union-find does not minimize the final correction weight.
 * This implementation scans edges at each growth event; it does not claim the
 * almost-linear complexity of specialized frontier-based implementations.
 * One instance per worker: decode() reuses mutable scratch storage.
 */
class UnionFind
{
public:
    struct config_type
    {
        uint64_t weight_scale{1024};
    };

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

    explicit UnionFind(const stim::DetectorErrorModel&);
    UnionFind(const stim::DetectorErrorModel&, config_type);

    Result decode(SyndromeRef, ObsRef);

    Statistics statistics() const;
    void print_stats(std::ostream&) const;
    void mpi_accumulate();

private:
    struct edge_type
    {
        size_t u;
        size_t v;
        uint64_t length;
        std::vector<size_t> observables;
    };

    std::vector<edge_type> edges_;
    size_t num_vertices_;
    SyndromeType baseline_syndrome_;
    ObsType baseline_observables_;
    Statistics stats_;

    // Disjoint sets hold parity and a boundary flag at their roots.
    std::vector<size_t> parent_;
    std::vector<size_t> cluster_size_;
    std::vector<uint8_t> parity_;
    std::vector<uint8_t> boundary_;
    std::vector<uint64_t> growth_;
    std::vector<size_t> forest_edges_;

    size_t find_root(size_t);
    void merge_clusters(size_t edge_id);
    void grow_clusters();
    void peel_forest(std::vector<uint8_t>& syndrome, ObsRef correction) const;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec

#endif // DECODER_UF_h
