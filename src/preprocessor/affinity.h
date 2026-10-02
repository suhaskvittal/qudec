/*
 *  author: Suhas Vittal
 *  date:   24 September 2026
 * */

#ifndef PP_AFFINITY_h
#define PP_AFFINITY_h

#include "hypergraph.h"

#include <zlib.h>

#include <vector>

namespace pp
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * `Affinity` measures what we call *affinity*, or the probability that two
 * detectors will ever be activated by the same error chain. Affinity is
 * pairwise.
 *
 * Affinity calculations require using a Hypergraph to represent a DEM
 * efficiently. We operate on this graph when computing the pairwise affinity
 * of a syndrome. THis also allows for more flexilibility across codes: for
 * example, a color code has three types of boundaries, so the user can provide
 * the specific boundaries to consider in the syndrome.
 * */

using AffinityResult = std::vector<double>;

/*
 * Reports the pairwise affinity for all detectors in the input syndrome. The input
 * syndrome should contain nodes in the input hypergraph. This should include boundary
 * nodes if needed.
 *
 * A minimum requirement is that `G::e_data_type` has a member `pr` that corresponds
 * to the error probability for the hyperedge. `G::v_data_type` also needs to have
 * a member `is_boundary`
 * */
template <class G>
AffinityResult measure_affinity(const G&, const std::vector<hg::id_type>& syndrome);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * In practice, we want to maintain a lookup table for pairwise-affinity. This table
 * is quantized.
 * */
class AffinityTable
{
public:
    using value_type = int8_t;
    constexpr static size_t BIT_WIDTH{2};
    constexpr static value_type MASK = static_cast<value_type>( (1<<BIT_WIDTH)-1 );

private: 
    size_t num_detectors_;
    std::vector<value_type> data_;
public:
    AffinityTable(const stim::DetectorErrorModel&);
    AffinityTable(gzFile&);

    /*
     * Retrieves the appropriate affinity value from `data_`.
     * */
    value_type operator()(size_t d1, size_t d2) const;

    /*
     * File IO: we use GZIP to keep table sizes small on disk. Latency doesn't matter once
     * the table is in memory.
     * */
    void write_to_file(gzFile&);

    size_t num_detectors() const { return num_detectors_; }
private:
    /*
     * We also provide `set()` which should only be used when constructing
     * the table or reading it from storage.
     * */
    void set(size_t d1, size_t d2, value_type);

    /*
     * Index of entry: returns idx of word in `data_` and bit offset
     * */
    std::tuple<size_t, size_t> index(size_t, size_t) const;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace aff
{

/*
 * `measure_affinity()` implementations
 * */
template <class G>
AffinityResult impl_dijkstra(const G&, const std::vector<hg::id_type>& syndrome);

/*
 * Hypergraph types.
 *
 * General function that converts DEMs into a hypergraph. We assume a max order of 16 in
 * our construction. 
 * */
struct v_data_type { bool is_boundary{false}; }
struct e_data_type { double pr; }
constexpr size_t HG_ORDER{16};
using hg_type = Hypergraph<v_data_type, e_data_type, HG_ORDER>;

hg_type build_hypergraph_from_dem(const stim::DetectorErrorModel&);

} // namespace aff

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace pp

#include "affinity.tpp"

#endif
