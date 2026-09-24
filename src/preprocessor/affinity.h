/*
 *  author: Suhas Vittal
 *  date:   24 September 2026
 * */

#ifndef PP_AFFINITY_h
#define PP_AFFINITY_h

#include "hypergraph.h"

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
 * to the error probability for the hyperedge.
 * */
template <class G>
AffinityResult measure_affinity(const G&, const std::vector<hg::id_type>& syndrome);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * `measure_affinity()` implementations
 * */

namespace aff
{

template <class G>
// OpenAI GPT-6: Accept the const syndrome supplied by measure_affinity().
AffinityResult impl_dijkstra(const G&, const std::vector<hg::id_type>& syndrome);

} // namespace aff

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace pp

#include "affinity.tpp"

#endif
