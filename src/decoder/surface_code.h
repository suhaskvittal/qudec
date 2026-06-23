/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#ifndef DECODER_SURFACE_CODE_h
#define DECODER_SURFACE_CODE_h

#include "decoder/common.h"

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>

#include <iosfwd>
#include <vector>

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class PYMATCHING
{
public:
    const size_t num_detectors;
    const size_t num_observables;
private:
    pm::Mwpm mwpm_;
public:
    PYMATCHING(const stim::DetectorErrorModel&);

    result_type decode(syndrome_ref);

    void print_stats(std::ostream&) const {}
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * `BLOSSOMV` is a global Minimum-Weight Perfect Matching (MWPM) decoder built on
 * top of Blossom-5. Unlike `CLUSTER_MATCH`, it does not split the syndrome into
 * clusters: it matches all detection events at once. As such, it is an exact
 * software reference / ground-truth decoder whose accuracy should track `PYMATCHING`.
 *
 * The decoder works in three steps:
 *   (1) The constructor encaches the decoding graph (adjacency matrix + boundary
 *       adjacency) from the detector error model.
 *   (2) On each syndrome, pairwise distances between all detection events are
 *       computed via Dijkstra over the full decoding graph.
 *   (3) Blossom-5 computes the MWPM and the Pauli frame flips of the matched edges
 *       are applied as the correction.
 * */
class BLOSSOMV
{
public:
    using det_id_type = int64_t;

    /*
     * Decoding graph edge: the other detector `d`, its probability `pr`, and the
     * Pauli frame flips `frame_flips` associated with the edge.
     * */
    struct adj_entry_type
    {
        det_id_type d;
        double      pr;
        obs_type    frame_flips;
    };

    using adj_list_type = std::vector<adj_entry_type>;

    /*
     * An edge in the matching problem: both endpoints, a quantized weight, and the
     * accumulated Pauli frame flips along the shortest path between the endpoints.
     * */
    struct mwpm_edge_type
    {
        det_id_type d1;
        det_id_type d2;
        uint64_t    w_qu;
        obs_type    frame_flips;
    };

    const size_t num_detectors;
    const size_t num_observables;

    BLOSSOMV(const stim::DetectorErrorModel&);

    result_type decode(syndrome_ref);

    void print_stats(std::ostream&) const {}
private:
    /*
     * Decoding graph adjacency data:
     * */
    std::vector<adj_list_type> adj_matrix_;
    adj_list_type boundary_adjacency_;

    const adj_list_type& adj_matrix(det_id_type) const;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder

#include "decoder/surface_code/cluster_match.h"

#endif // DECODER_SURFACE_CODE_h
