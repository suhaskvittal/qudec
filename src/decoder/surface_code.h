/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#ifndef DECODER_SURFACE_CODE_h
#define DECODER_SURFACE_CODE_h

#include "decoder/common.h"
#include "decoder/logger.h"

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>

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

    result_type decode(syndrome_ref, LOGGER&);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class CLUSTER_MATCH
{
public:
    using det_id_type = int64_t;

    /*
     * Need to represent relevant graph information
     * when decoding:
     *  `d`: the other detector
     *  `pr`: probability
     *  `frame_flips`: Pauli frame flips
     * */
    struct adj_entry_type
    {
        det_id_type d; 
        double      pr;
        obs_type    frame_flips;
    };

    using adj_list_type = std::vector<adj_entry_type>;

    /*
     * `cluster_type` is a group of detectors that corresponds to a matching
     * problem. `all` contains all detectors in the cluster, whereas
     * `flipped` contains only the detectors flipped in the syndrome.
     * */
    struct cluster_type
    {
        std::vector<det_id_type> all;
        std::vector<det_id_type> flipped;
    };

    /*
     * We need to compute pairwise distances between all
     * detectors in the syndrome to run Astrea's matching
     * step. This is the data we need:
     *  (1) both detector ids,
     *  (2) a integer quantized weight (see `quantization_level` below)
     *  (3) a double that is the unquantized weight
     *  (4) the Pauli frame flip
     * */
    struct mwpm_edge_type
    {
        det_id_type d1;
        det_id_type d2;
        uint64_t  w_qu;
        obs_type frame_flips;
    };

    struct matching_problem_type
    {
        std::vector<det_id_type> detectors;
        std::vector<mwpm_edge_type> edges;
    };
    
    /*
     * Quantization level:
     *   b4...b16 is an integer with the given number of bits
     *   fp just uses a full 64-bit float (`double`)
     * */
    enum class quantization_level { b4, b8, b16, b32 };

    const size_t num_detectors;
    const size_t num_observables;

    /*
     * Code distance of surface code: determines maximum
     * amount of cluster growth (grow up-to d/2).
     * */
    const size_t code_distance;

    /*
     * Maximum size of a cluster is determined by `astrea_max_hw`
     * */
    const size_t astrea_hw_max;

    /*
     * Weight quantization for Astrea.
     * */
    const quantization_level astrea_weight_quantization;
private:
    /*
     * Decoding graph adjacency data:
     * */
    std::vector<adj_list_type> adj_matrix_;
public:
    CLUSTER_MATCH(const stim::DetectorErrorModel&, 
                    size_t code_distance, 
                    size_t astrea_hw_max,
                    quantization_level);

    result_type decode(syndrome_ref, LOGGER&);

    /*
     * These are some useful values that can inform RTL implementation.
     * */
    size_t max_degree() const;
private:
    std::vector<cluster_type> uf_compute_clusters(syndrome_ref, LOGGER&);
    matching_problem_type synthesize_matching_problem(cluster_type&&, LOGGER&);
    result_type solve_matching_problem(matching_problem_type&&, LOGGER&);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder

#endif // DECODER_SURFACE_CODE_h
