/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#ifndef DECODER_SURFACE_CODE_CLUSTER_MATCH_h
#define DECODER_SURFACE_CODE_CLUSTER_MATCH_h

#include "decoder/common.h"
#include "stats.h"

#include <stim.h>

#include <iosfwd>
#include <vector>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Verilator class forward declarations
 * */
class Vinitialize_neighbors;
class Vfilter;
class Vastrea;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


namespace decoder
{

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
     * detectors in the syndrome to run Astrea's matching step.
     * */
    using mwpm_edge_type = MATCHING_DATA::assignment_type;

    /*
     * This is the barebones information needed for matching: a list
     * of detectors and the edges between them.
     * */
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

    /*
     * 
     * */
    enum hw_emu_flag : uint8_t 
    { 
        filter = 0x1, 
        uf_cluster = 0x2, 
        synthesis = 0x4, 
        astrea = 0x8,
        validate = 0x10
    };

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

    /*
     * Hardware emulation setting (bitvector)
     * */
    const uint8_t hw_emu_enable;

    /*
     * Statistics:
     * */
    STATS_HISTOGRAM<uint64_t> s_clusters{0,32,4},
                                s_filtered{0,128,16},
                                s_hamming_weight{0,512,16},
                                s_post_filter_hamming_weight{0,512,16},
                                s_cluster_hamming_weight{0,12,2},
                                s_cluster_size{0,128,16},
                                s_growth_ticks{0,256,16},
                                s_synthesis_ticks{0,256,16},
                                s_synthesis_ticks_norm{0,256,16};

    /*
     * Hardware emulation statistics:
     * */
    STATS_HISTOGRAM<uint64_t> s_hw_filter_latency{0, 256, 16},
                                s_astrea_latency{0, 256, 16};
private:
    /*
     * Decoding graph adjacency data:
     * */
    std::vector<adj_list_type> adj_matrix_;
    adj_list_type boundary_adjacency_;
public:
    CLUSTER_MATCH(const stim::DetectorErrorModel&,
                    size_t code_distance,
                    size_t astrea_hw_max,
                    quantization_level,
                    uint8_t hw_emu_enable = 0);

    const adj_list_type& adj_matrix(det_id_type) const;

    result_type decode(syndrome_ref);

    void print_stats(std::ostream&) const;
private:
    /*
     * `filter_isolated_errors()` removes any isolated weight-1 errors from the syndrome.
     * */
    result_type filter_isolated_errors(syndrome_ref);

    /*
     * `uf_compute_clusters()` computes sub-clusters within the syndrome that correspond
     * to different matching problems. The size of a cluster is limited by `astrea_hw_max`,
     * and its width is limited by `code_distance`.
     * */
    std::vector<cluster_type> uf_compute_clusters(syndrome_ref);

    /*
     * `synthesize_matching_problem()` computes the pairwise distances for all detection
     * events within a cluster.
     * */
    matching_problem_type synthesize_matching_problem(cluster_type);

    /*
     * `solve_matching_problem()` computes the min-weight error for the given matching problem.
     * */
    result_type solve_matching_problem(matching_problem_type, int cluster_id);

    /*
     * Verilator emulation of the above functions.
     * */
    result_type v_filter_isolated_errors(syndrome_ref, Vinitialize_neighbors&, Vfilter&);
    result_type v_solve_matching_problem(matching_problem_type, Vastrea&);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder

#endif
