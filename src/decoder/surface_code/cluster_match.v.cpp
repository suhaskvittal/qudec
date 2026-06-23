/*
 *  author: Suhas Vittal
 *  date:   19 June 2026
 * */

#include "decoder/surface_code/cluster_match.h"
#include "verilator_utility.h"

#include "Vastrea.h"

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type
CLUSTER_MATCH::v_solve_matching_problem(matching_problem_type mp, Vastrea& astrea)
{
    constexpr size_t MAX_SUPPORTED_HW{6};
    constexpr size_t EDGE_COUNT = _get_mwpm_edge_count(MAX_SUPPORTED_HW);
    constexpr size_t WEIGHT_WIDTH{24};
    constexpr size_t EDGE_WIDTH{WEIGHT_WIDTH+1};

    if (num_observables > 1)
        std::cerr << "CLUSTER_MATCH::v_solve_matching_problem: Astrea is configured for 1-bit frames" << _die{};

    if (astrea_weight_quantization != quantization_level::b16)
        std::cerr << "CLUSTER_MATCH::v_solve_matching_problem: Astrea is configured for 8-bit quantization" << _die{};

    if (mp.detectors.size() > MAX_SUPPORTED_HW)
    {
        std::cerr << "CLUSTER_MATCH::v_solve_matching_problem: unsupported hamming weight: " 
                    << mp.detectors.size() << ", max allowed = " << MAX_SUPPORTED_HW << _die{};
    }

    // initialize `astrea` pins:
    astrea.edges = verilator_build_packed_array<EDGE_COUNT, EDGE_WIDTH>(mp.edges.begin(), mp.edges.end(),
                                    [] (const mwpm_edge_type& e)
                                    {
                                        constexpr uint64_t w_mask = (1ull << WEIGHT_WIDTH)-1;
                                        uint64_t x{0};
                                        x |= (e.w_qu & w_mask);  // bottom 8 bits
                                        if (e.frame_flips[0])  // 9th bit is frame
                                            x |= (1ull << WEIGHT_WIDTH);
                                        return x;
                                    });
    uint64_t edges_v{0};
    for (size_t i = 0; i < mp.edges.size(); i++)
        edges_v |= (1ull << i);
    astrea.edges_v = edges_v; 

    // run astrea:
    astrea.eval();

    // return output:
    syndrome_type frame_flips(num_observables);
    const size_t ff_words = (num_observables+7)/8;
    memcpy(frame_flips.u8, &astrea.frame_flips, ff_words);
    result_type out{ .flipped_obs=std::move(frame_flips), .matching_weight=astrea.m_weight };

    // if we need to validate, then also run `solve_matching_problem()` and check that the
    // results match
    auto sw_result = solve_matching_problem(mp);
    for (size_t i = 0; i < num_observables; i++)
    {
        if (out.flipped_obs[i] != sw_result.flipped_obs[i])
        {
            std::cerr << "CLUSTER_MATCH::v_solve_matching_problem: Astrea had mismatch with software" 
                        << "\n\thamming weight = " << mp.detectors.size() 
                        << "\n\tedge count = " << mp.edges.size()
                        << "\n\tHW matching weight = " << out.matching_weight
                        << "\n\tSW matching weight = " << sw_result.matching_weight;
            std::cerr << "\n\tedges:";
            for (const auto& e : mp.edges)
                std::cerr << "\n\t\t" << e.d1 << ", " << e.d2 << ", w = " << e.w_qu << ", f = " << e.frame_flips[0];
            std::cerr << _die{};
        }
    }

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder
