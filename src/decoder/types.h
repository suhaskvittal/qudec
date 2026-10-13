/*
 *  author: Suhas Vittal
 *  contributor: OpenAI Codex
 * */

#ifndef DECODER_TYPES_h
#define DECODER_TYPES_h

#include <stim.h>

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <vector>

namespace dec
{

using SyndromeType = stim::simd_bits<64>;
using SyndromeRef = stim::simd_bits_range_ref<64>;
using ObsType = SyndromeType;
using ObsRef = SyndromeRef;

/*
 * `MatchingData` contains information that is useful for debugging
 * matching decoders.
 */
struct MatchingData
{
    using det_id_type = int64_t;
    using weight_type = uint64_t;
    struct assignment_type
    {
        det_id_type d1;
        det_id_type d2;
        double      pr;
        weight_type w_qu;
        ObsType     frame_flips;

        int matching_step{0};
        int cluster_id{0};
    };

    weight_type                  total_weight{0};
    double                       probability{1.0};
    ObsType                      frame_flips{1};
    std::vector<assignment_type> assignments;

    double gap;

    void add(const assignment_type&);
    void merge(const MatchingData&);
};

/*
 * `Result` is a generic result output that must be returned by a decoder.
 * Feel free to add extra fields if they are useful.
 */
struct Result
{
    ObsType flipped_obs{1};

    // Decoder specific parameters:
    MatchingData matching_data;
};

void matching_show_diff(std::ostream& ostrm, MatchingData, MatchingData, size_t num_observables);

} // namespace dec

#endif
