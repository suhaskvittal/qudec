/*
 *  author: Suhas Vittal
 *  contributor: OpenAI Codex
 *  date:   13 October 2026
 * */

#ifndef PP_PDZ_h
#define PP_PDZ_h

#include "decoder/types.h"
#include "preprocessor/affinity.h"

#include <limits>
#include <unordered_map>
#include <optional>

namespace pp
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


/*
 * `PDZ` = Predecoderization. This is essentially the concept of
 * taking the outputs of an initial decode and using them to seed
 * decisions of another decoder.
 *
 * For surface codes, this can be done by taking the matchings produced
 * by a decoder and eliminating the "obvious" matchings to reduce the 
 * amount of work done by another decoder.
 * */
class PDZ
{
public:
    using value_type = AffinityTable::value_type;
    using det_id_type = dec::MatchingData::det_id_type;

    struct config_type
    {
        int64_t aff_keep_threshold{2};
        bool    ignore_opposing_basis_detectors{true};
        bool    check_affinity_of_compl_match{false};
    };

    /*
     * `aff_match` is the affinity of the matched detectors.
     * `aff_compl_*` is the minimum quantized distance of either of the detectors to the rest of the
     *              flipped detectors. We track affinity for same basis and opposing
     *              basis detectors separately.
     * `aff_opposite_basis_match` is the match affinity of the opposing basis detector considered
     *              in `aff_compl_opposite_basis`
     *
     * */
    struct score_type
    {
        value_type aff_match{},
                   aff_compl_same_basis =std::numeric_limits<value_type>::max(),
                   aff_compl_opp_basis =std::numeric_limits<value_type>::max();
    
        std::optional<det_id_type> compl_same_basis_det_id,
                                   compl_opp_basis_det_id;

        value_type aff_same_basis_match{},
                   aff_opp_basis_match{};

        void update_affinity(det_id_type, value_type aff_1, value_type aff_2, value_type aff_m, bool same_basis);
    };

    const config_type conf;
        
    const stim::DetectorErrorModel& dem;
    const size_t num_detectors,
                 num_observables;

    /*
     * Statistics:
     * */
private:
    AffinityTable aff_;
    std::vector<bool> opposing_basis_;

    mutable std::vector<det_id_type> same_basis_events_, 
                                     opp_basis_events_;
    /*
     * Maps detectors to the idx of the assignment in `prev_dec_result` used below in 
     * `apply()` and `score_match()`. Should be initialized in `apply()`
     * */
    mutable std::unordered_map<det_id_type, size_t> assignment_map_;

    /*
     * Stores affinities for the matched pairs
     * */
    mutable std::unordered_map<det_id_type, value_type> aff_match_map_;
public:
    PDZ(const stim::DetectorErrorModel&);
    PDZ(const stim::DetectorErrorModel&, config_type);
    PDZ(const stim::DetectorErrorModel&, gzFile&);
    PDZ(const stim::DetectorErrorModel&, gzFile&, config_type);
    PDZ(const stim::DetectorErrorModel&, AffinityTable);
    PDZ(const stim::DetectorErrorModel&, AffinityTable, config_type);

    /*
     * `apply()` removes accepted matching endpoints from the input syndrome and
     * returns their matching data and observable correction. Combine this correction
     * with the decoder result for the remaining syndrome. One instance is not thread safe.
     * */
    dec::Result apply(dec::SyndromeRef, const dec::Result& prev_dec_result) const;
    // Read-only basis classification used by pipeline statistics. "Same basis"
    // means the memory basis; an opposing-basis pair is not counted as same basis.
    bool is_opposing_basis_detector(det_id_type) const;
private:
    score_type score_match(det_id_type, det_id_type) const;

    // wrapper for `aff_` that translates signed IDs to unsigned IDs.
    value_type affinity(det_id_type, det_id_type) const;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace pp

#endif

