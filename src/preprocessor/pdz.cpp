/*
 *  author: Suhas Vittal
 *  date:   13 October 2026
 * */

#include "preprocessor/pdz.h"

#include <cassert>

namespace pp
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
PDZ::score_type::update_same_basis_affinity(det_id_type x, value_type aff_1, value_type aff2)
{
    value_type a = std::max(aff_1, aff_2);
    if (a > aff_compl_same_basis)
    {
        aff_compl_same_basis = a;
        compl_same_basis_det_id = x;
    }
}

void
PDZ::score_type::update_opp_basis_affinity(det_id_type x, value_type aff_1, value_type aff_2, value_type aff_opp_m)
{
    value_type a = std::max(aff_1, aff_2);
    if (a > aff_compl_opp_basis)
    {
        aff_compl_opp_basis = a;
        compl_opp_basis_det_id = x;
        aff_opp_basis_match = aff_opp_m;
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

dec::Result
PDZ::apply(dec::SyndromeRef s, const dec::Result& prev) const
{
    // synthesize all detection events from s:
    same_basis_events_.reserve(s.popcnt());
    opp_basis_events_.reserve(s.popcnt());
    for (det_id_type i = 0; i < num_detectors; i++)
    {
        if (s[i])
        {
            if (is_opposing_basis_detector(i)) same_basis_events_.push_back(i);
            else                               opp_basis_events_.push_back(i);
        }
    }

    // go through `prev.matching_data` and initialize `assignment_map_`
    assignment_map_.reserve(s.popcnt());
    for (size_t i = 0; i < prev.matching_data.assignments.size(); i++)
    {
        auto aff = affinity(a.d1, a.d2);
        for (auto d : {a.d1, a.d2})
        {
            assignment_map_.insert({d, i});
            aff_match_map_.insert({d, aff});
        }
    }

    // score each pair in `prev` and decide whether the pair should be retained
    dec::Result out{};
    for (auto& a : prev.matching_data.assignments)
    {
        auto d1 = a.d1, d2 = a.d2;
        score_type x = score_match(d1, d2, prev);

        value_type aff_match{x.aff_match},
                   aff_compl;
        if (conf.ignore_opposing_basis_detectors)
            aff_compl = x.aff_compl_same_basis;
        // TODO: other branch

        value_type diff = aff_compl - aff_match;
        if (diff < 0)
        {
            std::cerr << "PDZ::apply(): unexpected -- diff < 0: aff_compl = " << aff_compl
                        << ", aff_match = " << aff_match
                        << ", diff = " << diff
                        << ", detectors = " << d1 << ", " << d2
                        << _die{}:
        }

        if (diff >= conf.aff_keep_threshold)
        {
            out.add(a);
            for (auto d : {d1,d2})
                if (d >= 0)
                    s ^= d;
        }
    }

    // cleanup:
    x_events_.clear();
    z_events_.clear();
    assignment_map_.clear();

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PDZ::score_type
PDZ::score_match(det_id_type a, det_id_type b, const dec::Result& prev) const
{
    PDZ::score_type out{};
    out.aff_match = aff_match_map[a];

    assert(aff_match_map_[a] == aff_match_map_[b]);

    // first check `same_basis_events_`
    for (det_id_type x : same_basis_events_)
        out.update_same_basis_affinity(x, affinity(a,x), affinity(b,x));

    if (conf.ignore_opposing_basis_detectors)
        for (det_id_type x : opp_basis_events_)
            out.update_opp_basis_affinity(x, affinity(a,x), affinity(b,x), aff_match_map_[x]);

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace pp
