/*
 *  author: Suhas Vittal
 *  date:   22 June 2026
 * */

#include "decoder/surface_code/matching_data.h"
#include "globals.h"

#include <iostream>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

std::string _frames_to_string(ObsRef, size_t num_observables);

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
MatchingData::add(const assignment_type& a)
{
    total_weight += a.w_qu;
    probability *= a.pr;
    if (assignments.empty())
        frame_flips = a.frame_flips;
    else
        frame_flips ^= a.frame_flips;
    assignments.push_back(a);
}

void
MatchingData::merge(const MatchingData& other)
{
    total_weight += other.total_weight;
    probability *= other.probability;
    if (assignments.empty())
        frame_flips = other.frame_flips;
    else
        frame_flips ^= other.frame_flips;
    for (const auto& a : other.assignments)
        assignments.push_back(a);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
matching_show_diff(std::ostream& ostrm, MatchingData m1, MatchingData m2, size_t num_observables)
{
    auto ff_diff = m1.frame_flips ^ m2.frame_flips;
    bool any_mismatch{false};
    for (size_t i = 0; i < num_observables; i++)
        any_mismatch |= (ff_diff[i] > 0);
    if (!any_mismatch)
        return;

    ostrm << "DIFF--------------------------------------\n";
    ostrm << "frame flips: ";
    for (size_t i = 0; i < num_observables; i++)
    {
        if (ff_diff[i])
            ostrm << "!";
        else
            ostrm << ".";
    }
    ostrm << "\ntotal weight: M1 = " << m1.total_weight << ", M2 = " << m2.total_weight
            << "\nprobability: M1 = " << m1.probability << ", M2 = " << m2.probability;
    ostrm << "\nmatching:";
    for (const auto& a : m1.assignments)
    {
        auto ad1 = a.d1,
             ad2 = a.d2;
        if (ad1 > ad2)
            std::swap(ad1, ad2);

        ostrm << "\n\tM!: " << ad1 << " <---> " << ad2 
                << ", W = " << a.w_qu 
                << ", PR = " << a.pr 
                << ", step = " << a.matching_step
                << ", cid = " << a.cluster_id
                << ", f = " << _frames_to_string(a.frame_flips, num_observables)
                << "\tM2:";
        for (auto ad : {ad1,ad2})
        {
            // find corresponding matching in a2.
            auto b_it = std::find_if(m2.assignments.begin(), m2.assignments.end(),
                                    [ad] (const auto& b) { return b.d1 == ad || b.d2 == ad; });
            if (b_it == m2.assignments.end())
            {
                ostrm << "\tN/A";
                continue;
            }
            const auto& b = *b_it;
            auto bd1 = b.d1,
                 bd2 = b.d2;
            if (bd1 > bd2)
                std::swap(bd1, bd2);
            
            bool frame_mismatch{false};
            for (size_t i = 0; i < num_observables; i++)
                if (a.frame_flips[i] != b.frame_flips[i])
                    frame_mismatch = true;

            if (ad1 != bd1 || ad2 != bd2 || frame_mismatch)
            {
                ostrm << "\t" << bd1 << " <---> " << bd2 
                        << ", W = " << b.w_qu 
                        << ", PR = " << b.pr
                        << ", step = " << b.matching_step
                        << ", cid = " << b.cluster_id
                        << ", f = " << _frames_to_string(b.frame_flips, num_observables);
            }
        }
    }
    ostrm << "\n";
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

std::string
_frames_to_string(ObsRef frame_flips, size_t num_observables)
{
    std::stringstream ss;
    for (size_t i = 0; i < num_observables; i++)
    {
        if (frame_flips[i])
            ss << "1";
        else
            ss << "0";
    }
    return ss.str();
}

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec
