/*
 *  author: Suhas Vittal
 *  contributor: OpenAI Codex
 *  date:   13 October 2026
 * */

#include "preprocessor/pdz.h"
#include "globals.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace pp
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PDZ::PDZ(const stim::DetectorErrorModel& model)
    :PDZ(model, config_type{})
{}

PDZ::PDZ(const stim::DetectorErrorModel& model, config_type config)
    :PDZ(model, AffinityTable(model), config)
{}

PDZ::PDZ(const stim::DetectorErrorModel& model, gzFile& file)
    :PDZ(model, file, config_type{})
{}

PDZ::PDZ(const stim::DetectorErrorModel& model, gzFile& file, config_type config)
    :PDZ(model, AffinityTable(file), config)
{}

PDZ::PDZ(const stim::DetectorErrorModel& model, AffinityTable table)
    :PDZ(model, std::move(table), config_type{})
{}

PDZ::PDZ(const stim::DetectorErrorModel& model, AffinityTable table, config_type config)
    :conf(config), dem(model), num_detectors(model.count_detectors()),
    num_observables(model.count_observables()), aff_(std::move(table)),
    opposing_basis_(num_detectors, false)
{
    if (aff_.num_detectors() != num_detectors)
        std::cerr << "PDZ: affinity table does not match DEM detector count" << _die{};
    if (!conf.ignore_opposing_basis_detectors)
        std::cerr << "PDZ: opposing-basis scoring is not implemented yet" << _die{};
    if (conf.aff_keep_threshold < 0)
        std::cerr << "PDZ: affinity threshold must be nonnegative" << _die{};

    // The SI1000 generators place checks at even integer coordinates. The first
    // detector belongs to the memory basis; checkerboard parity identifies the other.
    const auto coordinates = dem.get_detector_coordinates({});
    bool have_coordinates = num_detectors > 0;
    for (size_t d = 0; d < num_detectors; d++)
    {
        const auto it = coordinates.find(d);
        if (it == coordinates.end() || it->second.size() < 2)
        {
            have_coordinates = false;
            break;
        }
        for (size_t k = 0; k < 2; k++)
            if (!std::isfinite(it->second[k]) ||
                    std::abs(it->second[k] / 2 - std::round(it->second[k] / 2)) > 1e-6)
                have_coordinates = false;
    }
    if (have_coordinates)
    {
        auto parity = [&] (size_t d)
        {
            const auto& c = coordinates.at(d);
            return std::fmod(std::abs(std::round(c[0] / 2) + std::round(c[1] / 2)), 2.0);
        };
        for (size_t d = 0; d < num_detectors; d++)
            opposing_basis_[d] = parity(d) != parity(0);
        return;
    }

    // Without coordinates, infer CSS components from a graphlike/decomposed DEM.
    // Boundary nodes must not unite the two bases. Logical-bearing components
    // identify the memory basis; an observable-free DEM treats all checks alike.
    std::vector<size_t> parent(num_detectors);
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&] (size_t d)
    {
        while (parent[d] != d)
        {
            parent[d] = parent[parent[d]];
            d = parent[d];
        }
        return d;
    };
    std::vector<size_t> logical_detectors;
    dem.iter_flatten_error_instructions([&] (const auto& instruction)
    {
        instruction.for_separated_targets([&] (const auto& group)
        {
            std::vector<size_t> detectors;
            bool logical = false;
            for (const auto& target : group)
                if (target.is_relative_detector_id()) detectors.push_back(target.val());
                else if (target.is_observable_id()) logical = true;
            if (detectors.size() > 2)
                std::cerr << "PDZ: basis inference needs SI1000 coordinates or a decomposed DEM" << _die{};
            if (detectors.size() == 2)
                parent[root(detectors[1])] = root(detectors[0]);
            if (logical)
                logical_detectors.insert(logical_detectors.end(), detectors.begin(), detectors.end());
        });
    });
    if (!logical_detectors.empty())
    {
        std::vector<bool> logical_components(num_detectors, false);
        for (auto d : logical_detectors)
            logical_components[root(d)] = true;
        for (size_t d = 0; d < num_detectors; d++)
            opposing_basis_[d] = !logical_components[root(d)];
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
PDZ::score_type::update_affinity(det_id_type x, value_type aff_1, value_type aff_2, value_type aff_m, bool same_basis)
{
    value_type a = std::min(aff_1, aff_2);
    
    auto& aff_curr =    same_basis ? aff_compl_same_basis : aff_compl_opp_basis;
    auto& det_id =      same_basis ? compl_same_basis_det_id : compl_opp_basis_det_id;
    auto& aff_curr_m =  same_basis ? aff_same_basis_match : aff_opp_basis_match;

    if (a < aff_curr)
    {
        aff_curr = a;
        det_id = x;
        aff_curr_m = aff_m;
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

dec::Result
PDZ::apply(dec::SyndromeRef s, const dec::Result& prev) const
{
    if (s.num_bits_padded() < num_detectors)
        std::cerr << "PDZ::apply(): syndrome is too short" << _die{};

    // Clear at entry as well: validation can throw during the previous call.
    same_basis_events_.clear();
    opp_basis_events_.clear();
    assignment_map_.clear();
    aff_match_map_.clear();
    same_basis_events_.reserve(s.popcnt());
    opp_basis_events_.reserve(s.popcnt());
    for (size_t i = 0; i < num_detectors; i++)
    {
        if (s[i])
        {
            if (is_opposing_basis_detector(i)) opp_basis_events_.push_back(i);
            else                              same_basis_events_.push_back(i);
        }
    }

    // Validate before modifying the syndrome. The boundary can occur many times;
    // each real detector must occur exactly once in the supplied matching.
    assignment_map_.reserve(s.popcnt());
    aff_match_map_.reserve(s.popcnt());
    for (size_t i = 0; i < prev.matching_data.assignments.size(); i++)
    {
        const auto& a = prev.matching_data.assignments[i];
        const auto aff = affinity(a.d1, a.d2);
        for (auto d : {a.d1, a.d2})
            if (d >= 0)
                aff_match_map_.emplace(d, aff);
    }

    dec::Result out{};
    out.flipped_obs = dec::ObsType(num_observables);
    out.matching_data.frame_flips = dec::ObsType(num_observables);
    for (const auto& a : prev.matching_data.assignments)
    {
        const auto x = score_match(a.d1, a.d2);
        // Preserve the existing scoring rule. A negative margin rejects the
        // assignment; it is not an invalid decoder result.
        auto diff = x.aff_compl_same_basis - x.aff_match;
        bool ok = (diff >= conf.aff_keep_threshold);
        if (ok && diff == 1 && conf.check_affinity_of_compl_match)
        {
            bool retain = (x.aff_same_basis_match != x.aff_match) || (x.aff_match >= 2);
            ok = ok && retain;
        }

        if (ok)
        {
            out.matching_data.add(a);
            out.flipped_obs ^= a.frame_flips;
            for (auto d : {a.d1, a.d2})
                if (d >= 0)
                    s[d] = false;
        }
    }

    same_basis_events_.clear();
    opp_basis_events_.clear();
    assignment_map_.clear();
    aff_match_map_.clear();

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PDZ::score_type
PDZ::score_match(det_id_type a, det_id_type b) const
{
    PDZ::score_type out{};
    out.aff_match = affinity(a, b);

    if (a >= 0 && b >= 0)
        assert(aff_match_map_.at(a) == aff_match_map_.at(b));

    // Compare against other fired detectors in the basis of this pair.
    const auto detector = a >= 0 ? a : b;
    const bool opposing = is_opposing_basis_detector(detector);
    const auto& same = opposing ? opp_basis_events_ : same_basis_events_;
    const auto& opposite = opposing ? same_basis_events_ : opp_basis_events_;
    for (det_id_type x : same)
        if (x != a && x != b)
            out.update_affinity(x, affinity(a,x), affinity(b,x), aff_match_map_[x], true);

    if (!conf.ignore_opposing_basis_detectors)
        for (det_id_type x : opposite)
            out.update_affinity(x, affinity(a,x), affinity(b,x), aff_match_map_[x], false);

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

bool
PDZ::is_opposing_basis_detector(det_id_type d) const
{
    if (d == -1)
        return false;
    return opposing_basis_[d];
}

PDZ::value_type
PDZ::affinity(det_id_type a, det_id_type b) const
{
    auto translate = [&] (det_id_type d) -> size_t
    {
        if (d == -1)
            return num_detectors;
        return static_cast<size_t>(d);
    };
    return aff_(translate(a), translate(b));
}

} // namespace pp
