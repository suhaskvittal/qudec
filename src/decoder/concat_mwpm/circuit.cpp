/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Full spacetime concatenated matching with local flag-readout likelihoods.
 *  A flag is an independent per-round detector, never a round difference.
 * */
#include "decoder/concat_mwpm/circuit.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace dec::concat
{
namespace
{
std::vector<CircuitFault>
read_faults(const stim::DetectorErrorModel& dem)
{
    std::vector<CircuitFault> faults;
    dem.iter_flatten_error_instructions([&](const auto& instruction)
    {
        CircuitFault fault{instruction.arg_data[0], {}, 0};
        for (auto target : instruction.target_data)
        {
            if (target.is_relative_detector_id()) fault.detectors.push_back(target.val());
            else if (target.is_observable_id()) fault.observable ^= uint64_t{1} << target.val();
            else throw std::invalid_argument("Need original, undecomposed DEM events");
        }
        faults.push_back(std::move(fault));
    });
    return faults;
}

std::map<uint64_t, std::vector<double>>
read_coordinates(const stim::DetectorErrorModel& dem)
{
    std::set<uint64_t> ids;
    for (size_t id = 0; id < dem.count_detectors(); id++) ids.insert(id);
    return dem.get_detector_coordinates(ids);
}
} // namespace

CircuitDecoder::CircuitDecoder(const stim::DetectorErrorModel& dem, bool memory_x,
                               bool use_flags, size_t cache_capacity)
    :num_detectors(dem.count_detectors()), faults(read_faults(dem)),
     coordinates(read_coordinates(dem)),
     measurement_p_(num_detectors, 0), use_flags_(use_flags),
     cache_capacity_(cache_capacity)
{
    if (dem.count_observables() != 1 || cache_capacity == 0)
        throw std::invalid_argument("Concatenated MWPM needs one observable and positive cache capacity");
    std::vector<size_t> local(num_detectors, SIZE_MAX);
    for (const auto& [global, c] : coordinates)
    {
        if (c.size() != 6) throw std::invalid_argument("Missing color/basis/flag detector metadata");
        if (c[5] == 1) flag_globals.push_back(global);
        else if (c[5] == 0)
        {
            if ((c[3] != 0 && c[3] != 1 && c[3] != 2) || (c[4] != 0 && c[4] != 1))
                throw std::invalid_argument("Invalid check color or basis metadata");
            if (c[4] == (memory_x ? 0 : 1))
            {
                local[global] = colors_.size();
                colors_.push_back(c[3]);
                check_globals.push_back(global);
            }
        }
        else throw std::invalid_argument("Invalid flag-detector metadata");
    }
    if (colors_.empty()) throw std::invalid_argument("No check detectors in the selected memory basis");
    // Author: OpenAI Codex, 22 October 2026. Infer each flag's readout-error
    // probability from isolated original DEM errors. Check BOTH CSS bases
    // before marginalization; a flag+opposite-basis-check fault is not an
    // isolated flag flip. Correlated multi-flag faults are not singletons.
    for (const auto& fault : faults)
        if (fault.detectors.size() == 1 && fault.observable == 0 &&
            coordinates.at(fault.detectors[0])[5] == 1)
        {
            auto& old = measurement_p_[fault.detectors[0]];
            old = old+fault.probability-2*old*fault.probability;
        }
    for (auto f : flag_globals)
    {
        if (!(measurement_p_[f] >= 0 && measurement_p_[f] < .5))
            throw std::invalid_argument("Invalid isolated-flag DEM probability");
    }
    // Marginalize opposite-basis CHECKS, retain flags of BOTH bases. Merge
    // independent variants only when checks, flags and logical effects agree.
    using Key = std::tuple<std::vector<size_t>, uint64_t, std::vector<size_t>>;
    std::map<Key, double> probabilities;
    for (const auto& fault : faults)
    {
        std::vector<size_t> checks, flags;
        for (auto v : fault.detectors)
            if (coordinates.at(v)[5] == 1) flags.push_back(v);
            else if (local[v] != SIZE_MAX) checks.push_back(local[v]);
        if (checks.empty())
        {
            if (fault.observable) throw std::invalid_argument("Undetectable single logical fault");
            continue;
        }
        std::sort(checks.begin(), checks.end());
        std::sort(flags.begin(), flags.end());
        auto& old = probabilities[{checks, fault.observable, flags}];
        old = old+fault.probability-2*old*fault.probability;
    }
    std::map<std::pair<std::vector<size_t>, uint64_t>, size_t> class_ids;
    for (const auto& [key, probability] : probabilities)
    {
        const auto& [checks, observable, flags] = key;
        auto [it, inserted] = class_ids.emplace(std::pair{checks, observable}, classes_.size());
        if (inserted) classes_.push_back({checks, observable, {}, {}});
        classes_[it->second].variants.push_back({flags, probability});
    }
    for (auto& group : classes_)
    {
        std::set<size_t> flags;
        for (const auto& variant : group.variants) flags.insert(variant.flags.begin(), variant.flags.end());
        group.flags.assign(flags.begin(), flags.end());
        for (int c = 0; c < 3; c++)
        {
            size_t own = 0;
            for (auto v : group.checks) own += colors_[v] == c;
            size_t other = group.checks.size()-own;
            if (other > 2 || (other && own > 1) || (!other && own > 2))
                throw std::invalid_argument("Circuit event does not fit Algorithm 1 for every color");
        }
    }
}

CircuitDecoder::Context
CircuitDecoder::make_context(const std::vector<uint8_t>& syndrome) const
{
    Model model{colors_, {}};
    Context out;
    size_t total_active = 0;
    double all_flag_factor = 1;
    std::vector<size_t> observed_flags;
    for (auto f : flag_globals)
        if (syndrome[f])
        {
            total_active++;
            all_flag_factor *= measurement_p_[f];
            observed_flags.push_back(f);
        }
    for (const auto& group : classes_)
    {
        size_t active = 0;
        double local_flag_factor = 1;
        for (auto f : group.flags)
            if (syndrome[f])
            {
                active++;
                local_flag_factor *= measurement_p_[f];
            }
        double unflagged = 0, effective = 0;
        for (const auto& variant : group.variants)
        {
            if (variant.flags.empty()) unflagged = variant.probability;
            if (!use_flags_)
            {
                effective += variant.probability-2*effective*variant.probability;
                continue;
            }
            if (active == 0)
            {
                if (variant.flags.empty()) effective = variant.probability;
                continue;
            }
            size_t overlap = 0;
            for (auto f : variant.flags) overlap += syndrome[f] != 0;
            // A flagged channel requires at least one associated flag to fire.
            if (!variant.flags.empty() && overlap == 0) continue;
            double likelihood = 1;
            for (auto f : group.flags)
            {
                bool expected = std::binary_search(variant.flags.begin(),variant.flags.end(),f);
                likelihood *= expected == bool(syndrome[f]) ? 1-measurement_p_[f] : measurement_p_[f];
            }
            effective += likelihood*variant.probability;
        }
        // Author: OpenAI Codex, 22 October 2026. Equation 9's flag-mismatch
        // factor uses the FULL observed flag set F, not only class-local F.
        // The mixture above accounts for mismatches within group.flags.
        // Every observed flag outside that set supplies one additional pm,
        // including when the class has no associated flags at all.
        // This adds the missing flag penalty while preserving the existing
        // probability mixture and event arity convention for concatenation.
        if (use_flags_ && total_active > active)
        {
            double outside_factor = all_flag_factor;
            if (local_flag_factor > 0) outside_factor /= local_flag_factor;
            else
            {
                outside_factor = 1;
                for (auto f : observed_flags)
                    if (!std::binary_search(group.flags.begin(),group.flags.end(),f))
                        outside_factor *= measurement_p_[f];
            }
            effective *= outside_factor;
        }
        // For one relevant fired flag this is precisely (1-pm)Y + pm*X.
        // Multiple relevant flags use the local independent readout likelihood
        // without dividing by the flag's marginal probability. This is the
        // requested reweighting heuristic, not an exact global posterior.
        out.weights.push_back({group.checks, group.observable, unflagged, effective,
                               group.flags.size(), active});
        if (effective > 0) model.events.push_back({effective, group.checks, group.observable});
    }
    out.decoder = std::make_unique<Decoder>(std::move(model));
    return out;
}

CircuitDecoder::Context&
CircuitDecoder::context(const std::vector<uint8_t>& syndrome)
{
    if (syndrome.size() != num_detectors) throw std::invalid_argument("Wrong full syndrome length");
    std::string key;
    for (auto f : flag_globals) key.push_back(use_flags_ && syndrome[f]);
    auto it = contexts_.find(key);
    if (it == contexts_.end())
    {
        if (contexts_.size() >= cache_capacity_) contexts_.clear();
        it = contexts_.emplace(std::move(key), make_context(syndrome)).first;
    }
    return it->second;
}

std::vector<uint8_t>
CircuitDecoder::memory_syndrome(const std::vector<uint8_t>& full) const
{
    std::vector<uint8_t> out;
    for (auto global : check_globals) out.push_back(full[global]);
    return out;
}

uint64_t
CircuitDecoder::decode(const std::vector<uint8_t>& full)
{
    return context(full).decoder->decode(memory_syndrome(full)).observable;
}

CircuitResult
CircuitDecoder::inspect(const std::vector<uint8_t>& full)
{
    auto& ctx = context(full);
    auto colors = ctx.decoder->decode_colors(memory_syndrome(full));
    Result best = colors[0];
    for (size_t c = 1; c < colors.size(); c++)
        if (colors[c].weight < best.weight-1e-9) best = colors[c];
    return {std::move(best), std::move(colors), ctx.weights, ctx.decoder->model,
            ctx.decoder->audit_weights()};
}

std::vector<uint8_t>
CircuitDecoder::fault_syndrome(const std::vector<size_t>& ids) const
{
    std::vector<uint8_t> out(num_detectors);
    for (auto id : ids) for (auto v : faults.at(id).detectors) out[v] ^= 1;
    return out;
}

uint64_t
CircuitDecoder::fault_observable(const std::vector<size_t>& ids) const
{
    uint64_t out = 0;
    for (auto id : ids) out ^= faults.at(id).observable;
    return out;
}
} // namespace dec::concat
