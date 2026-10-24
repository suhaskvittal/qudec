/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Full spacetime concatenated matching with flag-selected class representatives.
 *  A flag is an independent per-round detector, never a round difference.
 * */
#include "decoder/concat_mwpm/circuit.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <iomanip>

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
        // Contributor: OpenAI Codex. Stim preserves native generator tags in
        // the DEM, including after per-component noise rewriting.
        if (instruction.tag == "qudec_data") fault.boundary_hint = BoundaryHint::Spatial;
        else if (instruction.tag == "qudec_measurement") fault.boundary_hint = BoundaryHint::Temporal;
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
                               bool use_flags, size_t cache_capacity, bool boundary_gap, EnsembleConfig ensemble)
    :num_detectors(dem.count_detectors()), faults(read_faults(dem)),
     coordinates(read_coordinates(dem)),
     flag_p_(num_detectors,0), measurement_p_(num_detectors, 0), use_flags_(use_flags),
     cache_capacity_(cache_capacity), boundary_gap_(boundary_gap), ensemble_(ensemble)
{
    if (ensemble.size == 0 || !(ensemble.spread_scale >= 0) || !std::isfinite(ensemble.spread_scale))
        throw std::invalid_argument("Invalid concatenated matching ensemble configuration");
    if (dem.count_observables() != 1 || cache_capacity == 0)
        throw std::invalid_argument("Concatenated MWPM needs one observable and positive cache capacity");
    std::vector<size_t> local(num_detectors, SIZE_MAX);
    std::vector<bool> keep_flag(num_detectors,false);
    for (const auto& [global, c] : coordinates)
    {
        if (c.size() != 6) throw std::invalid_argument("Missing color/basis/flag detector metadata");
        if (c[5] == 1)
        {
            if (c[4] != 0 && c[4] != 1) throw std::invalid_argument("Invalid flag extraction basis metadata");
            // Contributor: OpenAI Codex. X-check flags detect X hooks relevant
            // to Z-memory checks, and conversely for Z-check flags in X memory.
            if (c[4] == (memory_x ? 1 : 0))
            {
                keep_flag[global] = true;
                flag_globals.push_back(global);
            }
        }
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
    // Contributor: OpenAI Codex. Every original channel that flips a retained
    // flag contributes to its parity marginal, including data-correlated faults.
    for (const auto& fault : faults)
        for (auto f : fault.detectors) if (keep_flag[f])
            flag_p_[f] += fault.probability-2*flag_p_[f]*fault.probability;
    // Author: OpenAI Codex, 22 October 2026. Infer each flag's readout-error
    // probability from isolated original DEM errors. Check BOTH CSS bases
    // before marginalization; a flag+opposite-basis-check fault is not an
    // isolated flag flip. Correlated multi-flag faults are not singletons.
    for (const auto& fault : faults)
        if (fault.detectors.size() == 1 && fault.observable == 0 &&
            keep_flag[fault.detectors[0]])
        {
            auto& old = measurement_p_[fault.detectors[0]];
            old = old+fault.probability-2*old*fault.probability;
        }
    for (auto f : flag_globals)
    {
        if (!(measurement_p_[f] >= 0 && measurement_p_[f] < .5))
            throw std::invalid_argument("Invalid isolated-flag DEM probability");
    }
    // Contributor: OpenAI Codex. Check the complete original detector signature
    // before discarding any opposite-basis checks. Different flag signatures
    // may distinguish different logical effects; identical full signatures cannot.
    std::map<std::vector<size_t>,std::pair<uint64_t,size_t>> signatures;
    for (size_t id = 0; id < faults.size(); id++)
    {
        const auto& fault = faults[id];
        if (fault.probability == 0) continue;
        auto support = fault.detectors;
        std::sort(support.begin(),support.end());
        const auto [it,inserted] = signatures.emplace(support,std::pair{fault.observable,id});
        if (!inserted && it->second.first != fault.observable)
        {
            std::ostringstream error;
            error << "Identical full detector support with different frame flips (distance <= 2): DEM events "
                  << it->second.second << " and " << id << ", frame masks " << it->second.first
                  << " and " << fault.observable << ", detectors={";
            for (auto v : support) error << 'D' << v << ' ';
            error << '}';
            throw std::invalid_argument(error.str());
        }
    }
    // Marginalize opposite-basis checks and same-basis extraction flags. Merge
    // independent variants only when checks, flags and logical effects agree.
    using Key = std::tuple<std::vector<size_t>, uint64_t, std::vector<size_t>>;
    std::map<Key, double> probabilities;
    std::map<Key, std::pair<BoundaryHint,bool>> provenance;
    for (const auto& fault : faults)
    {
        if (fault.probability == 0) continue;
        std::vector<size_t> checks, flags;
        for (auto v : fault.detectors)
            if (coordinates.at(v)[5] == 1)
            {
                if (keep_flag[v]) flags.push_back(v);
            }
            else if (local[v] != SIZE_MAX) checks.push_back(local[v]);
        if (checks.empty())
        {
            if (fault.observable) throw std::invalid_argument("Undetectable single logical fault");
            continue;
        }
        std::sort(checks.begin(), checks.end());
        std::sort(flags.begin(), flags.end());
        const Key key{checks,fault.observable,flags};
        if (!probabilities.contains(key)) provenance[key] = {fault.boundary_hint,fault.boundary_hint!=BoundaryHint::Unknown};
        else {
            auto& hint = provenance[key];
            if (hint.first != fault.boundary_hint) hint.first = BoundaryHint::Unknown;
            hint.second |= fault.boundary_hint != BoundaryHint::Unknown;
        }
        auto& old = probabilities[key];
        old = old+fault.probability-2*old*fault.probability;
    }
    // Contributor: OpenAI Codex. The observable is a variant label, not a key.
    std::map<std::vector<size_t>,size_t> class_ids;
    for (const auto& [key, probability] : probabilities)
    {
        const auto& [checks, observable, flags] = key;
        auto [it, inserted] = class_ids.emplace(checks,classes_.size());
        const auto [hint,explicit_hint] = provenance.at(key);
        if (inserted) {
            classes_.push_back({checks, {}, {}});
            classes_.back().boundary_hint = hint;
            classes_.back().explicit_boundary_provenance = explicit_hint;
        } else {
            auto& group = classes_[it->second];
            if (group.boundary_hint != hint) group.boundary_hint = BoundaryHint::Unknown;
            group.explicit_boundary_provenance |= explicit_hint;
        }
        classes_[it->second].variants.push_back({flags, probability, observable});
    }
    for (auto& group : classes_)
    {
        // Contributor: OpenAI Codex. For untagged circuits, only same-round
        // multi-color support supplies the spatial missing-color convention.
        // Explicit native-generator provenance takes precedence below.
        // Singleton or cross-round mixed events have ambiguous provenance and
        // keep colorless sink access in every boundary-parity sector.
        std::set<int> colors;
        const double time = coordinates.at(check_globals[group.checks[0]])[2];
        bool same_time = true;
        for (auto v : group.checks)
        {
            colors.insert(colors_[v]);
            same_time &= coordinates.at(check_globals[v])[2] == time;
        }
        if (!group.explicit_boundary_provenance) {
            if (same_time && colors.size() > 1) group.boundary_hint = BoundaryHint::Spatial;
            else if (colors.size() == 1 && group.checks.size() > 1)
                group.boundary_hint = BoundaryHint::Temporal;
        }
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

// Contributor: OpenAI Codex. A flag can identify several possible fault supports.
// Keep the union for each observed flag, without flattening detector times.
std::map<size_t,std::vector<size_t>>
CircuitDecoder::flag_neighborhoods(const std::vector<uint8_t>& syndrome) const
{
    std::map<size_t,std::vector<size_t>> neighborhoods;
    if (!use_flags_) return neighborhoods;
    for (auto f : flag_globals) if (syndrome[f]) neighborhoods[f];
    for (const auto& group : classes_)
        for (const auto& variant : group.variants)
            for (auto f : variant.flags) if (syndrome[f])
            {
                auto& detectors = neighborhoods.at(f);
                detectors.insert(detectors.end(),group.checks.begin(),group.checks.end());
            }
    for (auto& [flag,detectors] : neighborhoods)
    {
        std::sort(detectors.begin(),detectors.end());
        detectors.erase(std::unique(detectors.begin(),detectors.end()),detectors.end());
    }
    return neighborhoods;
}

CircuitDecoder::Context
CircuitDecoder::make_context(const std::vector<uint8_t>& syndrome) const
{
    Model model{colors_, {}};
    Context out;
    const auto neighborhoods = flag_neighborhoods(syndrome);
    std::vector<size_t> observed_flags;
    if (use_flags_) for (auto f : flag_globals) if (syndrome[f]) observed_flags.push_back(f);
    for (const auto& group : classes_)
    {
        // Contributor: OpenAI Codex. Select one actual member of the detector
        // equivalence class. Do not mix member frames or probabilities.
        const Variant* selected = nullptr;
        size_t selected_id = SIZE_MAX, best_distance = SIZE_MAX, active = 0;
        double unflagged = 0;
        for (auto f : group.flags) active += use_flags_ && syndrome[f];
        for (size_t id = 0; id < group.variants.size(); id++)
        {
            const auto& variant = group.variants[id];
            if (variant.flags.empty()) unflagged = std::max(unflagged,variant.probability);
            // Contributor: OpenAI Codex. Without any overlap between this
            // class's flags and the syndrome, only an unflagged representative
            // is eligible; unrelated observed flags cannot enable the class.
            if (active == 0 && !variant.flags.empty()) continue;
            size_t overlap = 0;
            if (use_flags_) for (auto f : variant.flags) overlap += syndrome[f] != 0;
            const size_t distance = observed_flags.size()+variant.flags.size()-2*overlap;
            if (!selected || distance < best_distance ||
                (distance == best_distance && variant.probability > selected->probability))
            {
                selected = &variant;
                selected_id = id;
                best_distance = distance;
            }
        }
        std::vector<size_t> incident_flags;
        for (const auto& [f,detectors] : neighborhoods)
            if (std::any_of(group.checks.begin(),group.checks.end(),[&](auto v)
                {return std::binary_search(detectors.begin(),detectors.end(),v);}))
                incident_flags.push_back(f);
        const uint64_t observable = selected ? selected->observable : 0;
        double probability = selected ? selected->probability : 0;
        double denominator = 1;
        if (selected && use_flags_)
            for (auto f : selected->flags) if (syndrome[f]) denominator *= flag_p_[f];
        if (selected)
        {
            if (!(denominator > 0)) throw std::invalid_argument("Observed representative flag has zero marginal probability");
            probability /= denominator;
            if (probability > 1+1e-12)
                throw std::invalid_argument("Conditional probability exceeds one: the product-of-flag-marginals approximation is invalid for this representative");
            // Round tiny floating-point overshoots at the probability endpoint.
            probability = std::min(probability,1.0);
        }
        out.weights.push_back({group.checks, observable, unflagged, probability,
            group.flags.size(),active,std::move(incident_flags),1.0,selected_id,best_distance,
            selected ? selected->flags : std::vector<size_t>{},selected ? selected->probability : 0,denominator});
        if (selected) model.events.push_back({probability,group.checks,observable,group.boundary_hint});
    }
    out.decoder = std::make_unique<Decoder>(std::move(model),false,boundary_gap_,ensemble_);
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
CircuitDecoder::decode(const std::vector<uint8_t>& full, std::string* debug_info)
{
    if (!debug_info) return context(full).decoder->decode(memory_syndrome(full)).observable;
    std::ostringstream out;
    out << std::setprecision(12) << "ConcatMWPM debug (natural-log weights W=-log(p))\n";
    auto label = [&](size_t id)
    {
        const auto& c = coordinates.at(id);
        std::ostringstream name;
        name << 'D' << id << '[' << (c[5] == 1 ? 'F' : "RGB"[static_cast<int>(c[3])])
             << ",t=" << c[2] << ',' << (c[4] == 0 ? 'X' : 'Z') << ']';
        return name.str();
    };
    DebugTrace trace{out,{}};
    for (auto v : check_globals) trace.detector_labels.push_back(label(v));
    const auto inspected = inspect(full,&trace);
    const auto& best = inspected.correction;
    out << "Winner: !" << "RGB"[best.color] << " member=" << best.member << " boundary_mask=" << best.boundary_mask
        << " weight=" << best.weight << " observable=" << best.observable << " correction:";
    for (auto id : best.events) out << " E" << id;
    out << "\nSelection: minimum model correction weight; first candidate wins ties.\n";
    *debug_info = out.str();
    return best.observable;
}

CircuitResult
CircuitDecoder::inspect(const std::vector<uint8_t>& full, DebugTrace* debug)
{
    auto& ctx = context(full);
    auto candidates = ctx.decoder->decode_candidates(memory_syndrome(full),true,debug);
    std::array<Result,3> colors;
    for (int c = 0; c < 3; c++) colors[c] = {0,std::numeric_limits<double>::infinity(),c,{},{},{}};
    Result best = candidates[0];
    for (const auto& candidate : candidates)
    {
        if (candidate.weight < colors[candidate.color].weight-1e-9) colors[candidate.color] = candidate;
        if (candidate.weight < best.weight-1e-9) best = candidate;
    }
    return {std::move(best), std::move(colors), ctx.weights, ctx.decoder->model,
            ctx.decoder->audit_weights(), std::move(candidates),ctx.decoder->ambiguous_boundary_events()};
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
