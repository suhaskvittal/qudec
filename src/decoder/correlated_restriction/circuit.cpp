/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Flag equivalence classes are inspired by Vittal et al., arXiv:2409.14283,
 *  Sec. VI. Measurement-error mixing follows Suhas Vittal's instructions.
 *  One synchronous fractional-weight correlation update of parent hyperedges.
 *  Adapted by OpenAI Codex from the earlier corrected flag-protocol study.
 * */
#include "decoder/correlated_restriction/circuit.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace dec::correlated_restriction
{
namespace
{
edge_type
make_edge(det_id_type x, det_id_type y)
{
    return std::minmax(x,y);
}

void
require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}
} // namespace

CircuitRestrictionDecoder::CircuitRestrictionDecoder(const stim::DetectorErrorModel& dem, Lattice lattice,
                                                   bool memory_x, double physical_p, bool use_flags,
                                                   size_t ensemble_size, bool corrected_flag_weights)
    :num_detectors_(dem.count_detectors()), physical_p_(physical_p), use_flags_(use_flags), ensemble_size_(ensemble_size),
     corrected_flag_weights_(corrected_flag_weights), flag_measurement_p_(num_detectors_,0),
     lifter_(std::move(lattice), 1), global_to_check_(num_detectors_, -1),
     canonical_to_check_(lifter_.lattice.vertices.size(), -1)
{
    require(physical_p > 0 && physical_p < 0.1, "Circuit decoder requires 0 < p < 0.1");
    require(ensemble_size > 0, "Spacetime ensemble must be nonempty");
    require(dem.count_observables() == 1, "Currently supports one logical observable");
    std::set<uint64_t> ids;
    for (size_t d = 0; d < num_detectors_; d++) ids.insert(d);
    auto coords = dem.get_detector_coordinates(ids);
    for (const auto& [id,c] : coords)
    {
        require(c.size() == 6, "Need (x,y,time,color,basis,is_flag) detector metadata");
        if (c[5] == 0)
            require(c[3] >= 0 && c[3] <= 2 && c[3] == std::floor(c[3]), "Invalid check color");
    }
    // Author: OpenAI Codex. Isolated original flag-only flips, before any
    // check-sector marginalization, define a separate readout prior per flag.
    dem.iter_flatten_error_instructions([&](const auto& instruction)
    {
        if (instruction.target_data.size() != 1) return;
        auto target = instruction.target_data[0];
        if (!target.is_relative_detector_id() || coords.at(target.val())[5] != 1) return;
        auto& old = flag_measurement_p_[target.val()];
        double p = instruction.arg_data[0];
        old = old+p-2*old*p;
    });
    double final_time = -1;
    int basis = memory_x ? 0 : 1;
    for (auto& [__, c] : coords)
    {
        require(c.size() == 6, "Need (x,y,time,color,basis,is_flag) detector metadata");
        if (c[5] == 0 && c[4] == basis) final_time = std::max(final_time, c[2]);
    }
    std::map<std::pair<double,double>, det_id_type> canonical;
    for (auto& [__, c] : coords)
        if (c[5] == 0 && c[4] == basis && c[2] == final_time)
        {
            auto v = static_cast<det_id_type>(canonical.size());
            canonical[{c[0], c[1]}] = v;
            require(lifter_.lattice.vertices[v].color == c[3], "Geometry/DEM color mismatch");
        }
    require(canonical.size() == lifter_.lattice.num_checks, "Geometry/DEM check mismatch");
    for (auto& [global, c] : coords)
    {
        if (c[5] == 1) flag_globals_.push_back(global);
        else if (c[4] == basis)
        {
            auto v = static_cast<det_id_type>(vertices_.size());
            auto base = canonical.at({c[0], c[1]});
            vertices_.push_back({static_cast<int>(c[3]), base, false});
            global_to_check_[global] = v;
            check_to_global_.push_back(global);
            // The representative is an actual detector of the SECOND
            // extraction round (time 1), where both check sectors exist.
            // Never use first-round temporal-boundary or final-readout nodes.
            if (c[2] == 1 && c[2] < final_time) canonical_to_check_[base] = v;
        }
    }
    for (size_t v = 0; v < lifter_.lattice.num_checks; v++)
        require(canonical_to_check_[v] >= 0, "Missing second-extraction-round canonical detector (need rounds >= 2)");
    for (int color = 0; color < 3; color++)
    {
        boundaries_[color] = vertices_.size();
        auto base = lifter_.lattice.boundary_id(color);
        canonical_to_check_[base] = vertices_.size();
        vertices_.push_back({color, base, true});
    }
    // Marginalize the other check sector, merging only identical
    // (memory checks, flags, logical frame) events by parity probability.
    using key_type = std::tuple<std::vector<det_id_type>, std::vector<size_t>, uint64_t>;
    std::map<key_type, double> probabilities;
    dem.iter_flatten_error_instructions([&](const auto& inst)
    {
        std::vector<det_id_type> checks;
        std::vector<size_t> flags;
        uint64_t obs = 0;
        for (auto target : inst.target_data)
        {
            require(!target.is_separator(), "Pass an UNDECOMPOSED detector error model");
            if (target.is_observable_id()) obs ^= uint64_t{1} << target.val();
            else if (target.is_relative_detector_id())
            {
                auto id = target.val();
                if (coords.at(id)[5] == 1) flags.push_back(id);
                else if (global_to_check_[id] >= 0) checks.push_back(global_to_check_[id]);
            }
        }
        std::sort(checks.begin(), checks.end());
        std::sort(flags.begin(), flags.end());
        auto& p = probabilities[{checks, flags, obs}];
        p = p*(1-inst.arg_data[0]) + (1-p)*inst.arg_data[0];
    });
    for (auto& [key, p] : probabilities)
    {
        auto& [checks, flags, obs] = key;
        errors_.push_back({p, checks, flags, obs, {}});
        if (!checks.empty()) classes_[checks].push_back(errors_.size()-1);
        else require(obs == 0, "DEM contains a single undetectable logical fault");
    }
    // Recover a spatial representative BEFORE projection erases the boundary
    // color. A single cc666 circuit fault can propagate to three data qubits.
    // Include its logical frame: equal detector supports need not be equivalent.
    const auto& spatial = lifter_.lattice;
    const auto xor_key = [](std::string a, const std::string& b)
    {
        for (size_t i = 0; i < a.size(); i++) a[i] ^= b[i];
        return a;
    };
    std::vector<std::string> face_keys;
    std::map<std::string, std::vector<size_t>> small_sections;
    small_sections.emplace(std::string(spatial.num_checks+1, 0), std::vector<size_t>{});
    for (size_t q = 0; q < spatial.faces.size(); q++)
    {
        std::string key(spatial.num_checks+1, 0);
        for (auto v : spatial.faces[q].vertices)
            if (!spatial.vertices[v].is_boundary) key[v] ^= 1;
        key.back() = spatial.faces[q].observable_mask;
        face_keys.push_back(key);
        small_sections.emplace(key, std::vector<size_t>{q});
    }
    for (size_t a = 0; a < face_keys.size(); a++)
        for (size_t b = a+1; b < face_keys.size(); b++)
            small_sections.emplace(xor_key(face_keys[a], face_keys[b]), std::vector<size_t>{a,b});
    for (auto& error : errors_)
    {
        std::string key(spatial.num_checks+1, 0);
        for (auto v : error.checks) key[vertices_[v].canonical_detector] ^= 1;
        key.back() = error.observable_mask;
        auto found = small_sections.find(key);
        if (found != small_sections.end()) error.spatial_faces = found->second;
        else
        {
            for (size_t q = 0; q < face_keys.size(); q++)
            {
                found = small_sections.find(xor_key(key, face_keys[q]));
                if (found == small_sections.end()) continue;
                error.spatial_faces = found->second;
                error.spatial_faces.push_back(q);
                break;
            }
            require(found != small_sections.end(), "DEM event needs more than three spatial data errors");
        }
    }
    projections_.resize(errors_.size());
    for (size_t id = 0; id < errors_.size(); id++)
        for (int c = 0; c < 3; c++)
        {
            std::vector<det_id_type> support;
            for (auto v : errors_[id].checks) if (vertices_[v].color != c) support.push_back(v);
            require(support.size() <= 2,"Unsupported non-graphlike projection");
            if (support.empty()) continue;
            auto x = support[0];
            auto path = event_spatial_path(errors_[id],c,vertices_[x].canonical_detector,
                support.size() == 2 ? vertices_[support[1]].canonical_detector : -1);
            auto y = support.size() == 2 ? support[1] : boundaries_[lifter_.lattice.vertices[path.back()].color];
            projections_[id][c] = {make_edge(x,vertices_[y].is_boundary ? -1 : y),{x,y},std::move(path)};
        }

}

syndrome_type
CircuitRestrictionDecoder::event_syndrome(size_t event) const
{
    syndrome_type out(num_detectors_);
    const auto& error = errors_.at(event);
    for (auto v : error.checks) out[check_to_global_[v]] ^= 1;
    for (auto f : error.flags) out[f] ^= 1;
    return out;
}

std::unique_ptr<CircuitRestrictionDecoder::Context>
CircuitRestrictionDecoder::make_context(const syndrome_type& syndrome, size_t trial) const
{
    auto context = std::make_unique<Context>();
    std::mt19937_64 rng(trial);
    std::normal_distribution<double> gaussian(0,2.0);
    std::vector<size_t> observed_flags;
    double all_flag_factor = 1;
    for (auto f : flag_globals_)
        if (syndrome[f])
        {
            observed_flags.push_back(f);
            all_flag_factor *= flag_measurement_p_[f];
        }
    for (const auto& [checks, events] : classes_)
    {
        // The representative choice can use class-local distance because the
        // outside mismatch is constant within a class. Its WEIGHT must also
        // account for every observed flag outside the class (corrected mode).
        std::set<size_t> relevant_flags;
        for (auto id : events)
            relevant_flags.insert(errors_[id].flags.begin(), errors_[id].flags.end());
        size_t active = 0;
        for (auto f : relevant_flags) active += use_flags_ && syndrome[f];
        size_t best = SIZE_MAX;
        int best_distance = std::numeric_limits<int>::max();
        double best_probability = -1;
        for (auto id : events)
        {
            const auto& e = errors_[id];
            size_t overlap = 0;
            for (auto f : e.flags) overlap += syndrome[f];
            // Flagged channels are gated off unless their flags are observed.
            if (use_flags_ && !e.flags.empty() && overlap == 0) continue;
            int distance = use_flags_ ? int(active + e.flags.size() - 2*overlap) : 0;
            if (distance < best_distance || (distance == best_distance && e.probability > best_probability))
            {
                best = id;
                best_distance = distance;
                best_probability = e.probability;
            }
        }
        if (best == SIZE_MAX) continue;
        double p = best_probability;
        double unflagged = 0, flagged = 0;
        for (auto id : events)
        {
            if (errors_[id].flags.empty()) unflagged += errors_[id].probability;
            else flagged += errors_[id].probability;
        }
        if (!use_flags_)
        {
            p = 0;
            for (auto id : events) p = p*(1-errors_[id].probability) + (1-p)*errors_[id].probability;
        }
        if (use_flags_ && active)
        {
            // For one relevant flag observed as 1 this is EXACTLY
            // (1-p_m)*Y + p_m*X. For multiple relevant flags, sum each
            // variant's probability times its independent readout likelihood.
            const double pm = 5*physical_p_;
            p = 0;
            for (auto id : events)
            {
                const auto& e = errors_[id];
                size_t overlap = 0;
                for (auto f : e.flags) overlap += syndrome[f];
                if (corrected_flag_weights_ && !e.flags.empty() && overlap == 0) continue;
                size_t mismatch = active + e.flags.size() - 2*overlap;
                double likelihood = std::pow(pm, mismatch)
                                  * std::pow(1-pm, relevant_flags.size()-mismatch);
                if (corrected_flag_weights_)
                {
                    likelihood = 1;
                    for (auto f : relevant_flags)
                    {
                        bool expected = std::binary_search(e.flags.begin(),e.flags.end(),f);
                        likelihood *= expected == bool(syndrome[f]) ? 1-flag_measurement_p_[f]
                                                                   : flag_measurement_p_[f];
                    }
                }
                p += likelihood*e.probability;
            }
        }
        if (use_flags_ && corrected_flag_weights_ && observed_flags.size() > active)
        {
            double local_factor = 1;
            for (auto f : relevant_flags) if (syndrome[f]) local_factor *= flag_measurement_p_[f];
            double outside_factor = all_flag_factor;
            if (local_factor > 0) outside_factor /= local_factor;
            else
            {
                outside_factor = 1;
                for (auto f : observed_flags)
                    if (!relevant_flags.contains(f)) outside_factor *= flag_measurement_p_[f];
            }
            p *= outside_factor;
        }
        require(p > 0 && p < 0.5, "Reweighted class probability must be between zero and one half");
        const double nominal_p = p;
        if (trial)
        {
            // Decoder-candidate diversity, not an alteration of the sampled
            // physical noise. The nominal graph remains ensemble member zero.
            double odds = p/(1-p)*std::exp(gaussian(rng));
            p = std::min(0.49, odds/(1+odds));
        }
        context->representatives[best] = p;
        context->weights.push_back({best, unflagged, flagged, nominal_p, active});
        context->enabled_flag_errors += !errors_[best].flags.empty();
        for (int ignored = 0; ignored < 3; ignored++)
        {
            const auto& projection = projections_[best][ignored];
            if (projection.key.first == -2) continue;
            const auto key = projection.key;
            auto& edge = context->edges[ignored][key];
            edge.probability = edge.probability*(1-p) + (1-edge.probability)*p;
            if (p > edge.representative_probability)
            {
                edge.representative_probability = p;
                edge.event = best;
                edge.endpoints = projection.endpoints;
                edge.spatial_path = projection.spatial_path;
            }
        }
    }
    for (int ignored = 0; ignored < 3; ignored++)
    {
        stim::DetectorErrorModel projected;
        for (const auto& [edge, data] : context->edges[ignored])
        {
            std::vector<stim::DemTarget> targets;
            if (edge.first >= 0) targets.push_back(stim::DemTarget::relative_detector_id(edge.first));
            if (edge.second >= 0) targets.push_back(stim::DemTarget::relative_detector_id(edge.second));
            projected.append_error_instruction(data.probability, targets, "");
        }
        projected += stim::DetectorErrorModel("detector D" + std::to_string(check_to_global_.size()-1));
        context->matchers[ignored] = pm::detector_error_model_to_mwpm(projected, 1 << 20, true, false);
    }
    if (trial == 0)
        for (size_t i = 1; i < ensemble_size_; i++) context->alternatives.push_back(make_context(syndrome,i));
    return context;
}

std::vector<det_id_type>
CircuitRestrictionDecoder::event_spatial_path(const ErrorEvent& error, int ignored,
                                             det_id_type x, det_id_type y) const
{
    if (x == y) return {x}; // A pure temporal edge disappears on flattening.
    const auto& lattice = lifter_.lattice;
    std::vector<std::set<det_id_type>> adjacency(lattice.vertices.size());
    for (auto q : error.spatial_faces)
    {
        const auto& face = lattice.faces[q];
        auto a = face.vertices[(ignored+1)%3], b = face.vertices[(ignored+2)%3];
        // Projection is a chain over GF(2), so shared edges cancel.
        if (!adjacency[a].erase(b)) adjacency[a].insert(b);
        if (!adjacency[b].erase(a)) adjacency[b].insert(a);
    }
    if (y < 0)
    {
        for (int c = 0; c < 3; c++)
        {
            auto boundary = lattice.boundary_id(c);
            if (adjacency[boundary].size() % 2) { y = boundary; break; }
        }
        require(y >= 0, "Spatial event is missing its colored boundary endpoint");
    }
    std::vector<det_id_type> previous(lattice.vertices.size(), -1);
    std::queue<det_id_type> queue;
    queue.push(x);
    previous[x] = x;
    while (!queue.empty() && previous[y] < 0)
    {
        auto u = queue.front(); queue.pop();
        for (auto v : adjacency[u])
            if (previous[v] < 0)
            {
                previous[v] = u;
                queue.push(v);
            }
    }
    require(previous[y] >= 0, "Spatial representative has a disconnected error projection");
    std::vector<det_id_type> path;
    for (auto v = y; ; v = previous[v])
    {
        path.push_back(v);
        if (v == x) break;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

CircuitRestrictionDecoder::MatchingOutput
CircuitRestrictionDecoder::matching_positive(Context& context, const syndrome_type& syndrome)
{
    MatchingOutput out;
    for (int ignored = 0; ignored < 3; ignored++)
    {
        std::vector<uint64_t> events;
        for (size_t v = 0; v < check_to_global_.size(); v++)
            if (syndrome[v] && vertices_[v].color != ignored) events.push_back(v);
        auto& matcher = context.matchers[ignored];
        pm::decode_detection_events_to_match_edges(matcher, events);
        const auto* nodes = matcher.flooder.graph.nodes.data();
        const auto* search_nodes = matcher.search_flooder.graph.nodes.data();
        for (const auto& m : matcher.flooder.match_edges)
        {
            auto from = m.loc_from-nodes;
            size_t to = m.loc_to ? m.loc_to-nodes : SIZE_MAX;
            std::vector<det_id_type> expanded{static_cast<det_id_type>(from)};
            matcher.search_flooder.iter_edges_on_shortest_path_from_source(from, to,
                [&](const pm::SearchGraphEdge& e)
                {
                    auto* u = e.detector_node;
                    auto* w = u->neighbors[e.neighbor_index];
                    auto x = static_cast<det_id_type>(u-search_nodes);
                    auto y = w ? static_cast<det_id_type>(w-search_nodes) : -1;
                    out.graph_edges[ignored].push_back(make_edge(x,y));
                    const auto& data = context.edges[ignored].at(make_edge(x,y));
                    if (y < 0) y = vertices_[data.endpoints.first].is_boundary
                                      ? data.endpoints.first : data.endpoints.second;
                    require(expanded.back() == x, "Shortest-path edges are not ordered");
                    auto bx = vertices_[x].canonical_detector, by = vertices_[y].canonical_detector;
                    if (bx == by) expanded.push_back(y);
                    else
                    {
                        auto path = data.spatial_path;
                        if (path.front() != bx) std::reverse(path.begin(), path.end());
                        require(path.front() == bx && path.back() == by,
                                "Stored event decomposition has inconsistent endpoints");
                        for (size_t i = 1; i < path.size(); i++)
                            expanded.push_back(i+1 == path.size() ? y : canonical_to_check_[path[i]]);
                    }
                });
            // Shortcuts expanded through colored boundaries must be split
            // BEFORE building the assignment multigraph.
            Assignment assignment{ignored, expanded[0], expanded[0], {}};
            bool interior = !vertices_[expanded[0]].is_boundary;
            for (size_t i = 1; i < expanded.size(); i++)
            {
                assignment.path.emplace_back(expanded[i-1], expanded[i]);
                assignment.y = expanded[i];
                interior |= !vertices_[expanded[i]].is_boundary;
                if (vertices_[expanded[i]].is_boundary)
                {
                    if (interior) out.assignments.push_back(assignment);
                    assignment = Assignment{ignored, expanded[i], expanded[i], {}};
                    interior = false;
                }
            }
            if (interior && !assignment.path.empty()) out.assignments.push_back(std::move(assignment));
        }
    }
    return out;
}

CircuitRestrictionDecoder::MatchingOutput
CircuitRestrictionDecoder::matching(Context& context, const syndrome_type& syndrome)
{
    bool negative = false;
    for (const auto& graph : context.edges)
        for (const auto& [key,data] : graph) negative |= data.probability > .5;
    if (!negative) return matching_positive(context,syndrome);
    MatchingOutput out;
    State selected;
    std::array<std::map<edge_type,edge_type>,3> colored;
    for (int c = 0; c < 3; c++)
    {
        std::vector<uint64_t> fired;
        for (size_t v = 0; v < check_to_global_.size(); v++)
            if (syndrome[v] && vertices_[v].color != c) fired.push_back(v);
        std::vector<int64_t> edges;
        pm::decode_detection_events_to_edges(context.matchers[c],fired,edges);
        for (size_t i = 0; i < edges.size(); i += 2)
        {
            auto e = make_edge(edges[i],edges[i+1]);
            if (!selected[c].erase(e)) selected[c].insert(e);
        }
        syndrome_type prediction(check_to_global_.size());
        for (auto e : selected[c])
        {
            const auto& data = context.edges[c].at(e);
            colored[c][e] = make_edge(data.endpoints.first,data.endpoints.second);
            if (e.first >= 0) prediction[e.first] ^= 1;
            if (e.second >= 0) prediction[e.second] ^= 1;
            out.graph_edges[c].push_back(e);
        }
        for (size_t v = 0; v < prediction.size(); v++)
            if (vertices_[v].color != c) require(prediction[v] == syndrome[v], "Restricted correction has wrong syndrome");
    }
    auto trails = decompose_edges(vertices_,syndrome,selected,colored);
    for (const auto& trail : trails)
    {
        const int ignored = trail.ignored_color;
        std::vector<det_id_type> expanded{trail.x};
        for (auto [x,y] : trail.path)
        {
            auto key = make_edge(vertices_[x].is_boundary ? -1 : x,
                                 vertices_[y].is_boundary ? -1 : y);
            const auto& data = context.edges[ignored].at(key);
            auto bx = vertices_[x].canonical_detector, by = vertices_[y].canonical_detector;
            if (bx == by) expanded.push_back(y);
            else
            {
                auto path = data.spatial_path;
                if (path.front() != bx) std::reverse(path.begin(),path.end());
                require(path.front() == bx && path.back() == by,"Inconsistent stored spatial path");
                for (size_t j = 1; j < path.size(); j++)
                    expanded.push_back(j+1 == path.size() ? y : canonical_to_check_[path[j]]);
            }
        }
        Assignment assignment{ignored,expanded[0],expanded[0],{}};
        bool interior = !vertices_[expanded[0]].is_boundary;
        for (size_t i = 1; i < expanded.size(); i++)
        {
            assignment.path.emplace_back(expanded[i-1],expanded[i]);
            assignment.y = expanded[i];
            interior |= !vertices_[expanded[i]].is_boundary;
            if (vertices_[expanded[i]].is_boundary)
            {
                if (interior) out.assignments.push_back(assignment);
                assignment = Assignment{ignored,expanded[i],expanded[i],{}};
                interior = false;
            }
        }
        if (interior && !assignment.path.empty()) out.assignments.push_back(std::move(assignment));
    }
    return out;
}

CircuitRestrictionDecoder::Context&
CircuitRestrictionDecoder::context_for(const syndrome_type& full_syndrome)
{
    require(full_syndrome.size() == num_detectors_, "Wrong full detector syndrome size");
    std::string flag_key;
    for (auto f : flag_globals_) flag_key.push_back(use_flags_ && full_syndrome[f]);
    auto it = contexts_.find(flag_key);
    if (it == contexts_.end())
    {
        // Bound memory in long Monte Carlo runs. Graphs are rebuilt exactly
        // from their flag context, so eviction cannot change decoder behavior.
        if (contexts_.size() >= 8) contexts_.clear();
        it = contexts_.emplace(flag_key, make_context(full_syndrome)).first;
    }
    return *it->second;
}

CircuitResult
CircuitRestrictionDecoder::decode(const syndrome_type& full_syndrome, std::string* debug_info)
{
    auto& context = context_for(full_syndrome);
    auto best = decode_context(context,full_syndrome,nullptr,debug_info);
    for (size_t i = 0; i < context.alternatives.size(); i++)
    {
        auto candidate = decode_context(*context.alternatives[i],full_syndrome,nullptr,debug_info);
        // Same physical correction criterion as the code-capacity ensemble.
        // Benchmark against one trial: temporal evidence can favor a different
        // logical class than the smallest flattened physical correction.
        if (candidate.lifted.correction.size() < best.lifted.correction.size())
        {
            candidate.ensemble_member = i+1;
            best = std::move(candidate);
        }
    }
    return best;
}

std::vector<CircuitResult>
CircuitRestrictionDecoder::decode_candidates(const syndrome_type& full_syndrome)
{
    auto& context = context_for(full_syndrome);
    std::vector<CircuitResult> out;
    decode_context(context,full_syndrome,&out);
    for (size_t i = 0; i < context.alternatives.size(); i++)
    {
        size_t start = out.size();
        decode_context(*context.alternatives[i],full_syndrome,&out);
        for (size_t j = start; j < out.size(); j++) out[j].ensemble_member = i+1;
    }
    return out;
}

CircuitResult
CircuitRestrictionDecoder::decode_context(Context& context, const syndrome_type& full_syndrome,
                                        std::vector<CircuitResult>* candidates, std::string* debug_info)
{
    syndrome_type residual(check_to_global_.size());
    for (size_t v = 0; v < residual.size(); v++) residual[v] = full_syndrome[check_to_global_[v]];
    // Stage 1: restricted matching on the reweighted graphs. Never consume
    // matched flag events, alter the syndrome, or apply a consensus frame.
    auto matched = matching(context, residual);
    if (debug_info) append_debug(context,matched,full_syndrome,"initial",*debug_info);
    return lift_matching(context,full_syndrome,matched,candidates);
}

CircuitResult
CircuitRestrictionDecoder::lift_matching(const Context& context, const syndrome_type& full_syndrome,
                                         const MatchingOutput& matched, std::vector<CircuitResult>* candidates)
{
    syndrome_type residual(check_to_global_.size());
    for (size_t v = 0; v < residual.size(); v++) residual[v] = full_syndrome[check_to_global_[v]];
    // Stage 2: connected components in SPACETIME, before flattening paths.
    std::array<std::vector<ConnectedComponent>, 3> components;
    for (int c = 0; c < 3; c++)
        components[c] = build_connected_components(vertices_, boundaries_, matched.assignments, c);
    // Stage 3: canonicalize each edge, retain multiplicities, and lift.
    auto flattened = matched.assignments;
    for (auto& a : flattened)
        for (auto& [x,y] : a.path)
        {
            x = vertices_[x].canonical_detector;
            y = vertices_[y].canonical_detector;
        }
    syndrome_type spatial(lifter_.lattice.num_checks);
    for (size_t v = 0; v < residual.size(); v++) spatial[vertices_[v].canonical_detector] ^= residual[v];
    restriction::Result best;
    double best_weight = std::numeric_limits<double>::infinity();
    std::array<uint64_t,3> frames{};
    std::array<size_t,3> weights{};
    for (int c = 0; c < 3; c++)
    {
        auto result = lifter_.lifting(spatial, flattened, components[c], c);
        frames[c] = result.observable_mask;
        weights[c] = result.correction.size();
        if (candidates)
            candidates->push_back({result.observable_mask,result,context.enabled_flag_errors,0,{},{}});
        double weight = result.correction.size();
        if (weight < best_weight)
        {
            best_weight = weight;
            best = std::move(result);
        }
    }
    return {best.observable_mask, std::move(best), context.enabled_flag_errors, 0, frames, weights};
}

ReweightingDiagnostics
CircuitRestrictionDecoder::inspect_reweighting(const syndrome_type& syndrome) const
{
    require(syndrome.size() == num_detectors_, "Wrong reweighting syndrome size");
    auto context = make_context(syndrome);
    double reported_pm = 5*physical_p_;
    if (corrected_flag_weights_)
    {
        reported_pm = 0;
        for (auto f : flag_globals_) reported_pm += flag_measurement_p_[f];
        if (!flag_globals_.empty()) reported_pm /= flag_globals_.size();
    }
    ReweightingDiagnostics out{reported_pm, context->weights, {}};
    for (int ignored = 0; ignored < 3; ignored++)
    {
        const auto& graph = context->matchers[ignored].flooder.graph;
        const auto& search = context->matchers[ignored].search_flooder.graph;
        const double scale = graph.normalising_constant;
        for (const auto& [edge, data] : context->edges[ignored])
        {
            auto x = edge.first < 0 ? edge.second : edge.first;
            auto y = edge.first < 0 ? -1 : edge.second;
            const auto& node = graph.nodes[x];
            const auto& search_node = search.nodes[x];
            auto index = std::find(node.neighbors.begin(), node.neighbors.end(),
                                   y < 0 ? nullptr : &graph.nodes[y])-node.neighbors.begin();
            auto search_index = std::find(search_node.neighbors.begin(), search_node.neighbors.end(),
                                          y < 0 ? nullptr : &search.nodes[y])-search_node.neighbors.begin();
            require(index < static_cast<ptrdiff_t>(node.neighbors.size()), "Graph edge was lost in PyMatching");
            require(search_index < static_cast<ptrdiff_t>(search_node.neighbors.size()), "Search edge was lost in PyMatching");
            double weight = std::log((1-data.probability)/data.probability);
            double actual = node.neighbor_weights[index]/scale;
            double actual_search = search_node.neighbor_weights[search_index]/scale;
            edge_type globals{edge.first < 0 ? -1 : check_to_global_[edge.first], check_to_global_[edge.second]};
            out.graph_edges.push_back({ignored, globals, data.probability, weight, actual, actual_search, 2/scale});
        }
    }
    return out;
}

MatchingDiagnostics
CircuitRestrictionDecoder::inspect_matching(const syndrome_type& full_syndrome)
{
    require(full_syndrome.size() == num_detectors_, "Wrong inspection syndrome size");
    auto context = make_context(full_syndrome);
    syndrome_type s(check_to_global_.size());
    for (size_t v = 0; v < s.size(); v++) s[v] = full_syndrome[check_to_global_[v]];
    auto matched = matching(*context,s);
    MatchingDiagnostics out;
    for (int c = 0; c < 3; c++)
        for (auto e : matched.graph_edges[c])
        {
            const double p = context->edges[c].at(e).probability;
            out.projected_weights[c] += std::log((1-p)/p);
            out.projected_edges[c].emplace_back(e.first < 0 ? -1 : check_to_global_[e.first], check_to_global_[e.second]);
        }
    for (size_t v = 0; v < lifter_.lattice.num_checks; v++)
        out.canonical_global_detectors.push_back(check_to_global_[canonical_to_check_[v]]);
    return out;
}

MatchingDiagnostics
CircuitRestrictionDecoder::inspect_one_pass_matching(const syndrome_type& full_syndrome)
{
    require(full_syndrome.size() == num_detectors_, "Wrong inspection syndrome size");
    auto context = make_context(full_syndrome);
    syndrome_type s(check_to_global_.size());
    for (size_t v = 0; v < s.size(); v++) s[v] = full_syndrome[check_to_global_[v]];
    auto initial = matching(*context,s);
    State evidence;
    for (int c = 0; c < 3; c++)
        for (auto e : initial.graph_edges[c])
            if (!evidence[c].erase(e)) evidence[c].insert(e);
    auto updated = update_context(*context,evidence);
    auto matched = matching(*updated,s);
    MatchingDiagnostics out;
    for (int c = 0; c < 3; c++)
        for (auto e : matched.graph_edges[c])
        {
            const double p = updated->edges[c].at(e).probability;
            out.projected_weights[c] += std::log((1-p)/p);
            out.projected_edges[c].emplace_back(e.first < 0 ? -1 : check_to_global_[e.first], check_to_global_[e.second]);
        }
    for (size_t v = 0; v < lifter_.lattice.num_checks; v++)
        out.canonical_global_detectors.push_back(check_to_global_[canonical_to_check_[v]]);
    return out;
}


std::unique_ptr<CircuitRestrictionDecoder::Context>
CircuitRestrictionDecoder::update_context(const Context& original, const State& previous) const
{
    auto updated = std::make_unique<Context>();
    for (const auto& [id,p] : original.representatives)
        for (int target = 0; target < 3; target++)
        {
            const auto& projection = projections_[id][target];
            if (projection.key.first == -2) continue;
            double q = p;
            const size_t k = errors_[id].checks.size();
            if (k >= 3)
            {
                std::set<edge_type> evidence;
                for (int other = 0; other < 3; other++)
                {
                    const auto key = projections_[id][other].key;
                    if (other != target && key.first >= 0 && key.second >= 0 && previous[other].contains(key))
                        evidence.insert(key);
                }
                if (!evidence.empty())
                {
                    const double prior = std::clamp(p,1e-300,1-1e-12);
                    const double weight = std::log1p(-prior)-std::log(prior);
                    const double pairs = double(k)*(k-1)/2;
                    const double updated_weight = (1-double(evidence.size())/pairs)*weight;
                    q = 1/(1+std::exp(updated_weight));
                }
            }
            auto& data = updated->edges[target][projection.key];
            data.probability = data.probability*(1-q)+(1-data.probability)*q;
            if (q > data.representative_probability)
            {
                data.representative_probability = q;
                data.event = id;
                data.endpoints = projection.endpoints;
                data.spatial_path = projection.spatial_path;
            }
        }
    updated->enabled_flag_errors = original.enabled_flag_errors;
    for (int c = 0; c < 3; c++)
    {
        stim::DetectorErrorModel dem;
        for (const auto& [key,data] : updated->edges[c])
        {
            std::vector<stim::DemTarget> targets;
            if (key.first >= 0) targets.push_back(stim::DemTarget::relative_detector_id(key.first));
            if (key.second >= 0) targets.push_back(stim::DemTarget::relative_detector_id(key.second));
            dem.append_error_instruction(std::clamp(data.probability,1e-300,1-1e-12),targets,"");
        }
        dem += stim::DetectorErrorModel("detector D"+std::to_string(check_to_global_.size()-1));
        updated->matchers[c] = pm::detector_error_model_to_mwpm(dem,1 << 20,true,false);
    }
    return updated;
}

std::array<Variant,2>
CircuitRestrictionDecoder::decode_single_pass(const syndrome_type& full_syndrome)
{
    auto& original = context_for(full_syndrome);
    syndrome_type residual(check_to_global_.size());
    for (size_t v = 0; v < residual.size(); v++) residual[v] = full_syndrome[check_to_global_[v]];
    auto first = matching(original,residual);
    State previous;
    for (int c = 0; c < 3; c++)
        for (auto e : first.graph_edges[c])
            if (!previous[c].erase(e)) previous[c].insert(e);
    auto updated = update_context(original,previous);
    auto second = matching(*updated,residual);
    std::array<Variant,2> out;
    auto evaluate = [&](const Context& context, const MatchingOutput& matched, Variant& result)
    {
        result.score = std::numeric_limits<double>::infinity();
        try
        {
            result.frame = lift_matching(context,full_syndrome,matched).observable_mask;
            result.score = 0;
            for (int c = 0; c < 3; c++)
                for (auto edge : matched.graph_edges[c])
                {
                    const double p = original.edges[c].at(edge).probability;
                    result.score += std::log((1-p)/p);
                }
        }
        catch (const std::runtime_error&) {}
    };
    evaluate(original,first,out[0]);
    evaluate(*updated,second,out[1]);
    out[0].termination = "baseline";
    out[1].termination = "one_pass";
    out[1].iterations = 1;
    return out;
}

CircuitResult
CircuitRestrictionDecoder::decode_correlated(const syndrome_type& full_syndrome, std::string* debug_info)
{
    require(full_syndrome.size() == num_detectors_, "Wrong syndrome size");
    auto& original = context_for(full_syndrome);
    syndrome_type residual(check_to_global_.size());
    for (size_t v = 0; v < residual.size(); v++) residual[v] = full_syndrome[check_to_global_[v]];
    const auto initial = matching(original,residual);
    if (debug_info) append_debug(original,initial,full_syndrome,"initial",*debug_info);
    State evidence;
    for (int c = 0; c < 3; c++)
        for (auto e : initial.graph_edges[c])
            if (!evidence[c].erase(e)) evidence[c].insert(e);
    auto updated = update_context(original,evidence);
    const auto second = matching(*updated,residual);
    if (debug_info) append_debug(*updated,second,full_syndrome,"correlated one-pass",*debug_info);
    return lift_matching(*updated,full_syndrome,second);
}

// Contributor: OpenAI Codex. Diagnostics from the actual matching pass.
void CircuitRestrictionDecoder::append_debug(const Context& context,const MatchingOutput& matched,
    const syndrome_type&,const char* stage,std::string& debug) const
{
    std::ostringstream out;
    out << std::setprecision(12) << "Restriction " << stage << " matching and weights\n";
    auto name = [&](det_id_type v)
    {
        if (v < 0) return std::string("BoundarySink");
        if (vertices_.at(v).is_boundary) return std::string(1,"RGB"[vertices_[v].color])+"Boundary";
        return "D"+std::to_string(check_to_global_.at(v));
    };
    for (int c = 0; c < 3; c++)
    {
        out << "  !" << "RGB"[c] << " (" << "RGB"[(c+1)%3] << "RGB"[(c+2)%3] << ") selected graph edges:\n";
        for (const auto& edge : matched.graph_edges[c])
        {
            const double p = context.edges[c].at(edge).probability;
            out << "    " << name(edge.first) << " -- " << name(edge.second)
                << " p=" << p << " W=" << std::log((1-p)/p) << '\n';
        }
        for (const auto& assignment : matched.assignments) if (assignment.ignored_color == c)
        {
            out << "    match " << name(assignment.x) << " <-> " << name(assignment.y) << " expanded edges:";
            for (const auto& [x,y] : assignment.path) out << " (" << name(x) << ',' << name(y) << ')';
            out << '\n';
        }
    }
    debug += out.str();
}

} // namespace dec::correlated_restriction
