/*
 *  author: OpenAI Codex
 *  date:   21 October 2026
 *
 *  Algorithm references: Kubica--Delfosse, arXiv:1905.07393;
 *  triangular-boundary assignment components follow the construction used
 *  by Suhas Vittal in qontra/src/qontra/decoder/restriction.cpp.
 *  This is an independent implementation, using vendored PyMatching.
 * */
#include "decoder/restriction/lattice.h"
#include "circuit_generator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>

namespace dec::restriction
{

namespace
{
edge_type
make_edge(det_id_type x, det_id_type y)
{
    return std::minmax(x, y);
}

void
require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
} // namespace

det_id_type
Lattice::boundary_id(int color) const
{
    return static_cast<det_id_type>(num_checks + color);
}

syndrome_type
Lattice::syndrome(const std::vector<size_t>& errors) const
{
    syndrome_type out(num_checks);
    for (auto q : errors)
        for (auto v : faces.at(q).vertices)
            if (!vertices[v].is_boundary) out[v] ^= 1;
    return out;
}

uint64_t
Lattice::observable(const std::vector<size_t>& errors) const
{
    uint64_t out = 0;
    for (auto q : errors) out ^= faces.at(q).observable_mask;
    return out;
}

Lattice
make_code_capacity_lattice(bool octagons, uint32_t distance, double p)
{
    require(std::isfinite(p) && p > 0 && p < 0.5, "Need 0 < p < 0.5");
    auto circuit = octagons ? cc488_si1000(distance, 1, 0, false)
                           : cc666_si1000(distance, 1, 0, false);
    const size_t n = octagons ? (distance*distance + 2*distance - 1)/2
                             : (3*distance*distance + 1)/4;
    const size_t m = (n - 1)/2;
    std::set<uint64_t> ids;
    for (uint64_t i = 0; i < circuit.count_detectors(); i++) ids.insert(i);
    auto coordinates = circuit.get_detector_coordinates(ids);
    Lattice out;
    out.num_checks = m;
    out.faces.resize(n);
    for (auto& f : out.faces)
    {
        f.vertices = {static_cast<det_id_type>(m), static_cast<det_id_type>(m+1),
                      static_cast<det_id_type>(m+2)};
        f.probability = p;
    }
    size_t detector = 0;
    circuit.for_each_operation([&](const stim::CircuitInstruction& op)
    {
        if (op.gate_type != stim::GateType::DETECTOR) return;
        const auto& c = coordinates.at(detector++);
        if (c[2] != 1 || c[5] != 0) return;
        auto v = static_cast<det_id_type>(out.vertices.size());
        int color = static_cast<int>(c[3]);
        out.vertices.push_back({color, v, false});
        for (auto t : op.targets)
            if (t.rec_offset() >= -static_cast<int>(n))
            {
                auto q = n + t.rec_offset();
                require(out.faces[q].vertices[color] == out.boundary_id(color),
                        "Two same-color checks touch a qubit");
                out.faces[q].vertices[color] = v;
            }
    });
    require(out.vertices.size() == m, "Missing check supports");
    for (int c = 0; c < 3; c++) out.vertices.push_back({c, out.boundary_id(c), true});
    out.incident_faces.resize(m+3);
    for (size_t q = 0; q < n; q++)
        for (auto v : out.faces[q].vertices) out.incident_faces[v].push_back(q);
    return out;
}

RestrictionDecoder::RestrictionDecoder(Lattice lat, size_t ensemble_size, uint64_t seed)
    :lattice(std::move(lat)), sections_(lattice.num_checks)
{
    require(ensemble_size > 0, "Ensemble must be nonempty");
    // First build the three restricted decoding graphs. A graph edge is the
    // projection of a physical face; duplicate projections combine by XOR.
    for (int ignored = 0; ignored < 3; ignored++)
    {
        std::map<edge_type, double> probabilities;
        for (const auto& f : lattice.faces)
        {
            auto x = f.vertices[(ignored+1)%3], y = f.vertices[(ignored+2)%3];
            bool bx = lattice.vertices[x].is_boundary, by = lattice.vertices[y].is_boundary;
            if (bx && by) continue;
            if (bx) std::swap(x, y);
            auto graph_edge = make_edge(x, (bx || by) ? -1 : y);
            if (bx || by) boundary_edges_[ignored][graph_edge] = make_edge(x, y);
            auto& p = probabilities[graph_edge];
            p = p*(1-f.probability) + (1-p)*f.probability;
        }
        stim::DetectorErrorModel dem;
        for (auto [edge, p] : probabilities)
        {
            std::vector<stim::DemTarget> targets;
            if (edge.first >= 0) targets.push_back(stim::DemTarget::relative_detector_id(edge.first));
            if (edge.second >= 0) targets.push_back(stim::DemTarget::relative_detector_id(edge.second));
            dem.append_error_instruction(p, targets, "");
        }
        // Keep detector IDs identical across projections, including isolates.
        dem += stim::DetectorErrorModel("detector D" + std::to_string(lattice.num_checks-1));
        matchers_[ignored] = pm::detector_error_model_to_mwpm(dem, 1 << 20, true, false);
    }
    // Precompute the minimum-weight solution for every realizable local
    // boundary. The two sections of a star differ by a stabilizer.
    for (size_t v = 0; v < lattice.num_checks; v++)
    {
        const auto& faces = lattice.incident_faces[v];
        require(faces.size() < 20, "Local star unexpectedly large");
        std::map<std::vector<edge_type>, LocalSection> best;
        for (uint64_t mask = 0; mask < (uint64_t{1} << faces.size()); mask++)
        {
            LocalSection section;
            section.weight = 0;
            std::set<edge_type> edges;
            for (size_t i = 0; i < faces.size(); i++)
                if ((mask >> i) & 1)
                {
                    size_t q = faces[i];
                    section.faces.push_back(q);
                    section.weight += std::log((1-lattice.faces[q].probability)/lattice.faces[q].probability);
                    for (auto u : lattice.faces[q].vertices)
                        if (u != static_cast<det_id_type>(v))
                        {
                            auto e = make_edge(v, u);
                            if (!edges.erase(e)) edges.insert(e);
                        }
                }
            section.boundary.assign(edges.begin(), edges.end());
            auto it = best.find(section.boundary);
            if (it == best.end() || section.weight < it->second.weight)
                best[section.boundary] = section;
        }
        for (auto& [__, section] : best) sections_[v].push_back(std::move(section));
    }
    // Candidate diversity is generated only through restricted matching and
    // local lifting. No exact decoder or global syndrome-solver fallback.
    // Scores at the end use the ORIGINAL physical error probabilities.
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> gaussian(0, 2.0);
    for (size_t trial = 1; trial < ensemble_size; trial++)
    {
        auto candidate = lattice;
        for (auto& face : candidate.faces)
        {
            double odds = face.probability/(1-face.probability);
            odds *= std::exp(gaussian(rng));
            face.probability = std::min(0.49, odds/(1+odds));
        }
        alternatives_.push_back(std::make_unique<RestrictionDecoder>(std::move(candidate), 1));
    }
}

std::vector<Assignment>
RestrictionDecoder::restricted_lattice_matching(const syndrome_type& syndrome)
{
    require(syndrome.size() == lattice.num_checks, "Wrong syndrome size");
    std::vector<Assignment> assignments;
    for (int ignored = 0; ignored < 3; ignored++)
    {
        std::vector<uint64_t> events;
        for (size_t v = 0; v < syndrome.size(); v++)
            if (syndrome[v] && lattice.vertices[v].color != ignored) events.push_back(v);
        auto& matcher = matchers_[ignored];
        pm::decode_detection_events_to_match_edges(matcher, events);
        const auto* nodes = matcher.flooder.graph.nodes.data();
        const auto* search_nodes = matcher.search_flooder.graph.nodes.data();
        for (const auto& m : matcher.flooder.match_edges)
        {
            size_t from = m.loc_from - nodes, to = m.loc_to ? m.loc_to - nodes : SIZE_MAX;
            Assignment assignment{ignored, static_cast<det_id_type>(from), -1, {}};
            matcher.search_flooder.iter_edges_on_shortest_path_from_source(from, to,
                [&](const pm::SearchGraphEdge& e)
                {
                    auto* u = e.detector_node;
                    auto* w = u->neighbors[e.neighbor_index];
                    auto edge = make_edge(u-search_nodes, w ? w-search_nodes : -1);
                    assignment.path.push_back(w ? edge : boundary_edges_[ignored].at(edge));
                });
            require(!assignment.path.empty(), "Empty matching path");
            if (m.loc_to) assignment.y = to;
            else
            {
                const auto& last = assignment.path.back();
                assignment.y = lattice.vertices[last.first].is_boundary ? last.first : last.second;
                require(lattice.vertices[assignment.y].is_boundary, "Missing colored boundary");
            }
            assignments.push_back(std::move(assignment));
        }
    }
    return assignments;
}

std::vector<ConnectedComponent>
RestrictionDecoder::connected_components(const std::vector<Assignment>& assignments, int reference_color) const
{
    return build_connected_components(lattice.vertices,
        {lattice.boundary_id(0), lattice.boundary_id(1), lattice.boundary_id(2)}, assignments, reference_color);
}

std::vector<ConnectedComponent>
build_connected_components(const std::vector<Vertex>& vertices, const std::array<det_id_type, 3>& boundaries,
                           const std::vector<Assignment>& assignments, int reference_color)
{
    require(reference_color >= 0 && reference_color < 3, "Invalid reference color");
    // Multigraph adjacency stores ASSIGNMENT IDs, never just endpoint pairs.
    // Every fired interior vertex has degree two (one match in each lattice
    // containing its color). Boundaries can have arbitrarily many matches.
    std::vector<std::vector<size_t>> adjacency(vertices.size());
    for (size_t i = 0; i < assignments.size(); i++)
    {
        adjacency[assignments[i].x].push_back(i);
        adjacency[assignments[i].y].push_back(i);
    }
    for (size_t v = 0; v < vertices.size(); v++)
        require(vertices[v].is_boundary || adjacency[v].empty() || adjacency[v].size() == 2,
                "Assignment multigraph violates degree-two interior invariant");
    std::vector<ConnectedComponent> components;
    std::set<size_t> visited;
    const auto red_boundary = boundaries[reference_color];
    for (size_t first : adjacency[red_boundary])
    {
        if (visited.contains(first)) continue;
        ConnectedComponent component;
        auto current = red_boundary;
        auto next_edge = first;
        do
        {
            require(visited.insert(next_edge).second, "Assignment reused during boundary traversal");
            component.assignments.push_back(next_edge);
            const auto& a = assignments[next_edge];
            current = a.x == current ? a.y : a.x;
            if (vertices[current].is_boundary) break;
            const auto& incident = adjacency[current];
            next_edge = incident[0] == next_edge ? incident[1] : incident[0];
        } while (true);
        const int other = vertices[current].color;
        component.color = other == reference_color ? (reference_color+1)%3 : 3-other-reference_color;
        components.push_back(std::move(component));
    }
    return components;
}

Result
RestrictionDecoder::lifting(const syndrome_type& syndrome, const std::vector<Assignment>& assignments,
                            const std::vector<ConnectedComponent>& components, int reference_color) const
{
    Result out;
    out.diagnostics.assignments = assignments.size();
    out.diagnostics.components = components.size();
    std::vector<uint8_t> correction(lattice.faces.size());
    std::vector<uint8_t> in_component(assignments.size());
    const auto lift_group = [&](const std::vector<size_t>& group, int color)
    {
        std::map<edge_type, size_t> incidence;
        for (auto i : group)
            for (auto [x, y] : assignments[i].path)
            {
                x = lattice.vertices[x].canonical_detector;
                y = lattice.vertices[y].canonical_detector;
                if (x == y) continue; // Pure time-like step.
                if (lattice.vertices[x].is_boundary && lattice.vertices[y].is_boundary) continue;
                if (lattice.vertices[x].color == color || lattice.vertices[y].color == color)
                    out.diagnostics.max_multiplicity = std::max(out.diagnostics.max_multiplicity,
                                                               ++incidence[make_edge(x,y)]);
            }
        std::map<det_id_type, std::vector<edge_type>> stars;
        for (auto [e, count] : incidence)
            if (count % 2)
            {
                auto v = lattice.vertices[e.first].color == color ? e.first : e.second;
                require(!lattice.vertices[v].is_boundary, "Attempted to lift a boundary star");
                stars[v].push_back(e);
            }
        for (const auto& [v, boundary] : stars)
        {
            const auto& sections = sections_[v];
            auto it = std::find_if(sections.begin(), sections.end(),
                                  [&](const auto& s) { return s.boundary == boundary; });
            require(it != sections.end(), "Local lifting boundary is not realizable");
            for (auto q : it->faces) correction[q] ^= 1;
            ++out.diagnostics.local_lifts;
        }
    };
    for (const auto& component : components)
    {
        for (auto i : component.assignments) in_component[i] = 1;
        lift_group(component.assignments, component.color);
    }
    std::vector<size_t> bulk;
    for (size_t i = 0; i < assignments.size(); i++)
        if (!in_component[i] && assignments[i].ignored_color != reference_color) bulk.push_back(i);
    lift_group(bulk, reference_color);
    // Overlapping component lifts can leave stabilizer factors. Remove local
    // weight-reducing stabilizers before comparing candidates. This changes
    // neither the syndrome nor the logical class of any lifted correction.
    bool changed;
    do
    {
        changed = false;
        for (size_t v = 0; v < lattice.num_checks; v++)
        {
            double delta = 0;
            for (auto q : lattice.incident_faces[v])
            {
                double w = std::log((1-lattice.faces[q].probability)/lattice.faces[q].probability);
                delta += correction[q] ? -w : w;
            }
            if (delta < -1e-9)
            {
                for (auto q : lattice.incident_faces[v]) correction[q] ^= 1;
                changed = true;
            }
        }
    } while (changed);
    for (size_t q = 0; q < correction.size(); q++)
        if (correction[q]) out.correction.push_back(q);
    out.observable_mask = lattice.observable(out.correction);
    require(lattice.syndrome(out.correction) == syndrome, "Lifted correction does not reproduce the syndrome");
    return out;
}

Result
RestrictionDecoder::decode(const syndrome_type& syndrome)
{
    auto assignments = restricted_lattice_matching(syndrome);
    std::array<std::vector<ConnectedComponent>, 3> components;
    for (int color = 0; color < 3; color++) components[color] = connected_components(assignments, color);
    Result best;
    double best_weight = std::numeric_limits<double>::infinity();
    for (int color = 0; color < 3; color++)
    {
        auto result = lifting(syndrome, assignments, components[color], color);
        double weight = 0;
        for (auto q : result.correction)
            weight += std::log((1-lattice.faces[q].probability)/lattice.faces[q].probability);
        if (weight < best_weight)
        {
            best_weight = weight;
            best = std::move(result);
        }
    }
    for (auto& alternative : alternatives_)
    {
        auto result = alternative->decode(syndrome);
        double weight = 0;
        for (auto q : result.correction)
            weight += std::log((1-lattice.faces[q].probability)/lattice.faces[q].probability);
        if (weight < best_weight)
        {
            best_weight = weight;
            best = std::move(result);
        }
    }
    return best;
}

} // namespace dec::restriction
