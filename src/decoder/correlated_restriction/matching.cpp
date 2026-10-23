/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Matching uses vendored PyMatching; components and lifting reuse qudec.
 *  Fractional log-odds discount from distinct pairs selected in other lattices.
 * */
#include "decoder/correlated_restriction/matching.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace dec::correlated_restriction
{
namespace
{
constexpr edge_type absent{-2,-2};
edge_type
edge(det_id_type x, det_id_type y)
{
    return std::minmax(x,y);
}

void
require(bool ok, const char* msg)
{
    if (!ok) throw std::runtime_error(msg);
}
}

Decoder::Decoder(uint32_t d, uint32_t r, double p, bool correlated)
    :lattice(make_code_capacity_lattice(false,d,p)), rounds(r),
     num_detectors(lattice.num_checks*(r == 1 ? 1 : r+1)), correlated_(correlated),
     lifter_(lattice,1)
{
    require(r == 1 || r == d, "Use one capacity round or d phenomenological rounds");
    const size_t m = lattice.num_checks;
    for (size_t i = 0; i < num_detectors; i++)
        vertices.push_back({lattice.vertices[i%m].color,static_cast<det_id_type>(i%m),false});
    for (int c = 0; c < 3; c++)
    {
        boundaries[c] = vertices.size();
        vertices.push_back({c,lattice.boundary_id(c),true});
    }
    auto append = [&](double probability, std::vector<det_id_type> checks, uint64_t obs)
    {
        Event event{probability,std::move(checks),obs,{absent,absent,absent}};
        for (int ignored = 0; ignored < 3; ignored++)
        {
            std::vector<det_id_type> projected;
            for (auto v : event.checks)
                if (vertices[v].color != ignored) projected.push_back(v);
            if (projected.empty()) continue;
            require(projected.size() <= 2, "Non-graphlike event");
            det_id_type x = projected[0], y;
            if (projected.size() == 2) y = projected[1];
            else
            {
                int missing = 3-ignored-vertices[x].color;
                y = boundaries[missing];
            }
            auto key = edge(x,vertices[y].is_boundary ? -1 : y);
            event.projection[ignored] = key;
            colored_edges_[ignored][key] = edge(x,y);
            incident_events_[ignored][key].push_back(events.size());
        }
        events.push_back(std::move(event));
    };
    for (uint32_t t = 0; t < rounds; t++)
        for (const auto& face : lattice.faces)
        {
            std::vector<det_id_type> checks;
            for (auto v : face.vertices)
                if (!lattice.vertices[v].is_boundary) checks.push_back(t*m+v);
            append(p,std::move(checks),face.observable_mask);
        }
    if (rounds > 1)
        for (uint32_t t = 0; t < rounds; t++)
            for (size_t v = 0; v < m; v++) append(p,{det_id_type(t*m+v),det_id_type((t+1)*m+v)},0);
    // Build the three initial graphs once. Updates always start from events[].probability.
    syndrome_type empty(num_detectors);
    match(empty,nullptr,true);
}

syndrome_type
Decoder::syndrome(const std::vector<size_t>& faults) const
{
    syndrome_type s(num_detectors);
    for (auto f : faults) for (auto v : events.at(f).checks) s[v] ^= 1;
    return s;
}

uint64_t
Decoder::observable(const std::vector<size_t>& faults) const
{
    uint64_t result = 0;
    for (auto f : faults) result ^= events.at(f).observable;
    return result;
}

std::vector<std::array<double,3>>
Decoder::conditional_probabilities(const State& previous) const
{
    std::vector<std::array<double,3>> probabilities(events.size());
    for (size_t i = 0; i < events.size(); i++)
        for (int target = 0; target < 3; target++)
        {
            probabilities[i][target] = events[i].probability;
            // Only parent hyperedges are updated. Measurement faults and
            // spatial boundary events with one/two checks keep their priors.
            if (events[i].checks.size() < 3) continue;
            std::set<edge_type> evidence;
            for (int other = 0; other < 3; other++)
            {
                const auto key = events[i].projection[other];
                if (other != target && key.first >= 0 && key.second >= 0 && previous[other].contains(key))
                    evidence.insert(key);
            }
            if (evidence.empty()) continue;
            const size_t k = events[i].checks.size();
            const double pairs = double(k)*(k-1)/2;
            const double p = std::clamp(events[i].probability,probability_floor,1-probability_floor);
            const double weight = std::log1p(-p)-std::log(p);
            const double updated = (1-double(evidence.size())/pairs)*weight;
            probabilities[i][target] = 1/(1+std::exp(updated));
        }
    return probabilities;
}

State
Decoder::match(const syndrome_type& syndrome,
               const std::vector<std::array<double,3>>* probabilities, bool initialize)
{
    State out;
    for (int c = 0; c < 3; c++)
    {
        std::unique_ptr<pm::Mwpm> updated;
        if (initialize || probabilities)
        {
            stim::DetectorErrorModel dem;
            for (const auto& [key,ids] : incident_events_[c])
            {
                double p = 0;
                for (auto id : ids)
                {
                    const double q = probabilities ? (*probabilities)[id][c] : events[id].probability;
                    p = p*(1-q)+(1-p)*q;
                }
                // Finite log odds for exact posterior 0/1. PyMatching's edge
                // extraction includes negative edges, unlike its match-pair API.
                p = std::clamp(p,probability_floor,1-probability_floor);
                std::vector<stim::DemTarget> targets;
                if (key.first >= 0) targets.push_back(stim::DemTarget::relative_detector_id(key.first));
                if (key.second >= 0) targets.push_back(stim::DemTarget::relative_detector_id(key.second));
                dem.append_error_instruction(p,targets,"");
            }
            dem += stim::DetectorErrorModel("detector D"+std::to_string(num_detectors-1));
            auto matcher = pm::detector_error_model_to_mwpm(dem,1 << 20,true,false);
            if (initialize) initial_matchers_[c] = std::move(matcher);
            else updated = std::make_unique<pm::Mwpm>(std::move(matcher));
        }
        if (initialize) continue;
        std::vector<uint64_t> fired;
        for (size_t v = 0; v < num_detectors; v++)
            if (syndrome[v] && vertices[v].color != c) fired.push_back(v);
        std::vector<int64_t> selected;
        auto& matcher = updated ? *updated : initial_matchers_[c];
        pm::decode_detection_events_to_edges(matcher,fired,selected);
        for (size_t i = 0; i < selected.size(); i += 2)
        {
            auto key = edge(selected[i],selected[i+1]);
            if (!out[c].erase(key)) out[c].insert(key);
        }
        syndrome_type predicted(num_detectors);
        for (auto [x,y] : out[c])
        {
            if (x >= 0) predicted[x] ^= 1;
            if (y >= 0) predicted[y] ^= 1;
        }
        for (size_t v = 0; v < num_detectors; v++)
            if (vertices[v].color != c) require(predicted[v] == syndrome[v], "Matching has wrong boundary");
    }
    return out;
}

std::vector<Assignment>
decompose_edges(const std::vector<Vertex>& vertices, const syndrome_type& syndrome,
                const State& state, const std::array<std::map<edge_type,edge_type>,3>& colored_edges)
{
    std::vector<Assignment> assignments;
    for (int c = 0; c < 3; c++)
    {
        // Decompose the parity edge set into trails. A defect has one open
        // stub; each colored boundary stub is open. Pair other stubs locally.
        const size_t first_assignment = assignments.size();
        std::vector<edge_type> edges;
        std::vector<std::vector<size_t>> adjacency(vertices.size());
        for (auto key : state[c])
        {
            auto e = colored_edges[c].at(key);
            adjacency[e.first].push_back(edges.size());
            adjacency[e.second].push_back(edges.size());
            edges.push_back(e);
        }
        std::map<std::pair<det_id_type,size_t>,size_t> next;
        std::vector<std::pair<det_id_type,size_t>> starts;
        for (size_t v = 0; v < adjacency.size(); v++)
        {
            const auto& incident = adjacency[v];
            size_t first = 0;
            if (vertices[v].is_boundary)
            {
                for (auto i : incident) starts.emplace_back(v,i);
                continue;
            }
            if (syndrome[v] && vertices[v].color != c)
            {
                require(incident.size()%2 == 1, "Defect has even degree");
                starts.emplace_back(v,incident[0]);
                first = 1;
            }
            require((incident.size()-first)%2 == 0, "Unexpected unpaired path stub");
            for (size_t i = first; i < incident.size(); i += 2)
            {
                next[{v,incident[i]}] = incident[i+1];
                next[{v,incident[i+1]}] = incident[i];
            }
        }
        std::vector<bool> used(edges.size());
        auto traverse = [&](det_id_type start, size_t first)
        {
            Assignment a{c,start,start,{}};
            auto v = start;
            size_t current = first;
            while (!used[current])
            {
                used[current] = true;
                auto e = edges[current];
                auto w = e.first == v ? e.second : e.first;
                a.path.emplace_back(v,w);
                a.y = w;
                auto it = next.find({w,current});
                if (it == next.end()) break;
                v = w;
                current = it->second;
            }
            assignments.push_back(std::move(a));
        };
        for (auto [v,i] : starts) if (!used[i]) traverse(v,i);
        for (size_t i = 0; i < edges.size(); i++)
            if (!used[i])
            {
                traverse(edges[i].first,i);
                auto cycle = std::move(assignments.back());
                assignments.pop_back();
                require(cycle.x == cycle.y, "Residual trail is not closed");
                bool spliced = false;
                // A cycle touching an open trail belongs in that trail, not
                // in a new assignment at an already matched syndrome vertex.
                // Otherwise the assignment multigraph would gain degree four.
                for (size_t a = first_assignment; a < assignments.size() && !spliced; a++)
                {
                    auto& path = assignments[a].path;
                    for (size_t k = 0; k <= path.size() && !spliced; k++)
                    {
                        auto v = k == path.size() ? path.back().second : path[k].first;
                        for (size_t j = 0; j < cycle.path.size(); j++)
                            if (cycle.path[j].first == v)
                            {
                                std::rotate(cycle.path.begin(),cycle.path.begin()+j,cycle.path.end());
                                path.insert(path.begin()+k,cycle.path.begin(),cycle.path.end());
                                spliced = true;
                                break;
                            }
                    }
                }
                if (!spliced) assignments.push_back(std::move(cycle));
            }
    }
    return assignments;
}

std::vector<Assignment>
Decoder::assignments(const syndrome_type& syndrome, const State& state) const
{
    return decompose_edges(vertices,syndrome,state,colored_edges_);
}

// Same time-history coordinate descent as the prior phenomenological study.
// It preserves a lifted logical class and validates the entire spacetime syndrome.
double
Decoder::score_history(const syndrome_type& syndrome, const restriction::Result& candidate) const
{
    const size_t n = lattice.faces.size(), m = lattice.num_checks;
    const size_t T = rounds;
    std::vector<uint8_t> data(n*T), measurement(m*T);
    for (auto q : candidate.correction) data[q*T] = 1;
    for (size_t v = 0; v < m; v++)
    {
        uint8_t previous = 0;
        for (size_t t = 0; t < T; t++)
        {
            uint8_t value = syndrome[t*m+v]^previous;
            for (auto q : lattice.incident_faces[v]) value ^= data[q*T+t];
            measurement[v*T+t] = value;
            previous = value;
        }
        require(previous == syndrome[T*m+v], "History cannot close final temporal boundary");
    }
    bool changed;
    do
    {
        changed = false;
        for (size_t q = 0; q < n; q++)
        {
            std::vector<uint8_t> old_prefix(T);
            std::vector<std::array<uint8_t,2>> parent(T);
            uint8_t parity = 0;
            double old_cost = 0;
            for (size_t t = 0; t < T; t++)
            {
                parity ^= data[q*T+t]; old_prefix[t] = parity;
                old_cost += data[q*T+t];
                for (auto v : lattice.faces[q].vertices)
                    if (!lattice.vertices[v].is_boundary) old_cost += measurement[v*T+t];
            }
            std::array<double,2> costs{0,std::numeric_limits<double>::infinity()};
            for (size_t t = 0; t < T; t++)
            {
                std::array<double,2> next{};
                for (int state = 0; state < 2; state++)
                {
                    double cost = 0;
                    for (auto v : lattice.faces[q].vertices)
                        if (!lattice.vertices[v].is_boundary) cost += measurement[v*T+t]^old_prefix[t]^state;
                    int prior = costs[state] <= costs[1-state]+1 ? state : 1-state;
                    parent[t][state] = prior;
                    next[state] = costs[prior]+(prior != state)+cost;
                }
                costs = next;
            }
            if (costs[parity] >= old_cost-1e-9) continue;
            std::vector<uint8_t> prefix(T);
            uint8_t state = parity;
            for (size_t t = T; t-- > 0; ) { prefix[t] = state; state = parent[t][state]; }
            for (size_t t = 0; t < T; t++)
            {
                data[q*T+t] = prefix[t]^(t ? prefix[t-1] : 0);
                if (prefix[t] != old_prefix[t])
                    for (auto v : lattice.faces[q].vertices)
                        if (!lattice.vertices[v].is_boundary) measurement[v*T+t] ^= 1;
            }
            changed = true;
        }
        for (size_t v = 0; v < m; v++)
            for (size_t t = 0; t < T; t++)
            {
                int delta = 0;
                for (auto q : lattice.incident_faces[v]) delta += data[q*T+t] ? -1 : 1;
                if (delta >= 0) continue;
                for (auto q : lattice.incident_faces[v]) data[q*T+t] ^= 1;
                changed = true;
            }
    } while (changed);
    uint64_t frame = 0;
    for (auto bit : data) frame ^= bit;
    require(frame == candidate.observable_mask, "History changed logical class");
    for (size_t v = 0; v < m; v++)
        for (size_t t = 0; t <= T; t++)
        {
            uint8_t value = 0;
            if (t < T)
            {
                value ^= measurement[v*T+t];
                for (auto q : lattice.incident_faces[v]) value ^= data[q*T+t];
            }
            if (t) value ^= measurement[v*T+t-1];
            require(value == syndrome[t*m+v], "History fails full syndrome validation");
        }
    return std::accumulate(data.begin(),data.end(),0.0)+std::accumulate(measurement.begin(),measurement.end(),0.0);
}

Outcome
Decoder::lift(const syndrome_type& syndrome, const State& state) const
{
    auto paths = assignments(syndrome,state);
    syndrome_type spatial(lattice.num_checks);
    for (size_t v = 0; v < num_detectors; v++) spatial[vertices[v].canonical_detector] ^= syndrome[v];
    Outcome best;
    best.weight = std::numeric_limits<double>::infinity();
    for (int c = 0; c < 3; c++)
    {
        try
        {
            auto components = build_connected_components(vertices,boundaries,paths,c);
            auto flattened = paths;
            for (auto& a : flattened) for (auto& [x,y] : a.path)
            {
                x = vertices[x].canonical_detector; y = vertices[y].canonical_detector;
            }
            auto result = lifter_.lifting(spatial,flattened,components,c);
            double weight = rounds == 1 ? result.correction.size() : score_history(syndrome,result);
            if (weight < best.weight)
            {
                best.observable = result.observable_mask;
                best.weight = weight;
            }
        }
        catch (const std::runtime_error& e)
        {
            ++best.lift_failures;
            best.lifting_errors.push_back(e.what());
        }
    }
    return best;
}

Outcome
Decoder::decode(const syndrome_type& syndrome)
{
    if (correlated_) return decode_single_pass(syndrome)[1];
    auto state = match(syndrome,nullptr,false);
    auto out = lift(syndrome,state);
    out.state = std::move(state);
    out.termination = "baseline";
    return out;
}

Outcome
Decoder::inspect_lift(const syndrome_type& syndrome, const State& state) const
{
    return lift(syndrome,state);
}

std::array<Outcome,2>
Decoder::decode_single_pass(const syndrome_type& syndrome)
{
    require(syndrome.size() == num_detectors, "Wrong syndrome size");
    const auto initial = match(syndrome,nullptr,false);
    const auto probabilities = conditional_probabilities(initial);
    const auto updated = match(syndrome,&probabilities,false);
    std::array<Outcome,2> out;
    out[0] = lift(syndrome,initial);
    out[0].state = initial;
    out[0].termination = "baseline";
    out[1] = updated == initial ? out[0] : lift(syndrome,updated);
    out[1].state = updated;
    out[1].termination = "one_pass";
    out[1].iterations = 1;
    return out;
}

} // namespace dec::correlated_restriction
