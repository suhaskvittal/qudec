/*
 *  author: Claude Sonnet 5.1
 *  date:   26 July 2026
 * */

#include "circuit_generator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

/*
 * A 2D coordinate on the surface-code lattice. Mirrors the helper used by
 * stim's rotated surface-code generator so the geometry and CX schedule match.
 * */
struct surface_coord
{
    float x;
    float y;

    surface_coord operator+(surface_coord o) const { return {x + o.x, y + o.y}; }
    surface_coord operator-(surface_coord o) const { return {x - o.x, y - o.y}; }
    bool operator==(surface_coord o) const { return x == o.x && y == o.y; }
    bool operator<(surface_coord o) const { return x != o.x ? x < o.x : y < o.y; }
};

/*
 * Returns the qubits in `all` that do not appear in `active` (both sorted).
 * These are the idle qubits during a gate layer.
 * */
std::vector<uint32_t>
idle_qubits(const std::vector<uint32_t>& all, const std::vector<uint32_t>& active)
{
    std::vector<uint32_t> active_sorted(active);
    std::sort(active_sorted.begin(), active_sorted.end());
    std::vector<uint32_t> out;
    std::set_difference(all.begin(), all.end(),
                        active_sorted.begin(), active_sorted.end(),
                        std::back_inserter(out));
    return out;
}

/*
 * The SI1000 noise model, holding the per-location error rates derived from `p`
 * and providing the noisy-gate appenders used to build the circuit.
 * */
struct si1000_model
{
    double p;

    // Anti-basis flip: for a Z-basis (de)coherent op an X error is relevant,
    // and vice-versa. `basis` is the measurement/reset basis ('Z' or 'X').
    void
    anti_basis_flip(stim::Circuit& c, const std::vector<uint32_t>& targets, double pr, char basis) const
    {
        if (pr <= 0 || targets.empty())
            return;
        c.safe_append_ua(basis == 'X' ? "Z_ERROR" : "X_ERROR", targets, pr);
    }

    // Single-qubit gate: gate, then DEPOLARIZE1(p/10).
    void
    append_unitary_1(stim::Circuit& c, std::string_view name, const std::vector<uint32_t>& targets) const
    {
        if (targets.empty())
            return;
        c.safe_append_u(name, targets);
        c.safe_append_ua("DEPOLARIZE1", targets, p / 10.0);
    }

    // Two-qubit CX layer: gate, DEPOLARIZE2(p) on the pairs, then idle
    // DEPOLARIZE1(p) on every qubit untouched by the layer.
    void
    append_cx_layer(stim::Circuit& c,
                    const std::vector<uint32_t>& targets,
                    const std::vector<uint32_t>& all_qubits) const
    {
        c.safe_append_u("CX", targets);
        c.safe_append_ua("DEPOLARIZE2", targets, p);

        auto idle = idle_qubits(all_qubits, targets);
        if (!idle.empty())
            c.safe_append_ua("DEPOLARIZE1", idle, p);
    }

    // Reset (always Z-basis for ancilla; data uses the memory basis): gate,
    // then X_ERROR(2p) (anti-basis flip).
    void
    append_reset(stim::Circuit& c, const std::vector<uint32_t>& targets, char basis) const
    {
        std::string gate("R");
        gate.push_back(basis);
        c.safe_append_u(gate, targets);
        anti_basis_flip(c, targets, 2.0 * p, basis);
    }

    // Combined measure+reset for the ancilla (Z basis): X_ERROR(5p) before,
    // MR, then X_ERROR(2p) after (reset flip). No post-measure depolarize.
    // Data qubits are idle during this layer -> DEPOLARIZE1(2p).
    void
    append_measure_reset(stim::Circuit& c,
                         const std::vector<uint32_t>& targets,
                         const std::vector<uint32_t>& data_qubits) const
    {
        anti_basis_flip(c, targets, 5.0 * p, 'Z');
        c.safe_append_u("MR", targets);
        anti_basis_flip(c, targets, 2.0 * p, 'Z');
        if (!data_qubits.empty())
            c.safe_append_ua("DEPOLARIZE1", data_qubits, 4.0 * p);
    }

    // Final destructive data measurement in the memory basis: anti-basis
    // flip(5p) before, M, then DEPOLARIZE1(p) after.
    void
    append_measure(stim::Circuit& c, const std::vector<uint32_t>& targets, char basis) const
    {
        anti_basis_flip(c, targets, 5.0 * p, basis);
        std::string gate("M");
        gate.push_back(basis);
        c.safe_append_u(gate, targets);
        c.safe_append_ua("DEPOLARIZE1", targets, p);
    }
};

/*
 * A 2D coordinate on the LxL toric lattice, with both components taken
 * mod `2*L` (each axis is periodic). Mirrors `surface_coord` above but wraps
 * on addition so every stabilizer/qubit is well-defined with no boundary
 * special cases.
 * */
struct toric_coord
{
    int x;
    int y;
    int mod; // = 2*L

    static int
    wrap(int v, int m) { return ((v % m) + m) % m; }

    toric_coord
    operator+(std::pair<int, int> delta) const
    {
        return {wrap(x + delta.first, mod), wrap(y + delta.second, mod), mod};
    }
    bool operator==(toric_coord o) const { return x == o.x && y == o.y; }
    bool operator<(toric_coord o) const { return x != o.x ? x < o.x : y < o.y; }
};

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

stim::Circuit
sc_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x)
{
    using namespace stim;

    if (rounds < 1)
        throw std::invalid_argument("sc_si1000: need rounds >= 1.");
    if (distance < 2)
        throw std::invalid_argument("sc_si1000: need distance >= 2.");
    if (p < 0 || p > 1)
        throw std::invalid_argument("sc_si1000: need 0 <= p <= 1.");

    const uint32_t d = distance;
    const si1000_model model{p};

    // Place data qubits and collect the logical observable supports.
    std::set<surface_coord> data_coords;
    std::vector<surface_coord> x_observable;
    std::vector<surface_coord> z_observable;
    for (float x = 0.5; x <= d; x++)
    {
        for (float y = 0.5; y <= d; y++)
        {
            surface_coord q{x * 2, y * 2};
            data_coords.insert(q);
            if (y == 0.5)
                z_observable.push_back(q);
            if (x == 0.5)
                x_observable.push_back(q);
        }
    }

    // Place X/Z measurement qubits.
    std::set<surface_coord> x_measure_coords;
    std::set<surface_coord> z_measure_coords;
    for (size_t x = 0; x <= d; x++)
    {
        for (size_t y = 0; y <= d; y++)
        {
            surface_coord q{(float)x * 2, (float)y * 2};
            bool on_boundary_1 = x == 0 || x == d;
            bool on_boundary_2 = y == 0 || y == d;
            bool parity = x % 2 != y % 2;
            if (on_boundary_1 && parity)
                continue;
            if (on_boundary_2 && !parity)
                continue;
            if (parity)
                x_measure_coords.insert(q);
            else
                z_measure_coords.insert(q);
        }
    }

    // Interaction orders so hook errors run against the error grain.
    std::vector<surface_coord> z_order{{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
    std::vector<surface_coord> x_order{{1, 1}, {-1, 1}, {1, -1}, {-1, -1}};

    auto coord_to_index = [&](surface_coord q) -> uint32_t
    {
        q = q - surface_coord{0, fmodf(q.x, 2)};
        return (uint32_t)(q.x + q.y * (d + 0.5));
    };

    // Forward/reverse index every qubit.
    std::map<surface_coord, uint32_t> p2q;
    for (auto q : data_coords)
        p2q[q] = coord_to_index(q);
    for (auto q : x_measure_coords)
        p2q[q] = coord_to_index(q);
    for (auto q : z_measure_coords)
        p2q[q] = coord_to_index(q);
    std::map<uint32_t, surface_coord> q2p;
    for (const auto& kv : p2q)
        q2p[kv.second] = kv.first;

    // Target lists for the various qubit kinds.
    std::vector<uint32_t> data_qubits;
    std::vector<uint32_t> measurement_qubits;
    std::vector<uint32_t> x_measurement_qubits;
    std::vector<uint32_t> all_qubits;
    for (auto q : data_coords)
        data_qubits.push_back(p2q[q]);
    for (auto q : x_measure_coords)
    {
        measurement_qubits.push_back(p2q[q]);
        x_measurement_qubits.push_back(p2q[q]);
    }
    for (auto q : z_measure_coords)
        measurement_qubits.push_back(p2q[q]);
    all_qubits.insert(all_qubits.end(), data_qubits.begin(), data_qubits.end());
    all_qubits.insert(all_qubits.end(), measurement_qubits.begin(), measurement_qubits.end());
    std::sort(all_qubits.begin(), all_qubits.end());
    std::sort(data_qubits.begin(), data_qubits.end());
    std::sort(measurement_qubits.begin(), measurement_qubits.end());
    std::sort(x_measurement_qubits.begin(), x_measurement_qubits.end());

    // Measurement order used to compute detector record offsets.
    std::map<surface_coord, uint32_t> data_coord_to_order;
    std::map<surface_coord, uint32_t> measure_coord_to_order;
    for (auto q : data_qubits)
    {
        auto i = data_coord_to_order.size();
        data_coord_to_order[q2p[q]] = i;
    }
    for (auto q : measurement_qubits)
    {
        auto i = measure_coord_to_order.size();
        measure_coord_to_order[q2p[q]] = i;
    }

    // CX gate targets for each of the 4 sub-rounds.
    std::array<std::vector<uint32_t>, 4> cnot_targets;
    for (size_t k = 0; k < 4; k++)
    {
        for (auto measure : x_measure_coords)
        {
            auto data = measure + x_order[k];
            if (p2q.find(data) != p2q.end())
            {
                cnot_targets[k].push_back(p2q[measure]);
                cnot_targets[k].push_back(p2q[data]);
            }
        }
        for (auto measure : z_measure_coords)
        {
            auto data = measure + z_order[k];
            if (p2q.find(data) != p2q.end())
            {
                cnot_targets[k].push_back(p2q[data]);
                cnot_targets[k].push_back(p2q[measure]);
            }
        }
    }

    const auto& chosen_basis_observable = is_memory_x ? x_observable : z_observable;
    const auto& chosen_basis_measure_coords = is_memory_x ? x_measure_coords : z_measure_coords;

    // Repeated syndrome-extraction cycle.
    Circuit cycle_actions;
    cycle_actions.safe_append_u("TICK", {});
    model.append_unitary_1(cycle_actions, "H", x_measurement_qubits);
    for (const auto& targets : cnot_targets)
    {
        cycle_actions.safe_append_u("TICK", {});
        model.append_cx_layer(cycle_actions, targets, all_qubits);
    }
    cycle_actions.safe_append_u("TICK", {});
    model.append_unitary_1(cycle_actions, "H", x_measurement_qubits);
    cycle_actions.safe_append_u("TICK", {});
    model.append_measure_reset(cycle_actions, measurement_qubits, data_qubits);

    // Head: coords, reset, first cycle, and the first-round detectors.
    Circuit head;
    for (const auto& kv : q2p)
        head.safe_append_u("QUBIT_COORDS", {kv.first}, {kv.second.x, kv.second.y});
    model.append_reset(head, data_qubits, is_memory_x ? 'X' : 'Z');
    model.append_reset(head, measurement_qubits, 'Z');
    head += cycle_actions;
    for (auto measure : chosen_basis_measure_coords)
    {
        head.safe_append_u(
            "DETECTOR",
            {(uint32_t)(measurement_qubits.size() - measure_coord_to_order[measure]) | TARGET_RECORD_BIT},
            {measure.x, measure.y, 0});
    }

    // Body: cycle + detectors comparing this round to the previous one.
    Circuit body = cycle_actions;
    uint32_t m = measurement_qubits.size();
    body.safe_append_u("SHIFT_COORDS", {}, {0, 0, 1});
    for (auto m_coord : chosen_basis_measure_coords)
    {
        auto k = (uint32_t)measurement_qubits.size() - measure_coord_to_order[m_coord] - 1;
        body.safe_append_u(
            "DETECTOR", {(k + 1) | TARGET_RECORD_BIT, (k + 1 + m) | TARGET_RECORD_BIT}, {m_coord.x, m_coord.y, 0});
    }

    // Tail: destructive data readout, final detectors, and the observable.
    Circuit tail;
    model.append_measure(tail, data_qubits, is_memory_x ? 'X' : 'Z');
    for (auto measure : chosen_basis_measure_coords)
    {
        std::vector<uint32_t> detectors;
        for (auto delta : (is_memory_x ? x_order : z_order))
        {
            auto data = measure + delta;
            if (p2q.find(data) != p2q.end())
                detectors.push_back((data_qubits.size() - data_coord_to_order[data]) | TARGET_RECORD_BIT);
        }
        detectors.push_back(
            (data_qubits.size() + measurement_qubits.size() - measure_coord_to_order[measure]) | TARGET_RECORD_BIT);
        std::sort(detectors.begin(), detectors.end());
        tail.safe_append_u("DETECTOR", detectors, {measure.x, measure.y, 1});
    }
    std::vector<uint32_t> obs_inc;
    for (auto q : chosen_basis_observable)
        obs_inc.push_back((data_qubits.size() - data_coord_to_order[q]) | TARGET_RECORD_BIT);
    std::sort(obs_inc.begin(), obs_inc.end());
    tail.safe_append_ua("OBSERVABLE_INCLUDE", obs_inc, 0);

    return head + body * (rounds - 1) + tail;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

stim::Circuit
toric_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x)
{
    using namespace stim;

    if (rounds < 1)
        throw std::invalid_argument("toric_si1000: need rounds >= 1.");
    if (distance < 2)
        throw std::invalid_argument("toric_si1000: need distance >= 2.");
    if (p < 0 || p > 1)
        throw std::invalid_argument("toric_si1000: need 0 <= p <= 1.");

    const uint32_t L = distance;
    const int mod = 2 * (int)L;
    const si1000_model model{p};

    // Place data qubits (odd/even mixed parity) and ancilla qubits: X-check
    // vertices at (even, even), Z-check plaquettes at (odd, odd).
    std::set<toric_coord> data_coords;
    std::set<toric_coord> x_measure_coords;
    std::set<toric_coord> z_measure_coords;
    for (int x = 0; x < mod; x++)
    {
        for (int y = 0; y < mod; y++)
        {
            toric_coord q{x, y, mod};
            bool x_even = (x % 2) == 0;
            bool y_even = (y % 2) == 0;
            if (x_even && y_even)
                x_measure_coords.insert(q);
            else if (!x_even && !y_even)
                z_measure_coords.insert(q);
            else
                data_coords.insert(q);
        }
    }

    // Logical observable supports, per memory basis (two independent
    // observables on the torus).
    std::vector<toric_coord> observable_0, observable_1;
    if (is_memory_x)
    {
        for (int y = 0; y < mod; y += 2)
            observable_0.push_back(toric_coord{1, y, mod});
        for (int x = 0; x < mod; x += 2)
            observable_1.push_back(toric_coord{x, 1, mod});
    }
    else
    {
        for (int x = 1; x < mod; x += 2)
            observable_0.push_back(toric_coord{x, 0, mod});
        for (int y = 1; y < mod; y += 2)
            observable_1.push_back(toric_coord{0, y, mod});
    }

    // Interaction orders so hook errors run against the error grain.
    std::vector<std::pair<int, int>> x_order{{0, 1}, {-1, 0}, {1, 0}, {0, -1}};
    std::vector<std::pair<int, int>> z_order{{0, 1}, {1, 0}, {-1, 0}, {0, -1}};

    // Forward/reverse index every qubit (index = y * 2L + x).
    std::map<toric_coord, uint32_t> p2q;
    auto index_of = [&](toric_coord q) -> uint32_t { return (uint32_t)(q.y * mod + q.x); };
    for (auto q : data_coords)
        p2q[q] = index_of(q);
    for (auto q : x_measure_coords)
        p2q[q] = index_of(q);
    for (auto q : z_measure_coords)
        p2q[q] = index_of(q);
    std::map<uint32_t, toric_coord> q2p;
    for (const auto& kv : p2q)
        q2p[kv.second] = kv.first;

    // Target lists for the various qubit kinds.
    std::vector<uint32_t> data_qubits;
    std::vector<uint32_t> measurement_qubits;
    std::vector<uint32_t> x_measurement_qubits;
    std::vector<uint32_t> all_qubits;
    for (auto q : data_coords)
        data_qubits.push_back(p2q[q]);
    for (auto q : x_measure_coords)
    {
        measurement_qubits.push_back(p2q[q]);
        x_measurement_qubits.push_back(p2q[q]);
    }
    for (auto q : z_measure_coords)
        measurement_qubits.push_back(p2q[q]);
    all_qubits.insert(all_qubits.end(), data_qubits.begin(), data_qubits.end());
    all_qubits.insert(all_qubits.end(), measurement_qubits.begin(), measurement_qubits.end());
    std::sort(all_qubits.begin(), all_qubits.end());
    std::sort(data_qubits.begin(), data_qubits.end());
    std::sort(measurement_qubits.begin(), measurement_qubits.end());
    std::sort(x_measurement_qubits.begin(), x_measurement_qubits.end());

    // Measurement order used to compute detector record offsets.
    std::map<toric_coord, uint32_t> data_coord_to_order;
    std::map<toric_coord, uint32_t> measure_coord_to_order;
    for (auto q : data_qubits)
    {
        auto i = data_coord_to_order.size();
        data_coord_to_order[q2p[q]] = i;
    }
    for (auto q : measurement_qubits)
    {
        auto i = measure_coord_to_order.size();
        measure_coord_to_order[q2p[q]] = i;
    }

    // CX gate targets for each of the 4 sub-rounds. Every stabilizer is
    // weight-4 with periodic boundaries, so no membership guard is needed.
    std::array<std::vector<uint32_t>, 4> cnot_targets;
    for (size_t k = 0; k < 4; k++)
    {
        for (auto measure : x_measure_coords)
        {
            auto data = measure + x_order[k];
            cnot_targets[k].push_back(p2q[measure]);
            cnot_targets[k].push_back(p2q[data]);
        }
        for (auto measure : z_measure_coords)
        {
            auto data = measure + z_order[k];
            cnot_targets[k].push_back(p2q[data]);
            cnot_targets[k].push_back(p2q[measure]);
        }
    }

    const auto& chosen_basis_measure_coords = is_memory_x ? x_measure_coords : z_measure_coords;
    const auto& chosen_order = is_memory_x ? x_order : z_order;

    // Repeated syndrome-extraction cycle.
    Circuit cycle_actions;
    cycle_actions.safe_append_u("TICK", {});
    model.append_unitary_1(cycle_actions, "H", x_measurement_qubits);
    for (const auto& targets : cnot_targets)
    {
        cycle_actions.safe_append_u("TICK", {});
        model.append_cx_layer(cycle_actions, targets, all_qubits);
    }
    cycle_actions.safe_append_u("TICK", {});
    model.append_unitary_1(cycle_actions, "H", x_measurement_qubits);
    cycle_actions.safe_append_u("TICK", {});
    model.append_measure_reset(cycle_actions, measurement_qubits, data_qubits);

    // Head: coords, reset, first cycle, and the first-round detectors.
    Circuit head;
    for (const auto& kv : q2p)
        head.safe_append_u("QUBIT_COORDS", {kv.first}, {(float)kv.second.x, (float)kv.second.y});
    model.append_reset(head, data_qubits, is_memory_x ? 'X' : 'Z');
    model.append_reset(head, measurement_qubits, 'Z');
    head += cycle_actions;
    for (auto measure : chosen_basis_measure_coords)
    {
        head.safe_append_u(
            "DETECTOR",
            {(uint32_t)(measurement_qubits.size() - measure_coord_to_order[measure]) | TARGET_RECORD_BIT},
            {(float)measure.x, (float)measure.y, 0});
    }

    // Body: cycle + detectors comparing this round to the previous one.
    Circuit body = cycle_actions;
    uint32_t m = measurement_qubits.size();
    body.safe_append_u("SHIFT_COORDS", {}, {0, 0, 1});
    for (auto m_coord : chosen_basis_measure_coords)
    {
        auto k = (uint32_t)measurement_qubits.size() - measure_coord_to_order[m_coord] - 1;
        body.safe_append_u(
            "DETECTOR", {(k + 1) | TARGET_RECORD_BIT, (k + 1 + m) | TARGET_RECORD_BIT}, {(float)m_coord.x, (float)m_coord.y, 0});
    }

    // Tail: destructive data readout, final detectors, and both observables.
    Circuit tail;
    model.append_measure(tail, data_qubits, is_memory_x ? 'X' : 'Z');
    for (auto measure : chosen_basis_measure_coords)
    {
        std::vector<uint32_t> detectors;
        for (auto delta : chosen_order)
        {
            auto data = measure + delta;
            detectors.push_back((data_qubits.size() - data_coord_to_order[data]) | TARGET_RECORD_BIT);
        }
        detectors.push_back(
            (data_qubits.size() + measurement_qubits.size() - measure_coord_to_order[measure]) | TARGET_RECORD_BIT);
        std::sort(detectors.begin(), detectors.end());
        tail.safe_append_u("DETECTOR", detectors, {(float)measure.x, (float)measure.y, 1});
    }
    for (size_t obs_id = 0; obs_id < 2; obs_id++)
    {
        std::vector<uint32_t> obs_inc;
        for (auto q : (obs_id == 0 ? observable_0 : observable_1))
            obs_inc.push_back((data_qubits.size() - data_coord_to_order[q]) | TARGET_RECORD_BIT);
        std::sort(obs_inc.begin(), obs_inc.end());
        tail.safe_append_ua("OBSERVABLE_INCLUDE", obs_inc, (double)obs_id);
    }

    return head + body * (rounds - 1) + tail;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
