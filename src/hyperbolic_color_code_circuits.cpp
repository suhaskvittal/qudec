#include "circuit_generator.h"

#include <algorithm>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using Rows = std::vector<std::vector<uint32_t>>;

Rows read_checks(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open hyperbolic checks: " + path);
    Rows rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '%') continue;
        std::istringstream stream(line);
        std::vector<uint32_t> row;
        uint32_t q;
        while (stream >> q) {
            row.push_back(q);
        }
        if (row.empty()) continue;
        auto sorted = row;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            throw std::runtime_error("Repeated qubit in a hyperbolic face");
        rows.push_back(std::move(row));
    }
    if (rows.empty()) throw std::runtime_error("No check supports found in " + path);
    return rows;
}

Rows read_logicals(const std::string& path, uint32_t n) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open hyperbolic logical matrix: " + path);
    std::string line;
    if (!std::getline(in, line) || line.rfind("%%MatrixMarket matrix coordinate integer general", 0) != 0)
        throw std::runtime_error("Expected an integer coordinate MatrixMarket file: " + path);
    do {
        if (!std::getline(in, line)) throw std::runtime_error("Missing MatrixMarket dimensions: " + path);
    } while (line.empty() || line[0] == '%');
    uint32_t nrows = 0, ncols = 0, nnz = 0;
    std::istringstream(line) >> nrows >> ncols >> nnz;
    if (ncols != n) throw std::runtime_error("Logical matrix qubit count does not match H.txt");
    Rows rows(nrows);
    for (uint32_t i = 0; i < nnz; ++i) {
        uint32_t r, c;
        int value;
        if (!(in >> r >> c >> value) || r == 0 || r > nrows || c == 0 || c > n)
            throw std::runtime_error("Invalid logical matrix entry: " + path);
        if (value & 1) rows[r - 1].push_back(c - 1);
    }
    // MatrixMarket entries are additive; repeated coordinates cancel over GF(2).
    for (auto& row : rows) {
        std::sort(row.begin(), row.end());
        std::vector<uint32_t> reduced;
        for (size_t i = 0; i < row.size();) {
            size_t j = i + 1;
            while (j < row.size() && row[j] == row[i]) ++j;
            if ((j - i) & 1) reduced.push_back(row[i]);
            i = j;
        }
        row = std::move(reduced);
    }
    return rows;
}

std::vector<uint32_t> read_colors(const std::string& path, size_t faces) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open hyperbolic face colors: " + path);
    std::vector<uint32_t> colors;
    uint32_t color;
    while (in >> color) {
        if (color > 2) throw std::runtime_error("Invalid hyperbolic face color");
        colors.push_back(color);
    }
    if (colors.size() != faces) throw std::runtime_error("Face color count does not match H.txt");
    return colors;
}

// Recover each face's cyclic vertex order: in these trivalent tilings two face
// vertices are adjacent precisely when they also occur in another face. Each
// Bell ancilla contacts one contiguous half of this cycle in cyclic order.
// This keeps every propagated hook on a contiguous face arc.
Rows ordered_ancilla_paths(const Rows& checks, uint32_t n) {
    std::vector<std::vector<uint32_t>> qubit_faces(n);
    for (uint32_t f = 0; f < checks.size(); ++f)
        for (auto q : checks[f]) qubit_faces[q].push_back(f);
    Rows paths;
    paths.reserve(2 * checks.size());
    for (uint32_t f = 0; f < checks.size(); ++f) {
        const auto& row = checks[f];
        if (row.size() < 4 || (row.size() & 1))
            throw std::runtime_error("Hyperbolic faces need even weight at least four");
        std::vector<std::vector<uint32_t>> neighbors(row.size());
        for (size_t i = 0; i < row.size(); ++i)
            for (size_t j = i + 1; j < row.size(); ++j) {
                bool adjacent = false;
                for (auto other : qubit_faces[row[i]])
                    if (other != f && std::find(checks[other].begin(), checks[other].end(), row[j]) != checks[other].end())
                        adjacent = true;
                if (adjacent) {
                    neighbors[i].push_back(row[j]);
                    neighbors[j].push_back(row[i]);
                }
            }
        for (const auto& list : neighbors)
            if (list.size() != 2) throw std::runtime_error("Face supports do not define a cycle");
        const auto start = std::min_element(row.begin(), row.end()) - row.begin();
        std::vector<uint32_t> cycle{row[start]};
        cycle.push_back(std::min(neighbors[start][0], neighbors[start][1]));
        while (cycle.size() < row.size()) {
            const auto at = std::find(row.begin(), row.end(), cycle.back()) - row.begin();
            const auto next = neighbors[at][0] == cycle[cycle.size() - 2]
                ? neighbors[at][1] : neighbors[at][0];
            if (std::find(cycle.begin(), cycle.end(), next) != cycle.end())
                throw std::runtime_error("Face support cycle closes prematurely");
            cycle.push_back(next);
        }
        const auto last = std::find(row.begin(), row.end(), cycle.back()) - row.begin();
        if (std::find(neighbors[last].begin(), neighbors[last].end(), cycle.front()) == neighbors[last].end())
            throw std::runtime_error("Face support cycle does not close");
        const auto half = cycle.size() / 2;
        paths.emplace_back(cycle.begin(), cycle.begin() + half);
        paths.emplace_back(cycle.begin() + half, cycle.end());
    }
    return paths;
}

// Schedule the next contact of each ancilla on the earliest collision-free
// layer. Prefer longer unfinished paths when several ancillas need one qubit.
// X and Z bases run in separate passes, so no cross-basis CNOT conflict exists.
std::vector<std::vector<uint32_t>> cnot_layers(const Rows& checks, uint32_t n) {
    auto paths = ordered_ancilla_paths(checks, n);
    std::vector<size_t> next(paths.size());
    size_t remaining = 0;
    for (const auto& path : paths) remaining += path.size();
    std::vector<std::vector<uint32_t>> layers;
    while (remaining) {
        std::vector<int> chosen(n, -1);
        std::vector<uint32_t> layer;
        for (size_t step = 0; step < paths.size(); ++step) {
            const auto anc = (layers.size() + step) % paths.size();
            if (next[anc] == paths[anc].size()) continue;
            const auto q = paths[anc][next[anc]];
            const auto previous = chosen[q];
            if (previous < 0 ||
                paths[anc].size() - next[anc] > paths[previous].size() - next[previous])
                chosen[q] = int(anc);
        }
        for (uint32_t q = 0; q < n; ++q) {
            const auto anc = chosen[q];
            if (anc < 0) continue;
            layer.insert(layer.end(), {uint32_t(anc), q});
            ++next[anc];
            --remaining;
        }
        if (layer.empty()) throw std::logic_error("No progress in hyperbolic CNOT scheduling");
        layers.push_back(std::move(layer));
    }
    const auto nf = checks.size();
    for (const auto& layer : layers) {
        std::vector<bool> used_ancilla(2 * nf), used_data(n);
        for (size_t i = 0; i < layer.size(); i += 2) {
            if (used_ancilla[layer[i]] || used_data[layer[i + 1]])
                throw std::logic_error("Conflicting CNOTs in a hyperbolic extraction layer");
            used_ancilla[layer[i]] = true;
            used_data[layer[i + 1]] = true;
        }
    }
    return layers;
}
}  // namespace

stim::Circuit hyperbolic_color_code_circuit(const std::string& directory, uint32_t rounds,
    const CircuitNoisePolicy& policy, bool memory_x, bool include_opposite_basis_detectors,
    const ComponentErrorDistribution& noise) {
    if (!rounds) throw std::invalid_argument("Hyperbolic color code needs at least one round");
    policy.validate(rounds);
    const bool bell_flagged = policy.check_implementation() == CheckImplementation::GateSequence;
    const auto base = directory.empty() || directory.back() != '/' ? directory + "/" : directory;
    auto checks = read_checks(base + "H.txt");
    uint32_t n = 0;
    for (const auto& row : checks) for (auto q : row) n = std::max(n, q + 1);
    const uint32_t nf = checks.size();
    const uint32_t m = (bell_flagged ? 4 : 2) * nf;
    auto colors = read_colors(base + "check_colors.txt", nf);
    auto logicals = read_logicals(base + (memory_x ? "Lx.mtx" : "Lz.mtx"), n);
    const auto layers = bell_flagged ? cnot_layers(checks, n) : std::vector<std::vector<uint32_t>>{};

    std::vector<uint32_t> data(n), all(n + (bell_flagged ? 2 * nf : 0)), a(nf), b(nf), pairs, readout;
    std::iota(data.begin(), data.end(), 0);
    std::iota(all.begin(), all.end(), 0);
    for (uint32_t f = 0; f < nf; ++f) {
        a[f] = n + 2 * f;
        b[f] = a[f] + 1;
        pairs.insert(pairs.end(), {a[f], b[f]});
    }
    readout = a;
    readout.insert(readout.end(), b.begin(), b.end());

    stim::Circuit head;
    for (auto q : data) head.safe_append_u("QUBIT_COORDS", {q}, {(double)q, 0});
    for (uint32_t f = 0; f < nf; ++f) {
        if (bell_flagged) {
            head.safe_append_u("QUBIT_COORDS", {a[f]}, {(double)f, 1});
            head.safe_append_u("QUBIT_COORDS", {b[f]}, {(double)f, 2});
        }
    }
    auto idle = [&](stim::Circuit& c, const std::vector<uint32_t>& active, CircuitOperationKind kind) {
        std::vector<bool> busy(all.size());
        for (auto q : active) busy[q] = true;
        std::vector<uint32_t> inactive;
        for (auto q : all) if (!busy[q]) inactive.push_back(q);
        if (!inactive.empty()) policy.append(c, {kind, inactive});
    };
    auto h = [&](stim::Circuit& c, const std::vector<uint32_t>& qs) {
        policy.append(c, {CircuitOperationKind::Unitary1, qs, "H"});
        idle(c, qs, CircuitOperationKind::GateIdle);
        c.safe_append_u("TICK", {});
    };
    auto cx = [&](stim::Circuit& c, const std::vector<uint32_t>& qs) {
        policy.append(c, {CircuitOperationKind::Unitary2, qs, "CX"});
        idle(c, qs, CircuitOperationKind::GateIdle);
        c.safe_append_u("TICK", {});
    };
    policy.append(head, {CircuitOperationKind::Reset, all});
    if (bell_flagged) head.safe_append_u("TICK", {});
    if (memory_x) h(head, data);

    stim::Circuit cycle;
    policy.append(cycle, {CircuitOperationKind::RoundStart, data});
    std::vector<uint32_t> check_order(2 * nf), flag_order(2 * nf);
    if (!bell_flagged) {
        for (uint32_t basis = 0; basis < 2; ++basis) {
            for (uint32_t f = 0; f < nf; ++f) {
                std::vector<uint32_t> product;
                for (auto q : checks[f]) {
                    if (!product.empty()) product.push_back(stim::GateTarget::combiner().data);
                    product.push_back((basis == 0 ? stim::GateTarget::x(q) : stim::GateTarget::z(q)).data);
                }
                policy.append(cycle, {CircuitOperationKind::ProductMeasurement, product});
                check_order[basis * nf + f] = basis * nf + f;
            }
        }
    }
    if (bell_flagged) {
        for (uint32_t basis = 0; basis < 2; ++basis) {
            h(cycle, a);
            cx(cycle, pairs);
            for (const auto& layer : layers) {
                std::vector<uint32_t> gates;
                gates.reserve(layer.size());
                for (size_t i = 0; i < layer.size(); i += 2) {
                    const auto anc = n + layer[i], q = layer[i + 1];
                    if (basis == 0) gates.insert(gates.end(), {anc, q});
                    else gates.insert(gates.end(), {q, anc});
                }
                cx(cycle, gates);
            }
            cx(cycle, pairs);
            h(cycle, a);
            policy.append(cycle, {CircuitOperationKind::MeasureReset, readout});
            idle(cycle, readout, CircuitOperationKind::GateIdle);
            idle(cycle, readout, CircuitOperationKind::MeasurementIdle);
            cycle.safe_append_u("TICK", {});
            for (uint32_t f = 0; f < nf; ++f) {
                const auto offset = basis * 2 * nf;
                check_order[basis * nf + f] = offset + (basis == 0 ? f : nf + f);
                flag_order[basis * nf + f] = offset + (basis == 0 ? nf + f : f);
            }
        }
    }
    auto rec = [](uint32_t back) { return back | stim::TARGET_RECORD_BIT; };
    auto detectors = [&](stim::Circuit& c, bool first) {
        for (uint32_t basis = 0; basis < 2; ++basis) {
            for (uint32_t f = 0; f < nf; ++f) {
                const auto index = basis * nf + f;
                const auto check_back = m - check_order[index];
                if ((!first || (basis == 0) == memory_x) &&
                    ((basis == 0) == memory_x || include_opposite_basis_detectors)) {
                    std::vector<uint32_t> records{rec(check_back)};
                    if (!first) records.push_back(rec(check_back + m));
                    c.safe_append_u("DETECTOR", records, {(double)f, 0, 0, (double)colors[f], (double)basis, 0});
                }
                if (bell_flagged)
                    c.safe_append_u("DETECTOR", {rec(m - flag_order[index])},
                                    {(double)f, 0, 0, (double)colors[f], (double)basis, 1});
            }
        }
    };
    head += cycle;
    detectors(head, true);
    if (!bell_flagged) head.safe_append_u("TICK", {});
    stim::Circuit body;
    body.safe_append_u("SHIFT_COORDS", {}, {0, 0, 1});
    body += cycle;
    detectors(body, false);
    if (!bell_flagged) body.safe_append_u("TICK", {});
    head += body * (rounds - 1);

    if (memory_x) h(head, data);
    policy.append(head, {CircuitOperationKind::Measure, data});
    if (bell_flagged) {
        policy.append(head, {CircuitOperationKind::PostMeasurement, data});
        idle(head, data, CircuitOperationKind::GateIdle);
        idle(head, data, CircuitOperationKind::MeasurementIdle);
    }
    const auto basis = memory_x ? 0u : 1u;
    for (uint32_t f = 0; f < nf; ++f) {
        std::vector<uint32_t> records{rec(n + m - check_order[basis * nf + f])};
        for (auto q : checks[f]) records.push_back(rec(n - q));
        head.safe_append_u("DETECTOR", records, {(double)f, 0, 1, (double)colors[f], (double)basis, 0});
    }
    for (uint32_t k = 0; k < logicals.size(); ++k) {
        std::vector<uint32_t> records;
        for (auto q : logicals[k]) records.push_back(rec(n - q));
        head.safe_append_ua("OBSERVABLE_INCLUDE", records, k);
    }
    return apply_component_error_distribution(head, noise);
}
