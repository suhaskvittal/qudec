// Implementation authored by OpenAI Codex.
// Square/octagon layout portions: Copyright 2023 Google LLC, Apache-2.0.
// Adapted for bell-flagged extraction; upstream folding walks are replaced
// with face supports. License: docs/color_code_layout_LICENSE.
#include "circuit_generator.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using Point = std::pair<int, int>;
struct Face {
    int color;
    std::vector<Point> vertices;
};

// Honeycomb vertices and face centers on a triangular integer grid. Physical
// coordinates are (x/2, sqrt(3)*y/2); x already includes the row offset.
std::vector<Face> layout666(uint32_t distance) {
    const int width = distance + (distance - 1) / 2;
    std::set<Point> data;
    std::vector<Point> centers;
    for (int y = 0; y < width; ++y)
        for (int x = 0; x < width - y; ++x)
            if ((x + 2 * y) % 3 == 2)
                centers.emplace_back(2 * x + y, y);
            else
                data.emplace(2 * x + y, y);
    std::vector<Face> faces;
    for (auto [x, y] : centers) {
        Face f{(x + y) % 3, {}};
        for (auto [dx, dy] : std::vector<Point>{{2, 0}, {1, 1}, {-1, 1}, {-2, 0}, {-1, -1}, {1, -1}})
            if (data.contains({x + dx, y + dy}))
                f.vertices.emplace_back(x + dx, y + dy);
        faces.push_back(std::move(f));
    }
    return faces;
}

// Smooth triangular square/octagon patch. The integer rectangular embedding
// below is adapted from quantumlib/chromobius's _color_code_layouts.py
// (Copyright 2023 Google LLC, Apache-2.0; see docs/color_code_layout_LICENSE).
// Repeated vertices in its folding walks are omitted here: extraction couples
// once to each vertex, instead of implementing an inline folding circuit.
std::vector<Face> layout488(uint32_t distance) {
    const int w = distance, half = w / 2;
    std::vector<Face> faces;
    auto add = [&](int color, int x, int y, std::vector<Point> offsets) {
        Face f{color, {}};
        for (auto [dx, dy] : offsets)
            f.vertices.emplace_back(x + dx, y + dy);
        faces.push_back(std::move(f));
    };
    for (int x = 1; x < w; x += 2)
        add(2, x, 0, {{0, 1}, {0, 0}, {1, 0}, {1, 1}});
    for (int x = 0; x < half; x += 2)
        add(1, x, 2 * x, {{0, 0}, {1, 0}, {1, 1}, {1, 2}});
    for (int x = 0; x < half - 1; x += 2)
        add(2, w - x - 2, 2 * x + 2, {{0, 0}, {0, 1}, {0, 2}, {1, 0}});
    for (int offset = 0; offset < w; offset += 2) {
        for (int x = 0; x < half - offset / 2; ++x)
            add(0, x + 1 + offset, 2 * x + 1, {{0, 0}, {0, 1}, {1, 1}, {1, 0}});
        for (int x = 0; x < half - offset / 2 - 1; x += 2)
            add(1, x + 2 + offset, 2 * x, {{0, 0}, {0, 1}, {0, 2}, {0, 3},
                                                   {1, 0}, {1, 1}, {1, 2}, {1, 3}});
        for (int x = 0; x < half - offset / 2 - 2; x += 2)
            add(2, x + 3 + offset, 2 * x + 2, {{0, 0}, {0, 1}, {0, 2}, {0, 3},
                                                       {1, 0}, {1, 1}, {1, 2}, {1, 3}});
    }
    return faces;
}

std::pair<double, double> position(Point q, bool octagons) {
    auto [x, y] = q;
    if (!octagons)
        return {x * 0.5, y * std::sqrt(3.0) * 0.5};
    // Expand alternate pairs of rows into the square/octagon embedding.
    int r = 2 * x;
    if (((y + 1) / 2) % 2 == 0)
        r += (x + 1) % 2;
    else
        r += x % 2;
    return {r + (r % 2 == 0 ? 0.1 : -0.1), y + (y % 2 == 0 ? 0.1 : -0.1)};
}

stim::Circuit make_color_memory(uint32_t distance, uint32_t rounds, const CircuitNoisePolicy& policy,
                                bool memory_x, bool octagons, bool include_opposite_basis_detectors) {
    if (distance < 3 || distance % 2 == 0 || rounds == 0)
        throw std::invalid_argument("color code: need odd distance >= 3 and rounds >= 1");
    policy.validate(rounds);
    const bool bell_flagged = policy.check_implementation() == CheckImplementation::GateSequence;
    auto faces = octagons ? layout488(distance) : layout666(distance);
    std::set<Point> points;
    for (const auto& f : faces)
        points.insert(f.vertices.begin(), f.vertices.end());
    std::map<Point, uint32_t> index;
    stim::Circuit head;
    for (auto q : points) {
        uint32_t i = index.size();
        index[q] = i;
        auto [x, y] = position(q, octagons);
        head.safe_append_u("QUBIT_COORDS", {i}, {x, y});
    }
    const uint32_t n = points.size(), nf = faces.size(), m = (bell_flagged ? 4 : 2) * nf;
    std::vector<uint32_t> data(n), all(n + (bell_flagged ? 2 : 0) * nf), a, b;
    std::iota(data.begin(), data.end(), 0);
    std::iota(all.begin(), all.end(), 0);
    std::vector<std::vector<uint32_t>> supports;
    std::vector<std::pair<double, double>> centers;
    for (uint32_t k = 0; k < nf; ++k) {
        auto& f = faces[k];
        double x = 0, y = 0;
        for (auto q : f.vertices) {
            auto pos = position(q, octagons);
            x += pos.first / f.vertices.size();
            y += pos.second / f.vertices.size();
        }
        centers.emplace_back(x, y);
        if (bell_flagged) std::sort(f.vertices.begin(), f.vertices.end(), [&](Point u, Point v) {
            auto pu = position(u, octagons), pv = position(v, octagons);
            return std::atan2(pu.second - y, pu.first - x) < std::atan2(pv.second - y, pv.first - x);
        });
        std::vector<uint32_t> support;
        for (auto q : f.vertices) support.push_back(index.at(q));
        supports.push_back(std::move(support));
        a.push_back(n + 2 * k);
        b.push_back(n + 2 * k + 1);
        // Fourth qubit coordinate identifies data/ancilla role: 0/1/2.
        if (bell_flagged) {
            head.safe_append_u("QUBIT_COORDS", {a.back()}, {x - 0.15, y, (double)f.color, 1});
            head.safe_append_u("QUBIT_COORDS", {b.back()}, {x + 0.15, y, (double)f.color, 2});
        }
    }
    auto idle = [&](stim::Circuit& c, const std::vector<uint32_t>& active, CircuitOperationKind kind) {
        std::set<uint32_t> busy(active.begin(), active.end());
        std::vector<uint32_t> qs;
        for (auto q : all) if (!busy.contains(q)) qs.push_back(q);
        if (!qs.empty()) policy.append(c, {kind, qs});
    };
    auto cx = [&](stim::Circuit& c, const std::vector<uint32_t>& pairs) {
        if (!pairs.empty()) policy.append(c, {CircuitOperationKind::Unitary2, pairs, "CX"});
        idle(c, pairs, CircuitOperationKind::GateIdle);
        c.safe_append_u("TICK", {});
    };
    auto h = [&](stim::Circuit& c, const std::vector<uint32_t>& qs) {
        policy.append(c, {CircuitOperationKind::Unitary1, qs, "H"});
        if (bell_flagged) {
            idle(c, qs, CircuitOperationKind::GateIdle);
            c.safe_append_u("TICK", {});
        }
    };
    auto slow_idle = [&](stim::Circuit& c, const std::vector<uint32_t>& active) {
        idle(c, active, CircuitOperationKind::GateIdle);
        idle(c, active, CircuitOperationKind::MeasurementIdle);
    };
    policy.append(head, {CircuitOperationKind::Reset, data});
    if (bell_flagged) {
        slow_idle(head, data);
        head.safe_append_u("TICK", {});
    }
    if (memory_x) h(head, data);

    stim::Circuit cycle;
    // Fixed measurement order: basis X/Z, color R/G/B, all a then all b.
    std::vector<uint32_t> check_order(2 * nf), flag_order(2 * nf);
    uint32_t order = 0;
    policy.append(cycle, {CircuitOperationKind::RoundStart, data});
    if (!bell_flagged) {
        for (int basis = 0; basis < 2; ++basis) for (uint32_t k = 0; k < nf; ++k) {
            std::vector<uint32_t> product;
            for (auto q : supports[k]) {
                if (!product.empty()) product.push_back(stim::GateTarget::combiner().data);
                product.push_back((basis == 0 ? stim::GateTarget::x(q) : stim::GateTarget::z(q)).data);
            }
            policy.append(cycle, {CircuitOperationKind::ProductMeasurement, product});
            check_order[basis * nf + k] = order++;
        }
    }
    if (bell_flagged) for (int basis = 0; basis < 2; ++basis) {
        for (int color = 0; color < 3; ++color) {
            std::vector<uint32_t> aa, bb, active, bell;
            size_t half_weight = 0;
            for (uint32_t k = 0; k < nf; ++k) if (faces[k].color == color) {
                aa.push_back(a[k]); bb.push_back(b[k]);
                active.insert(active.end(), {a[k], b[k]});
                bell.insert(bell.end(), {a[k], b[k]});
                half_weight = std::max(half_weight, supports[k].size() / 2);
            }
            if (aa.empty()) continue;
            policy.append(cycle, {CircuitOperationKind::Reset, active});
            slow_idle(cycle, active);
            cycle.safe_append_u("TICK", {});
            h(cycle, aa);
            cx(cycle, bell);
            for (size_t step = 0; step < half_weight; ++step) {
                std::vector<uint32_t> pairs;
                for (uint32_t k = 0; k < nf; ++k) if (faces[k].color == color) {
                    const auto& s = supports[k];
                    if (step >= s.size() / 2) continue;
                    for (int side = 0; side < 2; ++side) {
                        uint32_t anc = side == 0 ? a[k] : b[k];
                        uint32_t q = s[step + side * (s.size() / 2)];
                        if (basis == 0) pairs.insert(pairs.end(), {anc, q});
                        else pairs.insert(pairs.end(), {q, anc});
                    }
                }
                cx(cycle, pairs);
            }
            cx(cycle, bell);
            h(cycle, aa);
            // Preserve a single readout-flip channel on all active ancillas,
            // followed by the historical a-then-b measurement record ordering.
            std::vector<uint32_t> readout = aa;
            readout.insert(readout.end(), bb.begin(), bb.end());
            policy.append(cycle, {CircuitOperationKind::Measure, readout});
            policy.append(cycle, {CircuitOperationKind::PostMeasurement, active});
            slow_idle(cycle, active);
            cycle.safe_append_u("TICK", {});
            for (uint32_t j = 0; j < aa.size(); ++j) {
                uint32_t k = (aa[j] - n) / 2;
                check_order[basis * nf + k] = order + j + (basis == 0 ? 0 : aa.size());
                flag_order[basis * nf + k] = order + j + (basis == 0 ? aa.size() : 0);
            }
            order += active.size();
        }
    }
    auto rec = [](uint32_t back) { return back | stim::TARGET_RECORD_BIT; };
    auto detectors = [&](stim::Circuit& c, bool first) {
        for (int basis = 0; basis < 2; ++basis)
            for (uint32_t k = 0; k < nf; ++k) {
                auto [x, y] = centers[k];
                auto back = m - check_order[basis * nf + k];
                if ((!first || (basis == 0) == memory_x) &&
                    ((basis == 0) == memory_x || include_opposite_basis_detectors)) {
                    std::vector<uint32_t> records{rec(back)};
                    if (!first) records.push_back(rec(back + m));
                    c.safe_append_u("DETECTOR", records, {x, y, 0, (double)faces[k].color, (double)basis, 0});
                }
                if (bell_flagged) c.safe_append_u("DETECTOR", {rec(m - flag_order[basis * nf + k])},
                                {x, y, 0, (double)faces[k].color, (double)basis, 1});
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
    if (bell_flagged && memory_x) h(head, data);
    policy.append(head, {CircuitOperationKind::Measure, data, {}, (!bell_flagged && memory_x) ? 'X' : 'Z'});
    if (bell_flagged) {
        policy.append(head, {CircuitOperationKind::PostMeasurement, data});
        slow_idle(head, data);
    }
    int basis = memory_x ? 0 : 1;
    for (uint32_t k = 0; k < nf; ++k) {
        std::vector<uint32_t> records{rec(n + m - check_order[basis * nf + k])};
        for (auto q : supports[k]) records.push_back(rec(n - q));
        auto [x, y] = centers[k];
        head.safe_append_u("DETECTOR", records, {x, y, 1, (double)faces[k].color, (double)basis, 0});
    }
    std::vector<uint32_t> observable;
    for (auto q : data) observable.push_back(rec(n - q));
    head.safe_append_ua("OBSERVABLE_INCLUDE", observable, 0);
    return head;
}
} // namespace

// Author: OpenAI Codex. One layout/round/detector generator per color code.
stim::Circuit cc666_circuit(uint32_t d, uint32_t r, const CircuitNoisePolicy& policy,
    bool x, bool opposite, const ComponentErrorDistribution& noise) {
    return apply_component_error_distribution(make_color_memory(d,r,policy,x,false,opposite),noise);
}
stim::Circuit cc488_circuit(uint32_t d, uint32_t r, const CircuitNoisePolicy& policy,
    bool x, bool opposite, const ComponentErrorDistribution& noise) {
    return apply_component_error_distribution(make_color_memory(d,r,policy,x,true,opposite),noise);
}
// Legacy model-named entry points are compatibility adapters only.
stim::Circuit cc666_si1000(uint32_t d,uint32_t r,double p,bool x,bool opposite,
    const ComponentErrorDistribution& noise) {
    return cc666_circuit(d,r,SI1000NoisePolicy(p),x,opposite,noise);
}
stim::Circuit cc488_si1000(uint32_t d,uint32_t r,double p,bool x,bool opposite,
    const ComponentErrorDistribution& noise) {
    return cc488_circuit(d,r,SI1000NoisePolicy(p),x,opposite,noise);
}
stim::Circuit cc666_code_capacity(uint32_t d,double p,bool x,bool opposite,
    const ComponentErrorDistribution& noise) {
    return cc666_circuit(d,1,CodeCapacityNoisePolicy(p),x,opposite,noise);
}
stim::Circuit cc488_code_capacity(uint32_t d,double p,bool x,bool opposite,
    const ComponentErrorDistribution& noise) {
    return cc488_circuit(d,1,CodeCapacityNoisePolicy(p),x,opposite,noise);
}
stim::Circuit cc666_phenomenological(uint32_t d,uint32_t r,double p,double q,bool x,
    bool opposite,const ComponentErrorDistribution& noise) {
    return cc666_circuit(d,r,PhenomenologicalNoisePolicy(p,q),x,opposite,noise);
}
stim::Circuit cc488_phenomenological(uint32_t d,uint32_t r,double p,double q,bool x,
    bool opposite,const ComponentErrorDistribution& noise) {
    return cc488_circuit(d,r,PhenomenologicalNoisePolicy(p,q),x,opposite,noise);
}
