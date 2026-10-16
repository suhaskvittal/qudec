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

stim::Circuit make_color_memory(uint32_t distance, uint32_t rounds, double p,
                                bool memory_x, bool octagons) {
    if (distance < 3 || distance % 2 == 0 || rounds == 0)
        throw std::invalid_argument("color code: need odd distance >= 3 and rounds >= 1");
    if (!std::isfinite(p) || p < 0 || p > 0.2)
        throw std::invalid_argument("color code SI1000: need finite 0 <= p <= 0.2 (measurement error is 5p)");
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
    const uint32_t n = points.size(), nf = faces.size(), m = 4 * nf;
    std::vector<uint32_t> data(n), all(n + 2 * nf), a, b;
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
        std::sort(f.vertices.begin(), f.vertices.end(), [&](Point u, Point v) {
            auto pu = position(u, octagons), pv = position(v, octagons);
            return std::atan2(pu.second - y, pu.first - x) < std::atan2(pv.second - y, pv.first - x);
        });
        std::vector<uint32_t> support;
        for (auto q : f.vertices) support.push_back(index.at(q));
        supports.push_back(std::move(support));
        a.push_back(n + 2 * k);
        b.push_back(n + 2 * k + 1);
        // Fourth qubit coordinate identifies data/ancilla role: 0/1/2.
        head.safe_append_u("QUBIT_COORDS", {a.back()}, {x - 0.15, y, (double)f.color, 1});
        head.safe_append_u("QUBIT_COORDS", {b.back()}, {x + 0.15, y, (double)f.color, 2});
    }
    auto noise = [&](stim::Circuit& c, const char* gate, const std::vector<uint32_t>& qs, double rate) {
        if (rate > 0 && !qs.empty()) c.safe_append_ua(gate, qs, rate);
    };
    auto flip = [&](stim::Circuit& c, const std::vector<uint32_t>& qs, bool x_basis, double rate) {
        noise(c, x_basis ? "Z_ERROR" : "X_ERROR", qs, rate);
    };
    auto idle = [&](stim::Circuit& c, const std::vector<uint32_t>& active, double rate) {
        std::set<uint32_t> busy(active.begin(), active.end());
        std::vector<uint32_t> qs;
        for (auto q : all) if (!busy.contains(q)) qs.push_back(q);
        noise(c, "DEPOLARIZE1", qs, rate);
    };
    auto cx = [&](stim::Circuit& c, const std::vector<uint32_t>& pairs) {
        if (!pairs.empty()) c.safe_append_u("CX", pairs);
        noise(c, "DEPOLARIZE2", pairs, p);
        idle(c, pairs, p / 10);
        c.safe_append_u("TICK", {});
    };
    auto h = [&](stim::Circuit& c, const std::vector<uint32_t>& qs) {
        c.safe_append_u("H", qs);
        noise(c, "DEPOLARIZE1", qs, p / 10);
        idle(c, qs, p / 10);
        c.safe_append_u("TICK", {});
    };
    auto slow_idle = [&](stim::Circuit& c, const std::vector<uint32_t>& active) {
        idle(c, active, p / 10);
        idle(c, active, 2 * p);
    };
    head.safe_append_u("R", data);
    flip(head, data, false, 2 * p);
    slow_idle(head, data);
    head.safe_append_u("TICK", {});
    if (memory_x) h(head, data);

    stim::Circuit cycle;
    // Fixed measurement order: basis X/Z, color R/G/B, all a then all b.
    std::vector<uint32_t> check_order(2 * nf), flag_order(2 * nf);
    uint32_t order = 0;
    for (int basis = 0; basis < 2; ++basis) {
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
            cycle.safe_append_u("R", active); flip(cycle, active, false, 2 * p);
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
            flip(cycle, active, false, 5 * p);
            cycle.safe_append_u("M", aa);
            cycle.safe_append_u("M", bb);
            noise(cycle, "DEPOLARIZE1", active, p);
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
                if (!first || (basis == 0) == memory_x) {
                    std::vector<uint32_t> records{rec(back)};
                    if (!first) records.push_back(rec(back + m));
                    c.safe_append_u("DETECTOR", records, {x, y, 0, (double)faces[k].color, (double)basis, 0});
                }
                c.safe_append_u("DETECTOR", {rec(m - flag_order[basis * nf + k])},
                                {x, y, 0, (double)faces[k].color, (double)basis, 1});
            }
    };
    head += cycle;
    detectors(head, true);
    stim::Circuit body;
    body.safe_append_u("SHIFT_COORDS", {}, {0, 0, 1});
    body += cycle;
    detectors(body, false);
    head += body * (rounds - 1);
    if (memory_x) h(head, data);
    flip(head, data, false, 5 * p);
    head.safe_append_u("M", data);
    noise(head, "DEPOLARIZE1", data, p);
    slow_idle(head, data);
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

stim::Circuit cc666_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x,
                         const ComponentErrorDistribution& noise_distribution) {
    return apply_component_error_distribution(
        make_color_memory(distance, rounds, p, is_memory_x, false), noise_distribution);
}
stim::Circuit cc488_si1000(uint32_t distance, uint32_t rounds, double p, bool is_memory_x,
                         const ComponentErrorDistribution& noise_distribution) {
    return apply_component_error_distribution(
        make_color_memory(distance, rounds, p, is_memory_x, true), noise_distribution);
}
