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

#include <cmath>
#include <set>
#include <stdexcept>

namespace dec::restriction
{

namespace
{
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

} // namespace dec::restriction
