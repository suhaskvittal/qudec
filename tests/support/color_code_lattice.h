/*
 *  author: OpenAI Codex
 *  date:   21 October 2026
 *
 *  Color-code lattice fixture for matching tests.
 * */
#ifndef TESTS_COLOR_CODE_LATTICE_H
#define TESTS_COLOR_CODE_LATTICE_H

#include <stim.h>
#include <array>
#include <cstdint>
#include <vector>

namespace test_support
{

using det_id_type = int32_t;
using syndrome_type = std::vector<uint8_t>;
struct Vertex
{
    int color;
    det_id_type canonical_detector;
    bool is_boundary{false};
};

// A physical data qubit is a triangle of the dual lattice. Missing check
// colors are completed with the corresponding colored boundary vertex.
struct Face
{
    std::array<det_id_type, 3> vertices;
    uint64_t observable_mask{1};
    double probability;
};

struct Lattice
{
    size_t num_checks;
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<std::vector<size_t>> incident_faces;

    det_id_type boundary_id(int color) const;
    syndrome_type syndrome(const std::vector<size_t>& errors) const;
    uint64_t observable(const std::vector<size_t>& errors) const;
};

Lattice make_code_capacity_lattice(bool octagons, uint32_t distance, double p);

} // namespace test_support
#endif
