/*
 *  author: OpenAI Codex
 *  date:   21 October 2026
 *
 *  Spatial restriction matching, connected components, and local lifting.
 * */
#ifndef DECODER_RESTRICTION_LATTICE_h
#define DECODER_RESTRICTION_LATTICE_h

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <vector>

namespace dec::restriction
{

using det_id_type = int32_t;
using edge_type = std::pair<det_id_type, det_id_type>;
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

struct Assignment
{
    // Index in the assignment vector is the unique assignment identity.
    int ignored_color;
    det_id_type x, y;
    std::vector<edge_type> path;
};

struct ConnectedComponent
{
    int color;
    std::vector<size_t> assignments;
};

std::vector<ConnectedComponent> build_connected_components(
    const std::vector<Vertex>& vertices, const std::array<det_id_type, 3>& boundaries,
    const std::vector<Assignment>& assignments, int reference_color);

struct Diagnostics
{
    size_t assignments{0};
    size_t components{0};
    size_t local_lifts{0};
    size_t max_multiplicity{0};
};

struct Result
{
    uint64_t observable_mask{0};
    std::vector<size_t> correction;
    Diagnostics diagnostics;
};

class RestrictionDecoder
{
public:
    explicit RestrictionDecoder(Lattice, size_t ensemble_size = 16, uint64_t seed = 0);
    Result decode(const syndrome_type&);

    // Deliberately public stages for independent validation and inspection.
    std::vector<Assignment> restricted_lattice_matching(const syndrome_type&);
    std::vector<ConnectedComponent> connected_components(const std::vector<Assignment>&,
                                                        int reference_color = 0) const;
    Result lifting(const syndrome_type&, const std::vector<Assignment>&,
                   const std::vector<ConnectedComponent>&, int reference_color = 0) const;

    const Lattice lattice;
private:
    std::vector<std::unique_ptr<RestrictionDecoder>> alternatives_;
    std::array<pm::Mwpm, 3> matchers_;
    std::array<std::map<edge_type, edge_type>, 3> boundary_edges_;
    // Enumerated local sections are small: 2^4, 2^6, or 2^8.
    struct LocalSection
    {
        std::vector<edge_type> boundary;
        std::vector<size_t> faces;
        double weight;
    };
    std::vector<std::vector<LocalSection>> sections_;
};

} // namespace dec::restriction
#endif
