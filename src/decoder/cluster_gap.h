/* Contributor: OpenAI Codex. */
#ifndef DECODER_CLUSTER_GAP_h
#define DECODER_CLUSTER_GAP_h

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace stim { struct DetectorErrorModel; }

namespace dec
{

void validate_cluster_gap_dem(const stim::DetectorErrorModel&);

// Weighted metric graph with inequivalent physical boundaries kept separate.
// Supports one planar logical observable carried only by boundary edges.
class ClusterGapGraph
{
public:
    struct Edge
    {
        size_t u, v;
        double length;
    };
    std::vector<Edge> edges;
    std::vector<std::vector<std::pair<size_t, size_t>>> adjacency;
    std::vector<bool> relevant;
    size_t boundary_zero, boundary_one;

    explicit ClusterGapGraph(size_t num_detectors);

    void add_edge(size_t u, size_t v, double length);

    void finish();

    // Source radii are in the same units as edge lengths. The ball construction
    // allows growth through collapsed boundaries, conservatively enlarging the
    // covered region. Partial coverage from both ends of each edge is retained.
    const std::vector<double>& ball_coverage(const std::vector<std::pair<size_t, double>>& sources);

    // Scratch buffers belong to one decoder instance and retain their capacity.
    std::vector<double>& reset_coverage();

    double gap(const std::vector<double>& coverage);

private:
    using Item = std::pair<double, size_t>;
    std::vector<double> distance_, shortest_, coverage_;
    std::vector<Item> heap_;

    void push(double distance, size_t vertex);

    Item pop();
};

} // namespace dec
#endif
