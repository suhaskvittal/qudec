/* Contributor: OpenAI Codex. */

#include "circuit_generator.h"
#include "decoder/surface_code.h"
#include "decoder/uf.h"

#include <stim/simulators/dem_sampler.h>
#include <stim/util_top/circuit_to_dem.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{

void
require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void
compare_assignments(const dec::Result& ordinary, const dec::Result& clustered)
{
    const auto& a = ordinary.matching_data.assignments;
    const auto& b = clustered.matching_data.assignments;
    require(a.size() == b.size(), "Cluster output changed assignment count");
    dec::ObsType parity(ordinary.flipped_obs.num_bits_padded());
    for (size_t i = 0; i < a.size(); i++)
    {
        require(a[i].d1 == b[i].d1 && a[i].d2 == b[i].d2 &&
                a[i].frame_flips == b[i].frame_flips,
                "Cluster output changed detector pairing or assignment observable");
        parity ^= b[i].frame_flips;
    }
    require(parity == clustered.flipped_obs, "Assignment observable parity is inconsistent");
}

void
check_decoder_outputs(const stim::DetectorErrorModel& dem, size_t count, bool exhaustive)
{
    dec::PyMatching ordinary(dem), cluster(dem, {.enable_cluster_gap=true});
    dec::PyMatching exact(dem, {.enable_gap_estimation=true});
    dec::PyMatching endpoints(dem, {.report_matching_data=true});
    dec::PyMatching cluster_endpoints(dem, {.report_matching_data=true, .enable_cluster_gap=true});
    dec::UnionFind uf(dem), uf_cluster(dem, {.enable_cluster_gap=true});
    dec::UnionFind uf_endpoints(dem, {.report_matching_data=true});
    dec::UnionFind uf_cluster_endpoints(dem, {.report_matching_data=true, .enable_cluster_gap=true});

    stim::DemSampler<64> sampler(dem, std::mt19937_64(41), count);
    sampler.resample(false);
    auto syndromes = sampler.det_buffer.transposed();
    auto truths = sampler.obs_buffer.transposed();
    dec::SyndromeType synthetic(dem.count_detectors());
    dec::ObsType zero(dem.count_observables()), opposite_truth(dem.count_observables());
    for (size_t shot = 0; shot < count; shot++)
    {
        if (exhaustive)
            for (size_t v = 0; v < dem.count_detectors(); v++)
                synthetic[v] = (shot >> v) & 1;
        dec::SyndromeRef syndrome = exhaustive ? dec::SyndromeRef(synthetic) : syndromes[shot];
        dec::ObsRef truth = exhaustive ? dec::ObsRef(zero) : truths[shot];
        opposite_truth[0] = !bool(truth[0]);

        const auto a = ordinary.decode(syndrome, truth);
        const auto b = cluster.decode(syndrome, truth);
        const auto complement = exact.decode(syndrome, truth);
        const auto opposite = cluster.decode(syndrome, opposite_truth);
        require(!a.has_cluster_gap && b.has_cluster_gap, "MWPM gap availability is wrong");
        require(a.flipped_obs == b.flipped_obs, "Cluster output changed MWPM frame");
        require(b.cluster_gap_db == opposite.cluster_gap_db, "MWPM cluster output depends on truth");
        require(std::isfinite(b.cluster_gap_db) && b.cluster_gap_db >= 0,
                "Invalid MWPM cluster value");
        require(b.cluster_gap_db <= std::abs(complement.matching_data.gap) + 1e-8,
                "MWPM cluster value exceeds exact gap");
        require(b.decode_latency_ns > 0, "Cluster decoding latency is missing");

        const auto e = endpoints.decode(syndrome, truth);
        const auto f = cluster_endpoints.decode(syndrome, truth);
        require(e.flipped_obs == f.flipped_obs && f.cluster_gap_db == b.cluster_gap_db,
                "MWPM endpoint mode changed cluster output");
        compare_assignments(e, f);

        const auto u = uf.decode(syndrome, truth);
        const auto v = uf_cluster.decode(syndrome, truth);
        const auto opposite_uf = uf_cluster.decode(syndrome, opposite_truth);
        require(u.flipped_obs == v.flipped_obs && v.has_cluster_gap && !u.has_cluster_gap,
                "UF frame or gap availability changed");
        require(v.cluster_gap_db == opposite_uf.cluster_gap_db, "UF cluster output depends on truth");
        require(std::isfinite(v.cluster_gap_db) && v.cluster_gap_db >= 0,
                "Invalid UF cluster value");
        const auto ue = uf_endpoints.decode(syndrome, truth);
        const auto ve = uf_cluster_endpoints.decode(syndrome, truth);
        require(ue.flipped_obs == ve.flipped_obs && ve.cluster_gap_db == v.cluster_gap_db,
                "UF endpoint mode changed cluster output");
        compare_assignments(ue, ve);
    }

    synthetic.clear();
    const auto empty = cluster.decode(synthetic, zero);
    const auto empty_exact = exact.decode(synthetic, zero);
    require(std::abs(empty.cluster_gap_db - std::abs(empty_exact.matching_data.gap)) < 1e-8,
            "Empty syndrome cluster distance disagrees with logical-chain distance");
}

void
check_partial_edges()
{
    dec::ClusterGapGraph graph(3);
    graph.add_edge(graph.boundary_one, 0, 4);
    graph.add_edge(0, 1, 6);
    graph.add_edge(1, graph.boundary_zero, 3);
    graph.add_edge(2, graph.boundary_zero, 0.25);
    graph.finish();
    require(!graph.relevant[2], "Opposing basis leaked through the physical boundary");
    require(graph.gap(std::vector<double>(4)) == 13, "Empty cluster logical distance wrong");
    const auto& coverage = graph.ball_coverage({{0, 2}, {1, 1}, {2, 100}});
    require(coverage[0] == 2 && coverage[1] == 3 && coverage[2] == 1,
            "Partial weighted-ball coverage is wrong");
    require(graph.gap(coverage) == 7, "Weighted cluster gap is wrong");
    require(graph.gap({1.5, 2.5, 0.5, 0}) == 8.5, "UF partial growth rounded to whole edges");
    require(graph.gap(graph.ball_coverage({{0, 4}, {1, 3}})) == 0,
            "Fully covered logical path must have zero cluster gap");
}

} // namespace

int
main()
{
    check_partial_edges();
    for (const auto* text : {
            "error(0.1) D0 L0\nerror(0.2) D0 D1\nerror(0.3) D1\n",
            "error(0.5) D0 L0\nerror(0.2) D0 D1\nerror(0.3) D1\n",
            "error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\n",
            "error(0.1) D0 L0\nerror(0.2) D0 D1\nerror(0.3) D1 D2\n"
            "error(0.15) D2 D0\nerror(0.05) D2 D3\nerror(0.12) D3\n"})
    {
        stim::DetectorErrorModel dem(text);
        check_decoder_outputs(dem, size_t{1} << dem.count_detectors(), true);
    }
    for (double p : {0.001, 0.02})
        for (int d : {3, 5, 9})
            for (bool opposing : {false, true})
            {
                auto dem = stim::circuit_to_dem(sc_si1000(d, d, p, false, opposing),
                        {.decompose_errors=true});
                check_decoder_outputs(dem, 32, false);
            }

    for (const auto* text : {
            "error(0.1) D0 D1 L0\nerror(0.1) D1\n",
            "error(0.1) L0\nerror(0.1) D0\n",
            "error(0.7) D0 L0\nerror(0.1) D0\n",
            "error(0.1) D0 L0\nerror(0.1) D0 L1\n"})
    {
        stim::DetectorErrorModel dem(text);
        bool pm_rejected = false, uf_rejected = false;
        try { dec::PyMatching decoder(dem, {.enable_cluster_gap=true}); }
        catch (const std::invalid_argument&) { pm_rejected = true; }
        try { dec::UnionFind decoder(dem, {.enable_cluster_gap=true}); }
        catch (const std::invalid_argument&) { uf_rejected = true; }
        require(pm_rejected && uf_rejected, "Unsupported geometry was accepted");
    }
    std::cout << "Cluster gap tests passed\n";
}
