/*
 *  author: Claude Opus 5
 *  date:   21 September 2026
 * */

#include "decoder/bunchaluts/hypergraph.h"
#include "globals.h"

#include <algorithm>
#include <iostream>

namespace decoder
{
namespace bal
{
namespace hg
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

/*
 * Combines the independent probabilities of two error mechanisms that
 * produce the same detector pair into the probability that an odd number of
 * them fire, mirroring the merge used for the `BlossomV` decoding graph
 * (`decoder::BlossomV`'s `_update_adjacency_list` in `surface_code.cpp`).
 * */
double
_xor_combine(double pa, double pb)
{
    return pa*(1-pb) + pb*(1-pa);
}

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Hypergraph<basic_vdata, basic_edata, 2>
from_toric_code_dem(const stim::DetectorErrorModel& dem)
{
    using graph_type = Hypergraph<basic_vdata, basic_edata, 2>;
    using id_type = graph_type::id_type;

    graph_type g;

    const size_t num_detectors = dem.count_detectors();
    for (size_t i = 0; i < num_detectors; i++)
        g.add_vertex(basic_vdata{});

    dem.iter_flatten_error_instructions(
            [&g] (const auto& inst)
            {
                const double pr = inst.arg_data[0];
                inst.for_separated_targets(
                        [&g, pr] (const auto& grp)
                        {
                            std::vector<id_type> dets;
                            std::vector<size_t> frame_flips;
                            for (const auto& t : grp)
                            {
                                if (t.is_relative_detector_id())
                                    dets.push_back(static_cast<id_type>(t.val()));
                                else if (t.is_observable_id())
                                    frame_flips.push_back(static_cast<size_t>(t.val()));
                            }

                            if (dets.empty())
                                return;

                            if (dets.size() == 1)
                            {
                                std::cerr << "Hypergraph::from_toric_code_dem(): encountered a "
                                              "single-detector error term (detector " << dets[0] << "), "
                                              "but a toric-code DEM should have no boundary -- if this "
                                              "code family ever needs one, add a vertex with "
                                              "`is_boundary = true` and route boundary terms to it"
                                          << _die{};
                            }
                            if (dets.size() > 2)
                            {
                                std::cerr << "Hypergraph::from_toric_code_dem(): encountered an "
                                              "error term spanning " << dets.size() << " detectors, "
                                              "but this hypergraph has max order 2 -- pass a fully "
                                              "decomposed DEM (circuit_to_dem with decompose_errors=true)"
                                          << _die{};
                            }

                            std::sort(frame_flips.begin(), frame_flips.end());

                            auto existing = g.edges_incident_to({dets[0], dets[1]});
                            if (existing.empty())
                            {
                                g.add_edge({dets[0], dets[1]}, basic_edata{.error_probability=pr, .frame_flips=frame_flips});
                            }
                            else
                            {
                                id_type e = existing[0];
                                if (g.e(e).frame_flips != frame_flips)
                                {
                                    std::cerr << "Hypergraph::from_toric_code_dem(): warning: merging "
                                                  "two error terms between detectors " << dets[0] << " and "
                                              << dets[1] << " with differing frame flips\n";
                                }
                                g.e(e).error_probability = _xor_combine(g.e(e).error_probability, pr);
                            }
                        });
            });

    return g;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace hg
} // namespace bal
} // namespace decoder
