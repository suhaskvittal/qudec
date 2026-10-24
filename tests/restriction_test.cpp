/*
 *  author: OpenAI Codex
 *  date:   21 October 2026
 *
 *  Regression: a blue detector must lift toward the BLUE boundary when its
 *  circuit event is the pair of faces B-C-D and BlueBoundary-C-D.
 * */
#include "circuit_generator.h"
#include "decoder/restriction.h"

#include <stim/util_top/circuit_to_dem.h>

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>

int
main()
{
    // Contributor: OpenAI Codex. The saved fixture explicitly includes X checks
    // in a Z-memory experiment, so request the optional opposite-basis detectors.
    const auto dem = stim::circuit_to_dem(cc666_si1000(7,7,0.0005,false,true), {});
    std::set<uint64_t> ids;
    for (size_t v = 0; v < dem.count_detectors(); v++) ids.insert(v);
    const auto coordinates = dem.get_detector_coordinates(ids);
    dec::SyndromeType syndrome(dem.count_detectors());
    const auto flip = [&](double x, double y, int time, int color, int basis, int flag)
    {
        for (const auto& [v,c] : coordinates)
            if (std::abs(c[0]-x) < 1e-9 && std::abs(c[1]-y) < 1e-9 && c[2] == time
                && c[3] == color && c[4] == basis && c[5] == flag)
            {
                syndrome[v] ^= 1;
                return;
            }
        throw std::runtime_error("Missing regression detector");
    };
    // Raw event 6148: A_Z(3), B_X(4), A_X(4), C_X(4), A_Z(4).
    flip(3.5,4.330127018922193,3,0,1,0);
    flip(5,3.4641016151377544,4,2,0,0);
    flip(3.5,4.330127018922193,4,0,0,0);
    flip(5,5.196152422706632,4,1,0,0);
    flip(3.5,4.330127018922193,4,0,1,0);
    // Raw event 7005: flag_D,X(4), C_X(4), B_Z(4), D_X(5).
    flip(6.125,4.113620667976083,4,0,0,1);
    flip(5,5.196152422706632,4,1,0,0);
    flip(5,3.4641016151377544,4,2,1,0);
    flip(6.125,4.113620667976083,5,0,0,0);

    dec::restriction::CircuitRestrictionDecoder core(dem,
        dec::restriction::make_code_capacity_lattice(false,7,0.0005),false,0.0005,true);
    dec::restriction::syndrome_type bits(dem.count_detectors());
    for (size_t v = 0; v < bits.size(); v++) bits[v] = syndrome[v];
    const auto result = core.decode(bits);
    if (result.observable_mask != 0 || result.lifted.correction != std::vector<size_t>{24,29})
        throw std::runtime_error("Incorrect blue-boundary hook lift");
    dec::RestrictionDecoder decoder(dem,{.distance=7,.physical_p=0.0005,.use_flags=true});
    dec::ObsType truth(1);
    if (decoder.decode(syndrome,truth).flipped_obs[0])
        throw std::runtime_error("Incorrect integrated decoder frame");
    std::cout << "Blue-boundary hook regression passed: correction={24,29}, frame=0\n";
}
