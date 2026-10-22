/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Regressions for flag inference, full-F penalties, colorless flags,
 *  capacity matching, and the saved d7 Bell-flagged two-hook syndrome.
 * */
#include "circuit_generator.h"
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"
#include "decoder/restriction/lattice.h"
#include <stim/util_top/circuit_to_dem.h>

#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace
{
void
require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void
flag_weights(bool memory_x)
{
    const int basis = memory_x ? 0 : 1;
    std::ostringstream text;
    text << "error(0.0001) D0 D1 D2 L0\n"
         << "error(0.0003) D0 D1 D2 D3 L0\n"
         << "error(0.003) D3\nerror(0.001) D3\nerror(0.007) D4\n"
         << "error(0.02) D3 D5\n"; // Opposite check + flag is NOT isolated.
    for (int c = 0; c < 3; c++)
        text << "detector(0,0,1," << c << ',' << basis << ",0) D" << c << '\n';
    // Flag color metadata deliberately falls outside the check-color range.
    text << "detector(0,0,1,99," << 1-basis << ",1) D3\n"
         << "detector(0,0,1,-1," << basis << ",1) D4\n"
         << "detector(0,0,1,0," << 1-basis << ",0) D5\n";
    stim::DetectorErrorModel dem(text.str());
    dec::concat::CircuitDecoder core(dem,memory_x);
    const double pm = .003+.001-2*.003*.001;
    require(std::abs(core.flag_measurement_probabilities()[3]-pm) < 1e-15,
            "Flag prior includes non-isolated events or loses XOR merging");
    require(std::abs(core.flag_measurement_probabilities()[4]-.007) < 1e-15,
            "Per-flag prior was replaced by a common rate");
    const auto absent = core.inspect({0,0,0,0,0,0});
    const auto associated = core.inspect({1,1,1,1,0,0});
    const auto unrelated = core.inspect({0,0,0,0,1,0});
    const auto both = core.inspect({1,1,1,1,1,0});
    const double mixture = (1-pm)*.0003+pm*.0001;
    require(std::abs(absent.weights.at(0).effective_probability-.0001) < 1e-15,
            "Unflagged prior changed with empty flag syndrome");
    require(std::abs(associated.weights.at(0).effective_probability-mixture) < 1e-15,
            "Single-flag mixture is incorrect");
    require(std::abs(unrelated.weights.at(0).effective_probability-.0001*.007) < 1e-15,
            "Unrelated observed flag failed to penalize unflagged error");
    require(std::abs(both.weights.at(0).effective_probability-mixture*.007) < 1e-15,
            "Full observed flag set was not used");
    require(associated.audited_graph_edges == 6 && both.audited_graph_edges == 6,
            "Six matching graphs failed the actual weight audit");
    dec::ConcatMWPM decoder(dem,{.memory_x=memory_x,.flag_cache_capacity=1});
    dec::SyndromeType syndrome(6);
    syndrome[0] = syndrome[1] = syndrome[2] = syndrome[3] = syndrome[4] = 1;
    dec::ObsType truth(1);
    require(decoder.decode(syndrome,truth).flipped_obs[0],"Integrated flag decode failed");
    truth[0] = 1;
    require(decoder.decode(syndrome,truth).flipped_obs[0],"Decoder depends on supplied logical truth");
    dec::SyndromeType empty(6);
    require(!decoder.decode(empty,truth).flipped_obs[0],"Vacuum decode failed");
    require(decoder.decode(syndrome,truth).flipped_obs[0],"Cache eviction changed the result");
    require(decoder.statistics().shots == 4,"Statistics missed decoded shots");
}

void
capacity()
{
    const auto lattice = dec::restriction::make_code_capacity_lattice(false,5,.01);
    dec::concat::Model model;
    for (size_t v = 0; v < lattice.num_checks; v++) model.colors.push_back(lattice.vertices[v].color);
    for (const auto& face : lattice.faces)
    {
        dec::concat::Event event{face.probability,{},face.observable_mask};
        for (auto v : face.vertices) if (!lattice.vertices[v].is_boundary) event.detectors.push_back(v);
        model.events.push_back(event);
    }
    dec::concat::Decoder decoder(model);
    std::ostringstream text;
    for (const auto& event : model.events)
    {
        text << "error(" << event.probability << ")";
        for (auto v : event.detectors) text << " D" << v;
        if (event.observable) text << " L0";
        text << '\n';
    }
    for (size_t v = 0; v < model.colors.size(); v++)
        text << "detector(0,0,0," << model.colors[v] << ",1,0) D" << v << '\n';
    dec::ConcatMWPM integrated(stim::DetectorErrorModel(text.str()));
    auto check = [&](const std::vector<size_t>& errors)
    {
        const auto syndrome = lattice.syndrome(errors);
        const auto expected = lattice.observable(errors);
        require(decoder.decode(syndrome).observable == expected,"Capacity matching failed");
        dec::SyndromeType bits(syndrome.size());
        for (size_t v = 0; v < syndrome.size(); v++) bits[v] = syndrome[v];
        dec::ObsType unused_truth(1);
        require(bool(integrated.decode(bits,unused_truth).flipped_obs[0]) == bool(expected),
                "Integrated code-capacity decode failed");
    };
    for (size_t a = 0; a < lattice.faces.size(); a++)
    {
        check({a});
        for (size_t b = a+1; b < lattice.faces.size(); b++) check({a,b});
    }
    require(decoder.audit_weights() > 0,"Capacity graphs were not audited");
}

void
circuit_regression()
{
    const auto dem = stim::circuit_to_dem(cc666_si1000(7,7,.0005,false),{});
    std::set<uint64_t> ids;
    for (size_t v = 0; v < dem.count_detectors(); v++) ids.insert(v);
    const auto coordinates = dem.get_detector_coordinates(ids);
    dec::SyndromeType syndrome(dem.count_detectors());
    auto flip = [&](double x, double y, int time, int color, int basis, int flag)
    {
        for (const auto& [v,c] : coordinates)
            if (std::abs(c[0]-x) < 1e-9 && std::abs(c[1]-y) < 1e-9 && c[2] == time &&
                c[3] == color && c[4] == basis && c[5] == flag)
            {
                syndrome[v] ^= 1;
                return;
            }
        throw std::runtime_error("Missing saved circuit-regression detector");
    };
    // Original DEM events 3150 + 3214: surviving blue Z checks and two flags.
    flip(.875,.649519052838329,2,2,1,0);
    flip(5,3.4641016151377544,2,2,1,0);
    flip(2,1.7320508075688776,2,0,0,1);
    flip(3.5,2.598076211353316,2,1,0,1);
    flip(.875,.649519052838329,2,2,0,0);
    flip(3.5,.8660254037844386,2,2,0,0);
    flip(3.5,2.598076211353316,3,1,0,0);
    dec::concat::CircuitDecoder core(dem,false);
    std::vector<uint8_t> full(dem.count_detectors());
    for (size_t v = 0; v < full.size(); v++) full[v] = syndrome[v];
    const auto detail = core.inspect(full);
    for (const auto& candidate : detail.colors)
        require(candidate.observable == 0 && candidate.events.size() == 2,
                "Two-hook regression failed in a color subdecoder");
    dec::ConcatMWPM decoder(dem);
    dec::ObsType truth(1);
    auto result = decoder.decode(syndrome,truth);
    require(!result.flipped_obs[0] && !result.matching_data.frame_flips[0],
            "Integrated two-hook flag regression failed");
    bool rejected = false;
    try
    {
        dec::SyndromeType short_syndrome(64);
        decoder.decode(short_syndrome,truth);
    }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected,"Short syndrome was accepted");
}
} // namespace

int
main()
{
    flag_weights(false);
    flag_weights(true);
    capacity();
    circuit_regression();
    std::cout << "Concatenated MWPM regressions passed\n";
}
