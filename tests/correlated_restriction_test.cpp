/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Fractional-weight, low-arity, flag, and public decoder integration regressions.
 * */
#include "circuit_generator.h"
#include "decoder/correlated_restriction.h"
#include "decoder/correlated_restriction/circuit.h"
#include <stim/util_top/circuit_to_dem.h>
#include <stim/simulators/frame_simulator_util.h>

#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

namespace
{
void
require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void
model(uint32_t rounds)
{
    using namespace dec::correlated_restriction;
    Decoder core(5,rounds,.01,true);
    dec::CorrelatedRestrictionDecoder decoder(5,rounds,.01);
    State all;
    for (const auto& event : core.events)
        for (int c = 0; c < 3; c++) all[c].insert(event.projection[c]);
    const auto weights = core.conditional_probabilities(all);
    bool triple = false;
    for (size_t id = 0; id < core.events.size(); id++)
    {
        const auto& event = core.events[id];
        if (event.checks.size() < 3)
        {
            for (int c = 0; c < 3; c++)
                require(weights[id][c] == event.probability,"A low-arity event changed");
        }
        else if (!triple)
        {
            State evidence;
            evidence[2].insert(event.projection[2]);
            require(core.conditional_probabilities(evidence)[id][2] == event.probability,
                    "Target graph used its own matching as evidence");
            const double w = std::log((1-event.probability)/event.probability);
            evidence[0].insert(event.projection[0]);
            double q = core.conditional_probabilities(evidence)[id][2];
            require(std::abs(std::log((1-q)/q)-2*w/3) < 1e-12,"One support must give 2W/3");
            evidence[1].insert(event.projection[1]);
            q = core.conditional_probabilities(evidence)[id][2];
            require(std::abs(std::log((1-q)/q)-w/3) < 1e-12,"Two supports must give W/3");
            triple = true;
        }
        const auto syndrome = core.syndrome({id});
        dec::SyndromeType bits(core.num_detectors);
        for (size_t v = 0; v < syndrome.size(); v++) bits[v] = syndrome[v];
        dec::ObsType truth(1);
        truth[0] = !core.observable({id}); // Deliberately wrong: must not affect decoding.
        const auto out = decoder.decode(bits,truth);
        require(out.flipped_obs[0] == bool(core.observable({id})),"Integrated single-event decode failed");
        require(out.flipped_obs[0] == out.matching_data.frame_flips[0],"Frame reporting disagrees");
    }
    require(triple,"Missing triple test");
    dec::SyndromeType zero(core.num_detectors);
    dec::ObsType truth(1);
    require(!decoder.decode(zero,truth).flipped_obs[0],"Vacuum decode failed");
    require(decoder.statistics().shots == core.events.size()+1,"Statistics lost shots");
}

void
circuit(bool memory_x)
{
    const double p = .0005;
    auto c = cc666_si1000(5,5,p,memory_x);
    auto dem = stim::circuit_to_dem(c,{});
    dec::CorrelatedRestrictionDecoder decoder(dem,{.distance=5,.physical_p=p,.memory_x=memory_x});
    dec::correlated_restriction::CircuitRestrictionDecoder core(dem,
        dec::restriction::make_code_capacity_lattice(false,5,p),memory_x,p,true,1,true);
    std::mt19937_64 rng(20261022);
    auto [dets,obs] = stim::sample_batch_detection_events<64>(c,128,rng);
    for (size_t shot = 0; shot < 128; shot++)
    {
        dec::SyndromeType bits(dem.count_detectors());
        dec::restriction::syndrome_type s(dem.count_detectors());
        for (size_t v = 0; v < s.size(); v++) bits[v] = s[v] = dets[v][shot];
        const auto paired = core.decode_single_pass(s);
        require(std::isfinite(paired[1].score),"Paired lift failed");
        dec::ObsType supplied_truth(1);
        supplied_truth[0] = !obs[0][shot];
        require(decoder.decode(bits,supplied_truth).flipped_obs[0] == bool(paired[1].frame),
                "Public wrapper does not use exactly one correlated pass");
    }
    auto clean = cc666_si1000(5,5,0,memory_x);
    auto clean_dem = stim::circuit_to_dem(clean,{});
    dec::CorrelatedRestrictionDecoder clean_decoder(clean_dem,{.distance=5,.physical_p=0,.memory_x=memory_x});
    dec::SyndromeType zero(clean_dem.count_detectors());
    dec::ObsType truth(1);
    truth[0] = 1;
    require(!clean_decoder.decode(zero,truth).flipped_obs[0],"Zero-noise circuit decode failed");
}
}

int
main()
{
    try
    {
        model(1);
        model(5);
        circuit(false);
        circuit(true);
        std::cout << "Correlated restriction regressions passed\n";
    }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
