/* Author: OpenAI Codex. Independent connected-chain and parity regressions. */
#include "error_chain.h"
#include "error_chain_display.h"
#include "circuit_generator.h"
#include <stim/util_top/circuit_to_dem.h>
#include <set>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
void require(bool value,const char* message)
{
    if (!value) throw std::runtime_error(message);
}
void tests()
{
    require(error_chain_combination_count(100,50) == "100891344545564193334812497256","Large binomial count");
    require(error_chain_combination_count(5,0) == "1","Empty binomial count");
    require(error_chain_combination_count(5,6) == "0","Impossible binomial count");
    const stim::DetectorErrorModel dem(R"(
        error(0.1) D0 L0
        error(0.1) D0 D1
        error(0.1) D1
        error(0.1) D5 L1
        error(0) D1
    )");
    size_t calls = 0;
    auto decode = [&](auto syndrome,auto truth)
    {
        calls++;
        require(!syndrome[5],"Disconnected event entered a pair");
        require(syndrome[0] != syndrome[1],"Wrong pair syndrome");
        require(truth[0] == syndrome[1] && !truth[1],"Wrong pair observable parity");
        return dec::Result{.flipped_obs=dec::ObsType(2)};
    };
    std::ostringstream out;
    const auto exhaustive = evaluate_error_chains(dem,{.length=2,.exhaustive=true},decode,out);
    require(exhaustive.combinations_visited == 6 && exhaustive.tested == 2 && exhaustive.failures == 1,"Wrong exhaustive counts");
    require(calls == 2,"Disconnected combinations decoded");
    require(out.str().find("Event 0 (p=0.1) flips detectors {D0}, observables {L0}") != std::string::npos,"Missing individual fault support");
    require(out.str().find("Event 1 (p=0.1) flips detectors {D0, D1}, observables {}") != std::string::npos,"Canceled detector missing from individual event");
    require(out.str().find("Total combinations before connectivity filtering: 6") < out.str().find("Failure 1"),"Count must precede decoding");
    std::ostringstream first,second;
    const auto random = evaluate_error_chains(dem,{.length=2,.samples=100,.seed=42},decode,first);
    evaluate_error_chains(dem,{.length=2,.samples=100,.seed=42},decode,second);
    require(random.tested == 100 && random.failures > 0 && random.failures < 100,"Random sampling failed");
    require(first.str() == second.str(),"Seed does not reproduce chains");
    bool rejected = false;
    try {evaluate_error_chains(dem,{.length=4},decode,out);}
    catch (const std::invalid_argument&) {rejected = true;}
    require(rejected,"Impossible connected chain accepted");

    // A separator is a decomposition hint, not another independently sampled fault.
    const stim::DetectorErrorModel separated("error(0.1) D0 L0 ^ D1 L0\nerror(0.1) D1 L1\n");
    const auto correlated = evaluate_error_chains(separated,{.length=2,.exhaustive=true},
        [&](auto syndrome,auto truth)
        {
            require(syndrome[0] && !syndrome[1],"Separator split a DEM event");
            require(!truth[0] && truth[1],"Observable XOR not respected");
            dec::Result result{.flipped_obs=dec::ObsType(2)};
            result.flipped_obs[1] = 1;
            return result;
        },out);
    require(correlated.tested == 1 && !correlated.failures,"Correlated event failed");
    const stim::DetectorErrorModel repeat("repeat 2 {\nerror(0.1) D0 L0\nshift_detectors 1\n}\n");
    const auto singles = evaluate_error_chains(repeat,{.length=1,.exhaustive=true},
        [&](auto syndrome,auto truth)
        {
            require(syndrome[0] != syndrome[1] && truth[0],"Flattened repeat/shift incorrect");
            dec::Result result{.flipped_obs=dec::ObsType(1)};
            result.flipped_obs[0] = 1;
            return result;
        },out);
    require(singles.tested == 2 && !singles.failures,"Flattened events lost");
    // Author: OpenAI Codex. Shared flags persist while check/observable parity cancels.
    const stim::DetectorErrorModel flagged(R"(
        error(0.1) D0 D1 L0
        error(0.1) D0 D1 D2 L0
        detector(0,0,0,0,1,0) D0
        detector(0,0,0,-1,0,1) D1
        detector(1,0,0,1,1,0) D2
    )");
    auto verify_flags = [](auto syndrome,auto truth)
    {
        require(!syndrome[0] && syndrome[1] && syndrome[2] && !truth[0],
                "Chain flags must use OR, with XOR for checks and observables");
        return dec::Result{.flipped_obs=dec::ObsType(1)};
    };
    evaluate_error_chains(flagged,{.length=2,.exhaustive=true},verify_flags,out);
    evaluate_error_chains(flagged,{.length=2,.samples=10},verify_flags,out);
}
void display_tests()
{
    for (const std::string code : {"cc666","cc488"}) for (uint32_t distance : {3,7,11,15})
    {
        const CodeCapacityNoisePolicy noise(.001);
        const auto circuit = code == "cc666" ? cc666_circuit(distance,1,noise,false) : cc488_circuit(distance,1,noise,false);
        const auto dem = stim::circuit_to_dem(circuit,{});
        ErrorChainDisplay display(dem,code,distance);
        dec::SyndromeType syndrome(dem.count_detectors());
        syndrome[0] = 1;
        std::ostringstream drawing;
        display.draw(syndrome,drawing);
        require(drawing.str().find("Color-code lattice ("+code+")") != std::string::npos,"Missing lattice drawing");
        require(display.detector_label(0).find("(0, 0):Z") != std::string::npos,"Wrong check identity");
        require(drawing.str().find(display.detector_label(0)) != std::string::npos,"Missing detector legend");
    }
    const auto dem = stim::circuit_to_dem(cc666_phenomenological(3,3,.001,.001,false,true),{});
    std::set<uint64_t> ids;
    for (size_t id = 0; id < dem.count_detectors(); id++) ids.insert(id);
    const auto coordinates = dem.get_detector_coordinates(ids);
    ErrorChainDisplay display(dem,"cc666",3);
    dec::SyndromeType syndrome(dem.count_detectors());
    size_t flipped = 0;
    for (const auto& [id,c] : coordinates)
        if (c[0] == coordinates.at(0)[0] && c[1] == coordinates.at(0)[1] && c[4] == 1 && (c[2] == 0 || c[2] == 2))
        {
            syndrome[id] = 1;
            require(display.detector_label(id).find("(0, "+std::to_string(static_cast<int>(c[2]))+"):Z") != std::string::npos,"Round changes spatial check ID");
            flipped++;
        }
    require(flipped == 2,"Did not exercise two distinct rounds");
    std::ostringstream drawing;
    display.draw(syndrome,drawing);
    require(drawing.str().find("(0, 0):Z") != std::string::npos && drawing.str().find("(0, 2):Z") != std::string::npos,"Multi-round legend incomplete");

    const auto bell_dem = stim::circuit_to_dem(cc666_si1000(3,2,.001,false),{});
    ids.clear();
    for (size_t id = 0; id < bell_dem.count_detectors(); id++) ids.insert(id);
    ErrorChainDisplay bell_display(bell_dem,"cc666",3);
    dec::SyndromeType bell_syndrome(bell_dem.count_detectors());
    bool found_flag = false;
    for (const auto& [id,c] : bell_dem.get_detector_coordinates(ids)) if (c[5] == 1)
    {
        bell_syndrome[id] = 1; found_flag = true;
        require(bell_display.detector_label(id).find("=F(") != std::string::npos,"Flag assigned a color");
        break;
    }
    require(found_flag,"No Bell flag tested");
    std::ostringstream flags;
    bell_display.draw(bell_syndrome,flags);
    require(flags.str().find("(colorless flag)") != std::string::npos,"Missing flag legend");
}

}
int main()
{
    try {tests(); display_tests(); std::cout << "Connected DEM chain tests passed\n";}
    catch (const std::exception& error) {std::cerr << error.what() << '\n'; return 1;}
}
