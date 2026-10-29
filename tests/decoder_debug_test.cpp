/* Author: OpenAI Codex. Opt-in traces must preserve actual decoder decisions. */
#include "circuit_generator.h"
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"
#include "error_chain.h"
#include <stim/util_top/circuit_to_dem.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
void require(bool condition,const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
size_t occurrences(const std::string& text,const std::string& needle)
{
    size_t count = 0,offset = 0;
    while ((offset=text.find(needle,offset)) != std::string::npos) {count++; offset+=needle.size();}
    return count;
}
void flags()
{
    const stim::DetectorErrorModel dem(R"(
        error(0.0001) D0 D1 D2 L0
        error(0.0003) D0 D1 D2 D3 L0
        error(0.003) D3
        error(0.001) D3
        error(0.007) D4
        detector(0,0,1,0,1,0) D0
        detector(1,0,1,1,1,0) D1
        detector(2,0,1,2,1,0) D2
        detector(0,0,1,99,0,1) D3
        detector(0,0,1,-1,1,1) D4
    )");
    dec::ConcatMWPM baseline(dem), traced(dem,{.debug=true});
    dec::SyndromeType syndrome(5);
    dec::ObsType truth(1);
    syndrome[0] = syndrome[1] = syndrome[2] = syndrome[3] = syndrome[4] = 1;
    const auto original = baseline.decode(syndrome,truth);
    const auto result = traced.decode(syndrome,truth);
    require(original.debug_info.empty(),"Debug enabled by default");
    require(result.flipped_obs[0] == original.flipped_obs[0],"Trace changed flag correction");
    const auto& text = result.debug_info;
    require(occurrences(text,"Stage 1 matching:")==3 && occurrences(text,"Stage 2 matching:")==3,"Missing six matchings");
    for (const auto color : {"!R","!G","!B"}) require(text.find(color)!=std::string::npos,"Missing restricted lattice");
    require(text.find("match ")!=std::string::npos && text.find("path:")!=std::string::npos,"Missing matched paths");
    require(text.find("Virtual detector definitions")!=std::string::npos,"Missing virtual detector definitions");
    require(text.find("Winner:")!=std::string::npos,"Missing selected candidate");
    for (const auto* removed : {"Flags:","p_M=","outside_flag_penalty=","symmetric_difference=","selected_flags=","expected_flags=","[F,"})
        require(text.find(removed)==std::string::npos,"Flag metadata still present in debug_info");
    syndrome.clear();
    const auto empty = traced.decode(syndrome,truth);
    require(empty.debug_info.find("empty matching")!=std::string::npos,"Vacuum matching trace missing after cache change");
    require(!empty.flipped_obs[0],"Trace changed vacuum correction");
    require(traced.decode(syndrome,truth).debug_info == empty.debug_info,"Cached debug output is not deterministic");

    dec::ConcatMWPM ensemble(dem,{.ensemble_size=3,.debug=true});
    syndrome[0] = syndrome[1] = syndrome[2] = syndrome[3] = 1;
    const auto candidates = ensemble.decode(syndrome,truth);
    require(candidates.debug_info.find("member=2")!=std::string::npos,"Missing perturbed matching calls");
    require(candidates.debug_info.find("Coalesced:")!=std::string::npos,"Missing duplicate first-stage solution");
}
void native_models()
{
    for (bool octagons : {false,true}) for (bool pheno : {false,true})
    {
        const auto circuit = pheno
            ? (octagons ? cc488_phenomenological(3,3,.001,.001,false) : cc666_phenomenological(3,3,.001,.001,false))
            : (octagons ? cc488_code_capacity(3,.001,false) : cc666_code_capacity(3,.001,false));
        const auto dem = stim::circuit_to_dem(circuit,{});
        dec::ConcatMWPM plain(dem,{.ensemble_size=3});
        dec::ConcatMWPM traced(dem,{.ensemble_size=3,.debug=true});
        dec::concat::CircuitDecoder source(dem,false);
        dec::ObsType truth(1);
        for (size_t event = 0; event < source.faults.size(); event++)
        {
            const auto bits = source.fault_syndrome({event});
            dec::SyndromeType syndrome(bits.size());
            for (size_t v = 0; v < bits.size(); v++) syndrome[v] = bits[v];
            require(plain.decode(syndrome,truth).flipped_obs[0] == traced.decode(syndrome,truth).flipped_obs[0],"Trace changed native matching decision");
        }
    }
}
void chain_output()
{
    const stim::DetectorErrorModel dem("error(0.1) D0 L0");
    std::ostringstream out;
    evaluate_error_chains(dem,{.length=1,.samples=1,.print_chains=1},[](auto,auto)
    {
        dec::Result result;
        result.debug_info = "decoder-specific sentinel\n";
        return result;
    },out);
    require(out.str().find("Decoder decisions:\ndecoder-specific sentinel")!=std::string::npos,"Chain report drops debug_info");
}
}
int main()
{
    try {flags(); native_models(); chain_output(); std::cout << "Decoder debug tests passed\n";}
    catch (const std::exception& error) {std::cerr << error.what() << '\n'; return 1;}
}
