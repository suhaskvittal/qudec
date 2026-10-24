/* Author: OpenAI Codex. Detector-only classes and nearest-flag representative selection. */
#include "decoder/concat_mwpm/circuit.h"
#include "circuit_generator.h"
#include <stim/util_top/circuit_to_dem.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
void require(bool condition,const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
std::string metadata()
{
    return R"(
        detector(0,0,0,0,1,0) D0
        detector(1,0,0,1,1,0) D1
        detector(2,0,0,2,1,0) D2
        detector(0,0,0,99,0,1) D3
        detector(0,0,0,-1,0,1) D4
        detector(2,0,1,2,1,0) D5
        detector(0,0,0,99,0,1) D6
        detector(0,0,1,-1,0,1) D7
    )";
}
void mixed_class()
{
    dec::concat::CircuitDecoder decoder(stim::DetectorErrorModel(R"(
        error(0.01) D0 D1 D2 L0
        error(0.02) D0 D1 D2 D3
        error(0.03) D0 D1 D2 D4 L0
        error(0.001) D3
        error(0.002) D4
    )"+metadata()),false);
    const auto& classes = decoder.equivalence_classes();
    require(classes.size()==1 && decoder.class_count()==1,"Logical effects split a detector equivalence class");
    const auto& group = classes[0];
    require(group.checks==std::vector<size_t>({0,1,2}),"Wrong class support");
    require(group.flags==std::vector<size_t>({3,4}),"Class flag union incorrect");
    require(group.variants.size()==3,"Variants were discarded or combined across flag signatures");
    bool unflagged=false,flag3=false,flag4=false;
    for (const auto& variant : group.variants)
    {
        if (variant.flags.empty()) unflagged=variant.observable==1 && variant.probability==.01;
        if (variant.flags==std::vector<size_t>{3}) flag3=variant.observable==0 && variant.probability==.02;
        if (variant.flags==std::vector<size_t>{4}) flag4=variant.observable==1 && variant.probability==.03;
    }
    require(unflagged && flag3 && flag4,"Per-variant frame/flag/probability metadata not retained");
    std::vector<uint8_t> syndrome(decoder.num_detectors);
    syndrome[0]=syndrome[1]=syndrome[2]=1;
    auto verify = [&](double probability,uint64_t observable,std::vector<size_t> flags,size_t distance)
    {
        const auto detail=decoder.inspect(syndrome);
        require(detail.weights.size()==1,"Class split during representative selection");
        const auto& selected=detail.weights[0];
        require(selected.selected_flags==flags && selected.flag_distance==distance,"Incorrect flag representative/distance");
        double denominator=1;
        for(auto flag:flags) if(syndrome[flag]) denominator*=decoder.flag_flip_probabilities()[flag];
        const double conditional=std::min(probability/denominator,1.0);
        require(std::abs(selected.effective_probability-conditional)<1e-14 && selected.observable==observable,"Conditional probability/frame incorrect");
        require(selected.representative_probability==probability && selected.conditional_denominator==denominator,"Normalization diagnostics incorrect");
        require(detail.correction.observable==observable,"Matching did not use the representative's frame");
    };
    verify(.01,1,{},0);
    syndrome[3]=1;
    verify(.02,0,{3},0);
    syndrome[3]=0; syndrome[4]=1;
    verify(.03,1,{4},0);
    syndrome[3]=1;
    verify(.03,1,{4},1); // Equal flag distance: higher probability wins.
    syndrome[3]=syndrome[4]=0; syndrome[6]=1;
    verify(.01,1,{},1); // Unrelated observed flag favors the unflagged member.
    syndrome[3]=1;
    verify(.02,0,{3},1); // F includes a flag outside the class union.

}
void flagged_only_and_ties()
{
    dec::concat::CircuitDecoder decoder(stim::DetectorErrorModel(R"(
        error(0.01) D0 D1 D2 D3
        error(0.02) D0 D1 D2 D4 D6 L0
    )"+metadata()),false);
    std::vector<uint8_t> syndrome(decoder.num_detectors);
    const auto absent=decoder.inspect(syndrome);
    require(absent.weights[0].selected_variant==SIZE_MAX && absent.model.events.empty(),"Flagged-only class enabled with empty F");
    syndrome[3]=1;
    auto detail=decoder.inspect(syndrome);
    require(detail.weights[0].selected_flags==std::vector<size_t>{3},"Probability overrode a closer flag match");
    syndrome[3]=0;syndrome[6]=1;
    detail=decoder.inspect(syndrome);
    require(detail.weights[0].selected_flags==std::vector<size_t>({4,6}) && detail.weights[0].flag_distance==1,"Partial multi-flag match was gated off");
    syndrome[6]=0;syndrome[7]=1;
    detail=decoder.inspect(syndrome);
    require(detail.weights[0].selected_variant==SIZE_MAX && detail.model.events.empty(),"Unrelated observed flags enabled a flagged-only class");

    dec::concat::CircuitDecoder tied(stim::DetectorErrorModel(R"(
        error(0.02) D0 D1 D2 D3 L0
        error(0.02) D0 D1 D2 D4
    )"+metadata()),false);
    syndrome.assign(tied.num_detectors,0);syndrome[3]=syndrome[4]=1;
    const auto tie=tied.inspect(syndrome);
    require(tie.weights[0].selected_variant==0 && tie.weights[0].observable==tied.equivalence_classes()[0].variants[0].observable,"Equal-probability tie not stable");
    dec::concat::CircuitDecoder ignored(stim::DetectorErrorModel(R"(
        error(0.01) D0 D1 D2 L0
        error(0.4) D0 D1 D2 D3
    )"+metadata()),false,false);
    syndrome.assign(ignored.num_detectors,0);syndrome[3]=1;
    const auto disabled=ignored.inspect(syndrome);
    require(disabled.weights[0].selected_flags.empty() && disabled.weights[0].effective_probability==.01,"Disabled flags did not use the unflagged graph");
}

void extraction_basis_flags()
{
    // Author: OpenAI Codex. Same-extraction-basis flags remain irrelevant
    // under both conditional weights and the paper's p^m renormalization.
    for (bool paper : {false,true}) for (bool memory_x : {false,true})
    {
        const int basis=memory_x ? 0 : 1;
        std::ostringstream text;
        text << "error(0.01) D0 L0\nerror(0.02) D0 D1 L0\n"
             << "error(0.03) D0 D2 L0\nerror(0.04) D0 D1 D2 L0\n"
             << "error(0.004) D1\nerror(0.005) D2\n"
             << "detector(0,0,0,0," << basis << ",0) D0\n"
             << "detector(0,0,0,99," << 1-basis << ",1) D1\n"
             << "detector(0,0,0,-1," << basis << ",1) D2\n";
        dec::concat::CircuitDecoder core(stim::DetectorErrorModel(text.str()),memory_x,true,16,false,{}, {paper,.001});
        require(core.flag_globals==std::vector<size_t>{1},"Wrong check-extraction flag basis kept");
        require(core.equivalence_classes()[0].flags==std::vector<size_t>{1},"Ignored flag remained on a class member");
        require(core.flag_measurement_probabilities()[2]==0,"Ignored flag prior was retained");
        const auto absent=core.inspect({1,0,0}), ignored=core.inspect({1,0,1});
        const auto relevant=core.inspect({1,1,0}), both=core.inspect({1,1,1});
        const double unflagged=.01+.03-2*.01*.03;
        const double flagged=.02+.04-2*.02*.04;
        require(std::abs(absent.weights[0].effective_probability-unflagged)<1e-15,"Ignored flag channels not marginalized");
        require(ignored.weights[0].effective_probability==absent.weights[0].effective_probability &&
                ignored.weights[0].selected_flags.empty(),"Same-basis flag altered the empty-F model");
        const double expected = paper ? flagged : flagged/core.flag_flip_probabilities()[1];
        require(std::abs(relevant.weights[0].effective_probability-expected)<1e-14 &&
                relevant.weights[0].selected_flags==std::vector<size_t>{1},"Opposite extraction-block flag not used");
        require(both.weights[0].effective_probability==relevant.weights[0].effective_probability &&
                both.weights[0].flag_distance==relevant.weights[0].flag_distance,"Ignored flag altered representative distance");
        require(absent.correction.observable==ignored.correction.observable && relevant.correction.observable==both.correction.observable,
                "Ignored flag altered a decoding result");
    }
    for (bool octagons : {false,true}) for (bool memory_x : {false,true})
    {
        const auto circuit=octagons ? cc488_si1000(3,3,.001,memory_x) : cc666_si1000(3,3,.001,memory_x);
        dec::concat::CircuitDecoder core(stim::circuit_to_dem(circuit,{}),memory_x);
        size_t total_flags=0;
        for (const auto& [id,c] : core.coordinates) total_flags+=c[5]==1;
        require(!core.flag_globals.empty() && core.flag_globals.size()*2==total_flags,"Wrong count of retained native extraction-block flags");
        for (auto id : core.flag_globals)
            require(core.coordinates.at(id)[4]==(memory_x ? 1 : 0),"Native flag metadata interpreted as the wrong extraction block");
    }
}

void conditional_probabilities()
{
    // The multi-flag product approximation is explicitly checked, not silently
    // clipped into the probability range when strong correlations invalidate it.
    dec::concat::CircuitDecoder core(stim::DetectorErrorModel(R"(
        error(0.01) D0 D1 D2 D3 D4 L0
        error(0.2) D3
        error(0.3) D4
    )"+metadata()),false);
    std::vector<uint8_t> syndrome(core.num_detectors);syndrome[3]=syndrome[4]=1;
    const auto detail=core.inspect(syndrome);
    const double q3=.01+.2-2*.01*.2, q4=.01+.3-2*.01*.3;
    require(std::abs(core.flag_flip_probabilities()[3]-q3)<1e-15 && std::abs(core.flag_flip_probabilities()[4]-q4)<1e-15,"Incorrect flag marginals");
    require(std::abs(detail.weights[0].effective_probability-.01/(q3*q4))<1e-14,"Incorrect multiple-flag conditional approximation");
    dec::concat::CircuitDecoder invalid(stim::DetectorErrorModel("error(0.01) D0 D1 D2 D3 D4 L0\n"+metadata()),false);
    bool rejected=false;
    try {invalid.inspect(syndrome);}
    catch(const std::invalid_argument& e) {rejected=std::string(e.what()).find("exceeds one")!=std::string::npos;}
    require(rejected,"Invalid correlated-flag independence approximation was silently accepted");

    // p>1/2 remains a nonnegative cost and must be present in both matching
    // and search graphs without clipping it to one half.
    dec::concat::Model input{{0,1,2},{{.8,{0},1},{.1,{0,1,2},0},{.15,{1},0},{.12,{2},0}}};
    dec::concat::Decoder decoder(input);
    require(decoder.audit_weights()>0,"Conditional graph weights failed audit");
    for(size_t mask=0;mask<8;mask++)
    {
        std::vector<uint8_t> bits{uint8_t(mask&1),uint8_t((mask>>1)&1),uint8_t((mask>>2)&1)};
        const auto result=decoder.decode(bits);
        std::vector<uint8_t> actual(3);double weight=0;
        for(auto id:result.events)
        {
            const auto& e=decoder.model.events[id];
            for(auto v:e.detectors) actual[v]^=1;
            weight-=std::log(e.probability);
        }
        require(actual==bits && std::abs(weight-result.weight)<1e-12,"Conditional correction has incorrect syndrome/score");
    }

}

void conflicts()
{
    for (const std::string flags : {""," D3"," D3 D4"})
    {
        bool rejected=false;
        try
        {
            std::ostringstream text;
            text << "error(0.01) D0 D1 D2" << flags << " L0\nerror(0.02) D2 D0 D1"
                 << flags << '\n' << metadata();
            dec::concat::CircuitDecoder decoder(stim::DetectorErrorModel(text.str()),false);
        }
        catch (const std::invalid_argument& error)
        {
            const std::string text=error.what();
            rejected=text.find("different frame flips")!=std::string::npos
                && text.find("events 0 and 1")!=std::string::npos && text.find("distance <= 2")!=std::string::npos;
        }
        require(rejected,"Identical full detector signature with conflicting frames accepted");
    }
}
void merging_and_rounds()
{
    dec::concat::CircuitDecoder decoder(stim::DetectorErrorModel(R"(
        error(0.01) D0 D1 D2 D3 L0
        error(0.02) D2 D0 D1 D3 L0
        error(0.04) D0 D1 D5 D3 L0
        error(0) D0 D1 D2 D3
    )"+metadata()),false);
    require(decoder.equivalence_classes().size()==2,"Detector rounds were flattened in class key");
    const auto& first=decoder.equivalence_classes()[0];
    require(first.variants.size()==1,"Identical frame and flag signatures not merged");
    require(std::abs(first.variants[0].probability-(.01+.02-2*.01*.02))<1e-15,"Wrong parity probability for merged channels");
    require(first.variants[0].observable==1,"Merging changed the logical frame");
}
}
int main()
{
    try {mixed_class(); flagged_only_and_ties(); extraction_basis_flags(); conditional_probabilities(); conflicts(); merging_and_rounds(); std::cout << "Detector-only equivalence class construction passed\n";}
    catch (const std::exception& error) {std::cerr << error.what() << '\n'; return 1;}
}
