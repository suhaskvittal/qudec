/* Author: OpenAI Codex. Chamberland et al., NJP 22, 023019, Sec. III D. */
#include "decoder/concat_mwpm/circuit.h"
#include "decoder/concat_mwpm.h"
#include <cmath>
#include <iomanip>
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
    const stim::DetectorErrorModel dem(R"(
        error(0.01) D0 D1 D2 L0
        error(0.002) D0 D1 D2 D3 D4 L0
        error(0.03) D0 D1
        error(0.004) D1 D2 D5
        error(0.005) D0 D2 D3 D4
        error(0.001) D3
        error(0.001) D4
        error(0.001) D5
        detector(0,0,0,0,1,0) D0
        detector(1,0,0,1,1,0) D1
        detector(2,0,0,2,1,0) D2
        detector(0,0,0,-1,0,1) D3
        detector(0,0,0,-1,0,1) D4
        detector(0,0,1,-1,0,1) D5
    )");
    dec::concat::CircuitDecoder core(dem,false,true,16,{}, {true,.001});
    std::vector<uint8_t> syndrome(core.num_detectors);
    auto weight = [](const auto& detail,std::vector<size_t> support) -> const dec::concat::FlagWeight&
    {
        for (const auto& w : detail.weights) if (w.checks == support) return w;
        throw std::runtime_error("Missing class");
    };
    auto detail = core.inspect(syndrome);
    require(weight(detail,{0,1,2}).effective_probability == .01,"No-flag ordinary graph changed");
    require(weight(detail,{1,2}).selected_variant == SIZE_MAX,"Quiet flag-only class enabled");
    syndrome[3] = 1;
    detail = core.inspect(syndrome);
    require(std::abs(weight(detail,{0,1,2}).effective_probability-.00001)<1e-16,
            "Unflagged representative did not receive p^m suppression");
    require(weight(detail,{0,2}).selected_variant == SIZE_MAX,"Partially observed flag pattern enabled");
    syndrome[4] = 1;
    detail = core.inspect(syndrome);
    require(weight(detail,{0,1,2}).effective_probability == .002 &&
            weight(detail,{0,1,2}).conditional_denominator == 1,
            "Supported event was divided by a flag marginal");
    require(std::abs(weight(detail,{0,1}).effective_probability-.00003)<1e-16,
            "Two flags from one measurement circuit counted twice");
    syndrome[5] = 1;
    detail = core.inspect(syndrome);
    require(std::abs(weight(detail,{0,1}).effective_probability-.00000003)<1e-18,
            "Distinct flagged rounds did not produce p squared");
    require(weight(detail,{0,1,2}).effective_probability == .002,
            "Unrelated observed flag suppressed an edge in C_f");
    syndrome[0] = syndrome[1] = syndrome[2] = 1;
    std::ostringstream out;
    out << std::setprecision(12);
    dec::concat::DebugTrace trace{out,{}};
    detail = core.inspect(syndrome,&trace);
    require(detail.audited_graph_edges > 0,"No actual graph weights audited");
    require(out.str().find("p=0.002, W=") != std::string::npos,
            "Projected supported edge included unsupported contributions");
    dec::ConcatMWPM wrapper(dem,{.flag_renormalization=true,.physical_error_rate=.001});
    dec::SyndromeType input(core.num_detectors);
    for (size_t v=0;v<syndrome.size();v++) input[v]=syndrome[v];
    dec::ObsType truth(1);
    require(wrapper.decode(input,truth).flipped_obs[0] == static_cast<uint64_t>(detail.correction.observable),
            "Wrapper lost renormalization configuration");
    bool rejected = false;
    try { dec::concat::CircuitDecoder invalid(dem,false,true,16,{}, {true,-.1}); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected,"Invalid physical probability accepted");
}
}
int main()
{
    try { tests(); std::cout << "Paper flag renormalization passed\n"; }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
