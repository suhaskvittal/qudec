// Author: OpenAI Codex. Adaptive ensemble decisions and matching reuse.
#include "circuit_generator.h"
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"
#include <stim/simulators/dem_sampler.h>
#include <stim/util_top/circuit_to_dem.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>

namespace
{
void require(bool value,const char* message)
{if (!value) throw std::runtime_error(message);}
void equal(const dec::concat::Result& a,const dec::concat::Result& b)
{
    require(a.observable==b.observable && a.color==b.color && a.member==b.member &&
        a.events==b.events && std::abs(a.weight-b.weight)<1e-9,
        "Gating changed the expected baseline/full-ensemble correction");
}
auto samples(const stim::DetectorErrorModel& dem,size_t shots)
{
    stim::DemSampler<64> sampler(dem,std::mt19937_64(8726),shots);
    sampler.resample(false);
    return sampler.det_buffer.transposed();
}
std::vector<uint8_t> bits(const stim::simd_bit_table<64>& table,size_t shot,size_t n)
{
    std::vector<uint8_t> result(n);
    for (size_t j=0;j<n;j++) result[j]=table[shot][j];
    return result;
}
void capacity()
{
    for (bool octagons:{false,true}) for (double p:{.01,.05})
    {
        const auto dem=stim::circuit_to_dem(octagons ? cc488_code_capacity(7,p,false):cc666_code_capacity(7,p,false),{});
        dec::concat::CircuitDecoder adapter(dem,false);
        std::vector<uint8_t> vacuum(dem.count_detectors());
        const auto model=adapter.inspect(vacuum).model;
        dec::concat::Decoder single(model,false,false,{1,20261024});
        dec::concat::Decoder full(model,false,false,{20,20261024});
        dec::concat::Decoder gated(model,false,false,{20,20261024});
        const auto table=samples(dem,256);
        size_t expanded=0,skipped=0;
        for (size_t shot=0;shot<256;shot++)
        {
            const auto syndrome=bits(table,shot,vacuum.size());
            const auto baseline=single.decode(syndrome),ensemble=full.decode(syndrome);
            const auto expected_calls=full.statistics.second_calls;
            const auto never=gated.decode_gated(syndrome,0);
            require(!never.used_ensemble && gated.statistics.first_calls==3,"Never-expand gate ran extra members");
            equal(never.result,baseline);
            const auto always=gated.decode_gated(syndrome,1000000);
            require(always.used_ensemble && gated.statistics.first_calls==60,"Always-expand gate did not reuse member zero");
            require(gated.statistics.second_calls==expected_calls,"Expansion repeated coalesced second stages");
            equal(always.result,ensemble);
            const auto selected=gated.decode_gated(syndrome,20);
            equal(selected.result,selected.used_ensemble ? ensemble:baseline);
            require(selected.baseline_color==baseline.color,"Gap belongs to a different baseline color");
            require(std::abs(selected.gap_db-never.gap_db)<1e-8,"Threshold changed the measured gap");
            selected.used_ensemble ? ++expanded:++skipped;
        }
        const auto empty=gated.decode_gated(vacuum,20);
        require(std::abs(empty.gap_db+10*std::log10(p))<1e-8,"Vacuum did not retain its physical corner gap");
        if (p==.01)
            require(!empty.used_ensemble && expanded && skipped,"20dB tie or selective expansion is incorrect");
        else require(empty.used_ensemble && !skipped,"20dB did not expand below the 5% corner cap");
    }
}
void flags()
{
    for (bool octagons:{false,true}) {
    const auto dem=stim::circuit_to_dem(octagons ? cc488_si1000(5,5,.001,false):cc666_si1000(5,5,.001,false),{});
    const auto table=samples(dem,64);
    std::set<uint64_t> ids;
    for (size_t d=0;d<dem.count_detectors();d++) ids.insert(d);
    const auto coordinates=dem.get_detector_coordinates(ids);
    bool saw_flag=false;
    for (bool renormalize:{false,true})
    {
        const dec::concat::FlagRenormalization policy{renormalize,.001};
        dec::concat::CircuitDecoder single(dem,false,true,16,false,{1,20261024},policy);
        dec::concat::CircuitDecoder full(dem,false,true,16,false,{5,20261024},policy);
        dec::concat::CircuitDecoder gated(dem,false,true,16,false,{5,20261024},policy,true);
        for (size_t shot=0;shot<64;shot++)
        {
            const auto syndrome=bits(table,shot,dem.count_detectors());
            for (const auto& [id,c]:coordinates) saw_flag|=c[5]==1 && c[4]==0 && syndrome[id];
            const auto baseline=single.inspect(syndrome).correction;
            const auto ensemble=full.inspect(syndrome).correction;
            const auto never=gated.decode_gated(syndrome,0);
            equal(never.result,baseline);
            const auto selected=gated.decode_gated(syndrome,1000000);
            equal(selected.result,selected.used_ensemble ? ensemble:baseline);
        }
    }
    require(saw_flag,"Flag-update regression bank contained no relevant flags");
    }
}
void wrapper()
{
    const auto dem=stim::circuit_to_dem(cc666_code_capacity(3,.01,false),{});
    dec::SyndromeType vacuum(dem.count_detectors());dec::ObsType truth(1);
    dec::ConcatMWPM gated(dem,{.ensemble_size=20,.debug=true,.gap_gating=true});
    auto result=gated.decode(vacuum,truth);
    require(result.has_cluster_gap && std::abs(result.cluster_gap_db-20)<1e-8,"Missing result gap metadata");
    require(gated.statistics().gap_evaluations==1 && gated.statistics().gap_expansions==0,"Wrong expansion counters");
    require(result.debug_info.find("Stage 1 matching")!=std::string::npos &&
        result.debug_info.find("Stage 2 matching")!=std::string::npos &&
        result.debug_info.find("Baseline winner:")!=std::string::npos &&
        result.debug_info.find("Winner:")!=std::string::npos,"Gating discarded decoder decisions from debug output");
    dec::ConcatMWPM k1(dem,{.gap_gating=true});
    result=k1.decode(vacuum,truth);
    require(!result.has_cluster_gap && k1.statistics().gap_evaluations==0,"K1 unnecessarily calculated a gating gap");
    for (double threshold:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
    {
        bool rejected=false;
        try {dec::ConcatMWPM invalid(dem,{.gap_gating=true,.gap_gating_threshold_db=threshold});}
        catch (const std::invalid_argument&) {rejected=true;}
        require(rejected,"Invalid gap threshold accepted");
    }
    for (bool parallel:{false,true})
    {
        bool rejected=false;
        try {dec::ConcatMWPM invalid(dem,{.boundary_gap=!parallel,.decoding_cores=parallel ? 2u:1u,.gap_gating=true});}
        catch (const std::invalid_argument&) {rejected=true;}
        require(rejected,"Unsupported gap-gating combination accepted");
    }
}
}
int main()
{
    capacity();flags();wrapper();
    std::cout<<"Gap gating: baseline/ensemble equivalence, matching reuse, thresholds, flag updates and debug passed\n";
}
