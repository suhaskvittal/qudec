#include "circuit_generator.h"
#include "decoder/libra.h"
#include "decoder/libra_variant.h"
#include <stim.h>
#include <stim/util_top/circuit_to_dem.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        uint64_t checked=0;
        for(int d:{5,7,9}) {
            const auto dem=stim::circuit_to_dem(sc_si1000(d,d,1e-3,false,true),{.decompose_errors=true});
            const size_t nd=dem.count_detectors();
            for(bool force_ensemble:{false,true}) {
                dec::Libra::Config oldcfg;oldcfg.seed=20261001;oldcfg.validate=true;
                oldcfg.gap_threshold_db=force_ensemble?1000:20;
                dec::Libra original(dem,oldcfg);
                dec::LibraVariant::Config cfg;cfg.seed=oldcfg.seed;cfg.validate=true;cfg.gap_threshold_db=oldcfg.gap_threshold_db;
                dec::LibraVariant preserved(dem,cfg);
                cfg.large_gap_fallback=false;dec::LibraVariant nofallback(dem,cfg);
                cfg.correlated_matching=false;dec::LibraVariant uncorrelated(dem,cfg);
                std::mt19937_64 rng(42+d+force_ensemble);
                stim::DemSampler<64> sampler(dem,std::move(rng),2048);sampler.resample(false);
                const auto dets=sampler.det_buffer.transposed();
                const size_t shots=force_ensemble?32:2048;
                for(size_t i=0;i<shots;++i){
                    const auto a=original.decode_detailed(dets[i]);
                    const auto b=preserved.decode_detailed(dets[i]);
                    const auto c=nofallback.decode_detailed(dets[i]),u=uncorrelated.decode_detailed(dets[i]);
                    if(a.prediction!=b.prediction||a.errors!=b.errors||a.gap_db!=b.gap_db||a.weights[0]!=b.weights[0]
                        ||a.weights[1]!=b.weights[1]||a.log_masses[0]!=b.log_masses[0]||a.log_masses[1]!=b.log_masses[1])
                        throw std::runtime_error("Default variant differs from preserved Libra");
                    if(c.fallback_used)throw std::runtime_error("Disabled fallback called");
                    if(b.gap_db<cfg.gap_threshold_db&&(b.prediction!=c.prediction||b.errors!=c.errors))
                        throw std::runtime_error("Small-gap branch changed");
                    if(b.gap_db>=cfg.gap_threshold_db&&c.prediction!=c.initial_prediction)
                        throw std::runtime_error("Cheaper representative not returned");
                    // Validation inside decode reconstructs both classes; also check returned assignments.
                    for(auto* decoder:{&preserved,&nofallback,&uncorrelated}) {
                        const auto& result=decoder==&preserved?b:decoder==&nofallback?c:u;
                        auto [syn,obs]=decoder->syndrome_of(result.errors);
                        if(obs!=result.prediction)throw std::runtime_error("Returned parity wrong");
                        for(size_t k=0;k<nd;++k)if(syn[k]!=dets[i][k])throw std::runtime_error("Returned syndrome wrong");
                    }
                    ++checked;
                }
            }
        }
        // Exact test: no correlations means both APIs must produce identical candidates and synthesis.
        const stim::DetectorErrorModel independent(R"(
            error(0.025) D0 L0
            error(0.04) D1
            error(0.015) D0 D1
            error(0.035) D2
            error(0.022) D3
            error(0.013) D2 D3
        )");
        dec::LibraVariant::Config cfg;cfg.validate=true;cfg.large_gap_fallback=false;cfg.gap_threshold_db=1000;
        dec::LibraVariant corr(independent,cfg);cfg.correlated_matching=false;dec::LibraVariant unc(independent,cfg);
        for(unsigned mask=1;mask<16;++mask){dec::SyndromeType syn(4);for(unsigned k=0;k<4;++k)syn[k]=(mask>>k)&1;
            auto a=corr.decode_detailed(syn),b=unc.decode_detailed(syn);
            if(a.prediction!=b.prediction||a.errors!=b.errors||std::abs(a.weights[0]-b.weights[0])>1e-9||std::abs(a.weights[1]-b.weights[1])>1e-9)
                throw std::runtime_error("Uncorrelated DEM disagrees between matching APIs");}
        std::cout<<"PASS: "<<checked<<" paired surface-code shots, forced ensembles, both-class syndrome checks, 15 independent-DEM syndromes\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
