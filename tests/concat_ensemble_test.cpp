/*
 *  author: OpenAI Codex
 *  date:   23 October 2026
 *
 *  First-stage ensembling: candidate coalescing, unperturbed scoring,
 *  fixed-seed prefixes, full spacetime syndromes and flag-cache rebuilding.
 * */
#include "circuit_generator.h"
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"
#include <stim/util_top/circuit_to_dem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>

namespace
{
void require(bool value,const char* message)
{
    if (!value) throw std::runtime_error(message);
}
const dec::concat::Result& prefix(const std::vector<dec::concat::Result>& candidates,size_t k)
{
    const dec::concat::Result* out = nullptr;
    for (const auto& candidate : candidates)
        if (candidate.member < k && (!out || candidate.weight < out->weight-1e-9)) out = &candidate;
    require(out,"Empty ensemble prefix");
    return *out;
}
void model_tests(bool pheno,bool boundary)
{
    const auto circuit = pheno ? cc666_phenomenological(7,7,.01,.01,false)
                               : cc666_code_capacity(7,.03,false);
    const auto dem = stim::circuit_to_dem(circuit,{});
    dec::concat::CircuitDecoder source(dem,false,true,16,boundary);
    const auto model = source.inspect(std::vector<uint8_t>(source.num_detectors)).model;
    dec::concat::Decoder full(model,false,boundary,{60,20261024,1});
    dec::concat::Decoder small(model,false,boundary,{20,20261024,1});
    dec::concat::Decoder repeated(model,false,boundary,{60,20261024,1});
    dec::concat::Decoder identical(model,false,boundary,{20,20261024,0});
    require(full.audit_weights()>0,"No ensemble graph weights audited");
    std::mt19937_64 rng(1729);
    for (size_t shot = 0; shot < 40; shot++)
    {
        std::vector<uint8_t> measured(source.num_detectors);
        for (const auto& fault : source.faults)
            if (std::uniform_real_distribution<double>()(rng) < fault.probability)
                for (auto v : fault.detectors) measured[v] ^= 1;
        std::vector<uint8_t> local;
        for (auto v : source.check_globals) local.push_back(measured[v]);
        const auto candidates = full.decode_candidates(local);
        const auto smaller = small.decode_candidates(local);
        const auto rerun = repeated.decode_candidates(local);
        const auto unperturbed = source.inspect(measured).correction;
        const auto& baseline = prefix(candidates,1);
        require(baseline.observable==unperturbed.observable && baseline.events==unperturbed.events,
                "Unperturbed ensemble member changed baseline correction");
        require(prefix(candidates,20).events==prefix(smaller,20).events,
                "Seeded ensemble prefix changed correction");
        require(prefix(candidates,60).events==prefix(rerun,60).events,
                "Repeated seeded ensemble changed correction");
        require(prefix(candidates,60).weight <= baseline.weight+1e-8,
                "Ensemble increased best original correction weight");
        const auto unchanged = identical.decode_candidates(local);
        require(identical.decode(local).observable==unperturbed.observable,
                "Zero-spread ensemble changed logical result");
        require(std::all_of(unchanged.begin(),unchanged.end(),[](const auto& c){return c.member==0;}),
                "Zero-spread members were not coalesced");
        std::array<std::set<std::vector<size_t>>,3> keys;
        for (const auto& candidate : candidates)
        {
            require(keys[candidate.color].insert(candidate.virtual_detectors).second,
                    "Duplicate virtual syndrome was decoded twice");
            auto residual = local;
            double score = 0;
            uint64_t logical = 0;
            for (auto id : candidate.events)
            {
                const auto& event = full.model.events[id];
                for (auto v : event.detectors) residual[v] ^= 1;
                logical ^= event.observable;
                score -= std::log(event.probability);
            }
            require(std::none_of(residual.begin(),residual.end(),[](auto x){return x!=0;}),
                    "Ensemble correction fails full syndrome");
            require(std::abs(score-candidate.weight)<1e-8 && logical==candidate.observable,
                    "Ensemble used perturbed final score");
        }
    }
    // Individual measurement faults retain both time endpoints; no flattening.
    for (size_t id = 0; id < source.faults.size(); id++)
    {
        const auto full_syndrome = source.fault_syndrome({id});
        std::vector<uint8_t> local;
        for (auto v : source.check_globals) local.push_back(full_syndrome[v]);
        require(full.decode(local).observable==source.fault_observable({id}),"Single DEM fault failed");
    }
}
void additional_models()
{
    // Exercise the second tiling and real Bell-flag contexts, including eviction.
    const auto square_dem = stim::circuit_to_dem(cc488_code_capacity(5,.01,false),{});
    dec::concat::CircuitDecoder squares(square_dem,false,true,1,false,{20,20261024,1});
    for (size_t id = 0; id < squares.faults.size(); id++)
        require(squares.decode(squares.fault_syndrome({id}))==squares.fault_observable({id}),
                "cc488 ensemble single fault failed");
    const auto bell_dem = stim::circuit_to_dem(cc666_si1000(5,5,.0005,false,true),{});
    dec::concat::CircuitDecoder bell(bell_dem,false,true,1,false,{20,20261024,1});
    size_t checked = 0;
    for (size_t id = 0; id < bell.faults.size() && checked < 6; id++)
    {
        const auto& event = bell.faults[id];
        bool flag = false, check = false;
        for (auto v : event.detectors)
        {
            flag |= bell.coordinates.at(v)[5]==1;
            check |= bell.coordinates.at(v)[5]==0 && bell.coordinates.at(v)[4]==1;
        }
        if (!flag || !check) continue;
        const auto syndrome = bell.fault_syndrome({id});
        const auto expected = bell.fault_observable({id});
        require(bell.decode(syndrome)==expected,"Bell-flag ensemble single fault failed");
        require(bell.decode(std::vector<uint8_t>(bell.num_detectors))==0,"Bell vacuum decode failed");
        require(bell.decode(syndrome)==expected,"Bell flag-cache rebuild changed correction");
        checked++;
    }
    require(checked==6,"No Bell flag-conditioned faults tested");
}
void flag_cache()
{
    const stim::DetectorErrorModel dem(
        "error(0.0001) D0 D1 D2 L0\nerror(0.0003) D0 D1 D2 D3 L0\nerror(0.003) D3\n"
        "detector(0,0,1,0,1,0) D0\ndetector(1,0,1,1,1,0) D1\n"
        "detector(0,1,1,2,1,0) D2\ndetector(0,0,1,99,0,1) D3\n");
    dec::concat::CircuitDecoder core(dem,false,true,1,false,{20,20261024,1});
    const std::vector<uint8_t> fired{1,1,1,1}, empty{0,0,0,0};
    const auto first = core.inspect(fired).correction;
    require(core.decode(empty)==0 && core.inspect(fired).correction.events==first.events,
            "Flag-cache rebuild changed ensemble correction");
    require(first.observable==1,"Colorless flag conditioned ensemble incorrectly");
    dec::ConcatMWPM wrapped(dem,{.flag_cache_capacity=1,.ensemble_size=20,.ensemble_seed=20261024});
    dec::SyndromeType syndrome(4);
    syndrome[0]=syndrome[1]=syndrome[2]=syndrome[3]=1;
    dec::ObsType truth(1);
    require(wrapped.decode(syndrome,truth).flipped_obs[0],"Wrapper lost ensemble options");
    bool rejected = false;
    try { dec::concat::CircuitDecoder invalid(dem,false,true,1,false,{0,0,1}); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected,"Zero-sized ensemble accepted");
}
}
int main()
{
    try
    {
        for (bool pheno : {false,true}) for (bool boundary : {false,true}) model_tests(pheno,boundary);
        flag_cache();
        additional_models();
        std::cout << "PASS first-stage ensemble, spacetime, prefix, coalescing and flag-cache regressions\n";
    }
    catch (const std::exception& error) {std::cerr << error.what() << '\n'; return 1;}
}
