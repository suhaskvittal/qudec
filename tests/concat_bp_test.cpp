// Author: OpenAI Codex. Exact flag marginals and posterior-model integration.
#include "circuit_generator.h"
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"
#include <stim/simulators/dem_sampler.h>
#include <stim/util_top/circuit_to_dem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
void require(bool value,const char* message)
{ if (!value) throw std::runtime_error(message); }
void exact_stars()
{
    dec::concat::FlagBP bp(2,{{.1,{0}},{.2,{0}},{.3,{1}},{.07,{}}});
    auto odd=bp.posterior({1,0}),even=bp.posterior({0,1});
    require(bp.independent_stars(),"Expected independent flag stars");
    require(std::abs(odd[0]-.08/.26)<1e-12,"Incorrect odd-flag posterior");
    require(std::abs(odd[1]-.18/.26)<1e-12,"Incorrect competing posterior");
    require(std::abs(even[0]-.02/.74)<1e-12,"Quiet flags did not condition faults");
    require(std::abs(odd[3]-.07)<1e-12,"Unflagged prior changed");
    require(even[2]>1-1e-10,"Certain flag event was not enabled");
    // A two-flag fault links two factors. Exact enumeration gives P(a|11).
    dec::concat::FlagBP correlated(2,{{.1,{0,1}},{.2,{0}},{.3,{1}}});
    const auto posterior=correlated.posterior({1,1});
    const double denominator=.1*.8*.7+.9*.2*.3;
    require(!correlated.independent_stars(),"Multi-flag fallback not used");
    require(std::abs(posterior[0]-(.1*.8*.7)/denominator)<1e-6,"Multi-flag tree posterior incorrect");
}
void auxiliary_check_conditioning()
{
    using namespace dec::concat;
    // The .2 pure-auxiliary event competes with the .1 correlated event.
    // Quiet and fired opposite-check outcomes have analytic star posteriors.
    const stim::DetectorErrorModel dem(
        "error(0.1) D0 D1\nerror(0.2) D1\nerror(0.07) D0\nerror(0.3) D2\n"
        "detector(0,0,0,0,1,0) D0\ndetector(0,0,0,1,0,0) D1\n"
        "detector(0,0,0,-1,0,1) D2\nlogical_observable L0\n");
    CircuitDecoder bp(dem,false,true,16,false,{}, {},false,20,true);
    std::vector<uint8_t> quiet(3),fired{0,1,0};
    const auto even=bp.inspect(quiet).model,odd=bp.inspect(fired).model;
    require(even.events.size()==1 && odd.events.size()==1,"Pure auxiliary event became a matching edge");
    require(std::abs(even.events[0].probability-(.02/.74+.07))<1e-12,
        "Quiet opposite check or pure auxiliary fault ignored");
    require(std::abs(odd.events[0].probability-(.08/.26+.07))<1e-12,
        "Fired opposite check or pure auxiliary fault ignored");
    require(std::abs(bp.inspect(quiet).model.events[0].probability-even.events[0].probability)<1e-12,
        "Opposite-check context was stale after changing only an auxiliary check");
}
void integration(bool memory_x,bool opposite=false)
{
    using namespace dec::concat;
    const auto dem=stim::circuit_to_dem(cc666_si1000(5,5,.001,memory_x,opposite),{});
    CircuitDecoder legacy(dem,memory_x),bp(dem,memory_x,true,16,false,{3,20261024},{},false,20,true);
    std::vector<uint8_t> zero(dem.count_detectors());
    require(bp.decode(zero)==0,"Noiseless BP decoding failed");
    size_t x_flags=0,z_flags=0;
    for (auto flag:bp.flag_globals)
    {
        x_flags+=bp.coordinates.at(flag)[4]==0;
        z_flags+=bp.coordinates.at(flag)[4]==1;
    }
    require(x_flags && z_flags,"BP must retain both extraction flag bases");
    for (auto flag:legacy.flag_globals)
        require(legacy.coordinates.at(flag)[4]==(memory_x ? 1:0),"Legacy flag filtering changed");
    auto extra=zero;
    for (const auto& [id,c]:bp.coordinates) if (c[5]==1 && c[4]==(memory_x ? 0:1)) extra[id]=1;
    const auto vacuum=bp.inspect(zero).model,conditioned=bp.inspect(extra).model;
    require(vacuum.events.size()==conditioned.events.size(),"Flag conditioning changed model topology");
    bool changed=false;
    for (size_t i=0;i<vacuum.events.size();i++)
        changed |= std::abs(vacuum.events[i].probability-conditioned.events[i].probability)>1e-12;
    require(changed,"Same-basis flags did not condition BP probabilities");
    if (opposite)
    {
        auto auxiliary_only=zero;
        for (const auto& [id,c]:bp.coordinates)
            if (c[5]==0 && c[4]==(memory_x ? 1:0)) {auxiliary_only[id]=1;break;}
        const auto with_check=bp.inspect(auxiliary_only).model;
        bool different=false;
        for (size_t i=0;i<vacuum.events.size();i++)
            different |= std::abs(vacuum.events[i].probability-with_check.events[i].probability)>1e-12;
        require(different,"Opposite check-only change did not invalidate posterior context");
        const auto restored=bp.inspect(zero).model;
        for (size_t i=0;i<vacuum.events.size();i++)
            require(std::abs(restored.events[i].probability-vacuum.events[i].probability)<1e-12,
                "Returning to quiet auxiliary outcomes did not restore posterior context");
    }
    stim::DemSampler<64> sampler(dem,std::mt19937_64(20261100),64);
    sampler.resample(false);const auto table=sampler.det_buffer.transposed();
    for (size_t shot=0;shot<16;shot++)
    {
        auto syndrome=zero;for (size_t d=0;d<zero.size();d++) syndrome[d]=table[shot][d];
        const auto result=bp.inspect(syndrome);
        // Independently reconstruct original auxiliary-signature variables:
        // marginalizing opposite-basis checks before inference would be wrong.
        std::vector<FlagBPVariable> variables;
        const auto& classes=bp.equivalence_classes();
        std::map<std::vector<size_t>,size_t> ids;
        for (size_t i=0;i<classes.size();i++) ids[classes[i].checks]=i;
        std::vector<size_t> local(dem.count_detectors(),SIZE_MAX);
        for (size_t i=0;i<bp.check_globals.size();i++) local[bp.check_globals[i]]=i;
        using Key=std::tuple<size_t,uint64_t,std::vector<size_t>>;
        std::map<Key,double> priors;
        for (const auto& fault:bp.faults)
        {
            std::vector<size_t> checks,auxiliary;
            for (auto d:fault.detectors)
                if (local[d]!=SIZE_MAX) checks.push_back(local[d]);
                else auxiliary.push_back(d);
            if (checks.empty() && auxiliary.empty()) continue;
            std::sort(checks.begin(),checks.end());
            std::sort(auxiliary.begin(),auxiliary.end());
            auto& prior=priors[{checks.empty() ? SIZE_MAX:ids.at(checks),fault.observable,auxiliary}];
            prior+=fault.probability-2*prior*fault.probability;
        }
        for (const auto& [key,p]:priors)
        {
            const auto& [id,observable,auxiliary]=key;
            variables.push_back({p,auxiliary,id,observable});
        }
        FlagBP reference_bp(dem.count_detectors(),variables);
        const auto posterior=reference_bp.posterior(syndrome);
        std::vector<std::map<uint64_t,double>> frames(classes.size());
        for (size_t i=0;i<variables.size();i++)
        {
            const auto& v=variables[i];
            if (v.class_id==SIZE_MAX) continue;
            double p=posterior[i];
            // Independent parity formula remains a cross-check for star graphs.
            if (reference_bp.independent_stars() && !v.flags.empty())
            {
                double product=1;
                for (size_t j=0;j<variables.size();j++) if (i!=j && variables[j].flags==v.flags)
                    product*=1-2*variables[j].probability;
                const double prior=v.probability,all=(1-2*prior)*product;
                const double exact=prior*(1+(syndrome[v.flags[0]] ? 1:-1)*product)/(1+(syndrome[v.flags[0]] ? -1:1)*all);
                require(std::abs(p-exact)<1e-11,"Star inference differs from independent parity formula");
            }
            frames[v.class_id][v.observable]+=p;
        }
        Model reference;reference.colors=result.model.colors;
        for (size_t i=0;i<classes.size();i++)
        {
            const auto best=std::max_element(frames[i].begin(),frames[i].end(),[](const auto& a,const auto& b){return a.second<b.second;});
            reference.events.push_back({std::clamp(best->second,1e-12,1-1e-9),classes[i].checks,best->first,classes[i].boundary_hint});
            require(std::abs(reference.events.back().probability-result.model.events[i].probability)<1e-11,"Integrated posterior differs from exact reference");
        }
        Decoder direct(reference,false,false,{3,20261024});
        std::vector<uint8_t> checks;for (auto d:bp.check_globals) checks.push_back(syndrome[d]);
        const auto expected=direct.decode(checks);
        require(expected.observable==result.correction.observable && std::abs(expected.weight-result.correction.weight)<1e-8,
            "Posterior update did not reach both matching stages");
    }
    bool rejected=false;
    try {CircuitDecoder invalid(dem,memory_x,true,16,false,{}, {true,.001},false,20,true);}
    catch (const std::invalid_argument&) {rejected=true;}
    require(rejected,"BP and renormalization conflict accepted");
}
}
int main() try
{
    exact_stars();auxiliary_check_conditioning();integration(false);integration(true);integration(false,true);integration(true,true);
    std::cout<<"BP exact marginals and integrated matching passed\n";
}
catch (const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
