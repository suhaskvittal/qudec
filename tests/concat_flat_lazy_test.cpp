// Author: OpenAI Codex. Flat BP state and demand-driven ensemble preparation.
#include "circuit_generator.h"
#include "decoder/concat_mwpm/circuit.h"
#include <stim/util_top/circuit_to_dem.h>
#include <stim/simulators/dem_sampler.h>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

void require(bool ok,const char* message)
{ if (!ok) throw std::runtime_error(message); }
void same(const dec::concat::Result& a,const dec::concat::Result& b)
{
    require(a.observable==b.observable && a.events==b.events && a.color==b.color &&
        a.member==b.member && std::abs(a.weight-b.weight)<1e-9,"Lazy preparation changed correction");
    require(a.candidate_gap.has_value()==b.candidate_gap.has_value(),"Lazy preparation changed gap availability");
    if (a.candidate_gap) require(std::abs(*a.candidate_gap-*b.candidate_gap)<1e-9,"Lazy preparation changed gap");
}
int main() try
{
    const auto dem=stim::circuit_to_dem(cc666_si1000(5,5,.001,false),{});
    dec::concat::CircuitDecoder bp(dem,false,true,16,{5,20261024},{},true);
    stim::DemSampler<64> sampler(dem,std::mt19937_64(918783),64);
    sampler.resample(false);const auto shots=sampler.det_buffer.transposed();
    std::vector<uint8_t> zero(dem.count_detectors());
    const auto topology=bp.inspect(zero).model;
    dec::concat::Decoder lazy(topology,false,{5,20261024},true);
    dec::concat::Decoder full(topology,false,{5,20261024},true);
    dec::concat::Decoder single(topology,false,{},true);
    for (size_t i=0;i<64;i++)
    {
        std::vector<uint8_t> bits(zero.size()),memory;
        for (size_t d=0;d<bits.size();d++) bits[d]=shots[i][d];
        for (auto d:bp.check_globals) memory.push_back(bits[d]);
        const auto inspected=bp.inspect(bits);
        std::vector<double> probabilities;std::vector<dec::concat::Frame> frames;
        for (const auto& event:inspected.model.events)
        { probabilities.push_back(event.probability);frames.push_back(event.observable); }
        // Dense update must preserve full paths/corrections versus the general model path.
        lazy.update_dense(probabilities,frames);
        single.update(inspected.model);
        const auto first=lazy.decode(memory);
        require(lazy.statistics.ensemble_updates==12,"Used ensemble graphs were not prepared lazily");
        same(first,single.decode(memory));
        lazy.update_dense(probabilities,frames);full.update(inspected.model);
        const auto expanded=lazy.decode(memory);
        require(lazy.statistics.ensemble_updates==12,"Ensemble graphs were not refreshed after dense weight updates");
        same(expanded,full.decode(memory));
        // Repeat with unchanged context: no graph should need updating twice.
        same(lazy.decode(memory),full.decode(memory));
        require(lazy.statistics.ensemble_updates==0,"Repeated context rebuilt prepared ensemble graphs");
        require(inspected.weights.size()==inspected.model.events.size(),"On-demand BP diagnostics lost classes");
        lazy.audit_weights();
    }
    require(bp.inspect(zero).correction.observable==0,"Returning to quiet flags left stale BP state");
    std::cout<<"Flat BP and lazy ensemble graph updates: 64 changing-flag syndromes passed\n";
}
catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
