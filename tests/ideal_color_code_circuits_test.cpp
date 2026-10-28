/* Author: OpenAI Codex. Native code-capacity/phenomenological model regressions. */
#include "circuit_generator.h"
#include "decoder/concat_mwpm/circuit.h"
#include <stim/simulators/frame_simulator_util.h>
#include <stim/util_top/circuit_to_dem.h>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace
{
void require(bool value,const char* message) {if (!value) throw std::runtime_error(message);}
stim::Circuit generate(bool octagons,int d,int rounds,double p,double q,bool memory_x,bool opposite,
    const ComponentErrorDistribution& noise={})
{
    return octagons ? cc488_phenomenological(d,rounds,p,q,memory_x,opposite,noise)
                    : cc666_phenomenological(d,rounds,p,q,memory_x,opposite,noise);
}
void clean_models()
{
    for (bool octagons:{false,true}) for (int d:{3,5,7}) for (bool memory_x:{false,true}) for (bool opposite:{false,true})
    {
        const size_t n=octagons ? (d*d+2*d-1)/2:(3*d*d+1)/4,m=(n-1)/2;
        for (int rounds:{1,d})
        {
            const auto circuit=generate(octagons,d,rounds,0,0,memory_x,opposite);
            require(circuit.count_qubits()==n,"Ideal model has ancillary qubits");
            const size_t expected=(rounds+1)*m+(opposite ? (rounds-1)*m:0);
            require(circuit.count_detectors()==expected,"Wrong detector/closure count");
            std::mt19937_64 rng(923000+d+rounds);
            const auto [detectors,observables]=stim::sample_batch_detection_events<64>(circuit,4096,rng);
            require(!detectors.data.not_zero() && !observables.data.not_zero(),"Noiseless events");
            require(stim::circuit_to_dem(circuit,{}).count_errors()==0,"Noiseless DEM errors");
            std::set<uint64_t> ids;for (size_t id=0;id<circuit.count_detectors();id++) ids.insert(id);
            for (const auto& [id,c]:circuit.get_detector_coordinates(ids))
            {
                require(c.size()==6 && c[5]==0,"Missing metadata or unexpected flag");
                require(c[3]>=0 && c[3]<=2,"Invalid color");
                if (c[4]!=(memory_x ? 0:1)) require(c[2]>0 && c[2]<rounds,"Opposite boundary detector");
            }
        }
        const auto cap=octagons ? cc488_code_capacity(d,0,memory_x,opposite):cc666_code_capacity(d,0,memory_x,opposite);
        require(cap==generate(octagons,d,1,0,0,memory_x,opposite),"Capacity is not exactly one ideal round");
    }
}
void channels_and_provenance()
{
    for (bool octagons:{false,true}) for (bool memory_x:{false,true})
    {
        const int d=5;const size_t n=octagons ? 17:19,m=(n-1)/2;
        const int basis=memory_x ? 0:1;
        for (bool measurement_only:{false,true})
        {
            const double p=measurement_only ? 0:.01,q=measurement_only ? .02:0;
            const auto circuit=generate(octagons,d,d,p,q,memory_x,false);
            const auto dem=stim::circuit_to_dem(circuit,{});
            std::set<uint64_t> ids;for (size_t id=0;id<dem.count_detectors();id++) ids.insert(id);
            const auto coords=dem.get_detector_coordinates(ids);size_t count=0;
            dem.iter_flatten_error_instructions([&](const auto& event)
            {
                std::vector<size_t> detectors;uint64_t obs=0;
                for (auto target:event.target_data)
                    if (target.is_relative_detector_id()) detectors.push_back(target.val());
                    else if (target.is_observable_id()) obs^=uint64_t(1)<<target.val();
                require(std::abs(event.arg_data[0]-(measurement_only ? q:p))<1e-12,"Wrong event probability");
                require(event.tag==(measurement_only ? "qudec_measurement":"qudec_data"),"Noise provenance tag lost");
                if (measurement_only)
                {
                    require(detectors.size()==2 && obs==0,"Measurement event not temporal pair");
                    const auto& a=coords.at(detectors[0]);const auto& b=coords.at(detectors[1]);
                    require(a[0]==b[0] && a[1]==b[1] && std::abs(a[2]-b[2])==1 && a[3]==b[3],"Measurement endpoints wrong");
                }
                else
                {
                    require(obs==1 && !detectors.empty() && detectors.size()<=3,"Data event frame/support wrong");
                    for (auto v:detectors) require(coords.at(v)[2]==coords.at(detectors[0])[2],"Data fault spans time");
                }
                for (auto v:detectors) require(coords.at(v)[4]==basis,"Wrong basis sector");
                count++;
            });
            require(count==d*(measurement_only ? m:n),"Unexpected independent event count");
            dec::concat::CircuitDecoder decoder(dem,memory_x,true,16);
            const auto inspect=decoder.inspect(std::vector<uint8_t>(dem.count_detectors()));
            require(inspect.ambiguous_boundary_events==0,"Known generator has ambiguous boundary provenance");
            for (const auto& event:inspect.model.events)
                require(event.boundary_hint==(measurement_only ? dec::concat::BoundaryHint::Temporal:dec::concat::BoundaryHint::Spatial),"Wrong explicit event provenance");
            for (size_t id=0;id<decoder.faults.size();id++)
                require(decoder.decode(decoder.fault_syndrome({id}))==decoder.fault_observable({id}),"Integrated single fault failed");
        }
        const auto cap=octagons ? cc488_code_capacity(d,.01,memory_x):cc666_code_capacity(d,.01,memory_x);
        auto dem=stim::circuit_to_dem(cap,{});dec::concat::CircuitDecoder decoder(dem,memory_x,true,16);
        require(decoder.inspect(std::vector<uint8_t>(dem.count_detectors())).ambiguous_boundary_events==0,"Capacity corner provenance unknown");
    }
    const stim::DetectorErrorModel mixed("error[qudec_data](0.01) D0\nerror[qudec_measurement](0.02) D0\ndetector(0,0,0,0,1,0) D0\nlogical_observable L0\n");
    dec::concat::CircuitDecoder decoder(mixed,false,true,16);
    const auto inspect=decoder.inspect({0});
    require(inspect.model.events[0].boundary_hint==dec::concat::BoundaryHint::Unknown,"Mixed provenance not conservative");
    require(inspect.ambiguous_boundary_events==1,"Mixed provenance count wrong");
}
void disorder()
{
    const ComponentErrorDistribution noise{.4,942};
    for (bool octagons:{false,true})
    {
        const auto clean=generate(octagons,5,3,.01,.02,false,false);
        const auto perturbed=generate(octagons,5,3,.01,.02,false,false,noise);
        require(perturbed==apply_component_error_distribution(clean,noise),"Ideal generator drops nonuniform distribution");
        require(perturbed==generate(octagons,5,3,.01,.02,false,false,noise),"Nonuniform model not reproducible");
        std::map<std::string,double> rates;std::set<double> distinct;size_t repeats=0;
        perturbed.for_each_operation([&](const stim::CircuitInstruction& op)
        {
            if (op.gate_type!=stim::GateType::MPP) return;
            std::ostringstream key;for (auto target:op.targets) key << target.data << ',';
            const auto [it,inserted]=rates.emplace(key.str(),op.args[0]);
            if (!inserted) {require(it->second==op.args[0],"Readout rate changes between rounds");repeats++;}
            distinct.insert(op.args[0]);require(op.tag=="qudec_measurement","MPP transform lost tag");
        });
        require(repeats>0 && distinct.size()>1,"No fixed per-check disorder");
        const auto dem=stim::circuit_to_dem(perturbed,{});dec::concat::CircuitDecoder decoder(dem,false,true,16);
        require(decoder.inspect(std::vector<uint8_t>(dem.count_detectors())).ambiguous_boundary_events==0,"Nonuniform rewriting lost provenance");
    }
    const auto batched=apply_component_error_distribution(stim::Circuit("MPP[qudec_measurement](0.02) X0*X1 Z0*Z1\n"),noise);
    require(batched.operations.size()==2,"MPP products not split into independent components");
}
}
int main()
{
    clean_models();channels_and_provenance();disorder();
    std::cout << "Ideal color-code models: clean circuits, spacetime DEMs, provenance, singleton decoding and fixed disorder passed\n";
}
