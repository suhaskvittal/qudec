// Author: OpenAI Codex. Physical extraction and frame-parity regression tests.
#include "circuit_generator.h"
#include <stim/simulators/frame_simulator_util.h>
#include <stim/simulators/tableau_simulator.h>
#include <stim/util_top/circuit_to_dem.h>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <tuple>
#include <stdexcept>

namespace
{
void require(bool ok,const char* message) {if (!ok) throw std::runtime_error(message);}
stim::Circuit generate(bool oct,uint32_t d,uint32_t r,double p,bool x,bool opposite,int phase,bool native)
{
    return oct ? cc488_superdense(d,r,SI1000NoisePolicy(p),x,opposite,{},phase,native)
               : cc666_superdense(d,r,SI1000NoisePolicy(p),x,opposite,{},native);
}
void schedule(const stim::Circuit& c,bool native)
{
    std::set<uint32_t> active;
    c.for_each_operation([&](const auto& op)
    {
        if (op.gate_type==stim::GateType::TICK) {active.clear();return;}
        if (op.gate_type==stim::GateType::R || op.gate_type==stim::GateType::M ||
            op.gate_type==stim::GateType::H || op.gate_type==stim::GateType::CX || op.gate_type==stim::GateType::CZ)
            for (auto t:op.targets) require(active.insert(t.qubit_value()).second,"Conflicting physical moment");
        if (native) require(op.gate_type!=stim::GateType::CX,"Native circuit contains CX");
        else require(op.gate_type!=stim::GateType::CZ,"CNOT circuit contains CZ");
    });
}
void absolute_parities(const stim::Circuit& c)
{
    const auto reference=stim::TableauSimulator<64>::reference_sample_circuit(c);
    size_t measured=0;
    c.for_each_operation([&](const auto& op)
    {
        measured+=op.count_measurement_results();
        if (op.gate_type!=stim::GateType::DETECTOR && op.gate_type!=stim::GateType::OBSERVABLE_INCLUDE) return;
        bool parity=false;
        for (auto t:op.targets) parity^=reference[measured+t.rec_offset()];
        require(!parity,"Nonzero absolute noiseless detector/logical parity");
    });
}
// Optional independent validation against the downloaded source release.
// Compare gates and timings, ignoring annotation formatting and batching.
auto physical_moments(const stim::Circuit& c)
{
    using Atom=std::tuple<int,std::vector<uint32_t>>;
    std::vector<std::vector<Atom>> moments;
    std::vector<Atom> current;
    c.for_each_operation([&](const auto& op)
    {
        if (op.gate_type==stim::GateType::TICK)
        {
            std::sort(current.begin(),current.end());
            moments.push_back(std::move(current)); current.clear(); return;
        }
        if (op.gate_type!=stim::GateType::R && op.gate_type!=stim::GateType::M && op.gate_type!=stim::GateType::H &&
            op.gate_type!=stim::GateType::CX && op.gate_type!=stim::GateType::CZ) return;
        const size_t stride=(op.gate_type==stim::GateType::CX || op.gate_type==stim::GateType::CZ) ? 2:1;
        for (size_t j=0;j<op.targets.size();j+=stride)
        {
            std::vector<uint32_t> qs;
            for (size_t k=0;k<stride;k++) qs.push_back(op.targets[j+k].qubit_value());
            if (op.gate_type==stim::GateType::CZ) std::sort(qs.begin(),qs.end());
            current.emplace_back(int(op.gate_type),qs);
        }
    });
    if (!current.empty()) {std::sort(current.begin(),current.end());moments.push_back(std::move(current));}
    return moments;
}
}
int main(int argc,char** argv)
{
    std::mt19937_64 rng(418);
    require(cc666_superdense(3,1,SI1000NoisePolicy(0),false).str()==
            generate(false,3,1,0,false,false,0,false).str(),"6.6.6 API default is not CNOT");
    require(cc488_superdense(3,1,SI1000NoisePolicy(0),false).str()==
            generate(true,3,1,0,false,false,0,false).str(),"4.8.8 API default is not CNOT");
    size_t cases=0;
    for (bool oct:{false,true}) for (uint32_t d: {3u,5u,7u,9u,11u,15u})
    {
        if (oct && d%4!=3) continue;
        const size_t n=oct ? (d*d+2*d-1)/2:(3*d*d+1)/4;
        const size_t na=oct ? (d*d+3*d-4)/2:n-1;
        const size_t flags=oct ? (d-1)/2:0;
        for (bool x:{false,true}) for (bool native:{false,true})
            for (int phase=0;phase<(oct ? 2:1);phase++) for (bool opposite:{false,true})
                for (uint32_t rounds:{1u,2u,d})
                {
                    const auto c=generate(oct,d,rounds,0,x,opposite,phase,native);
                    require(c.count_qubits()==n+na,"Incorrect qubit census");
                    require(c.count_measurements()==rounds*na+n,"Incorrect measurement census");
                    const auto nf=(n-1)/2;
                    require(c.count_detectors()==nf*(rounds+1+(opposite ? rounds-1:0))+flags*rounds,"Incorrect detector census");
                    schedule(c,native);
                    std::set<uint64_t> ids;
                    for (size_t k=0;k<c.count_detectors();k++) ids.insert(k);
                    for (const auto& [id,p]:c.get_detector_coordinates(ids))
                    {
                        require(p.size()==6,"Missing six-coordinate metadata");
                        require(p[5]==1 ? p[3]==-1:(p[3]>=0 && p[3]<=2),"Incorrect detector color");
                    }
                    // DEM construction rejects nondeterministic detectors AND logicals;
                    // samples additionally exercise random opposite-basis preparation.
                    (void)stim::circuit_to_dem(c,{});
                    absolute_parities(c);
                    const auto [dets,obs]=stim::sample_batch_detection_events<64>(c,128,rng);
                    require(!dets.data.not_zero() && !obs.data.not_zero(),"Noiseless detection/logical events");
                    ++cases;
                }
        for (bool x:{false,true}) for (bool native:{false,true}) for (int phase=0;phase<(oct ? 2:1);phase++)
        {
            const auto noisy=generate(oct,d,d,.001,x,true,phase,native);
            const auto dem=stim::circuit_to_dem(noisy,{});
            require(dem.count_errors()>0,"Missing SI1000 faults");
            dem.iter_flatten_error_instructions([&](const auto& op)
            {
                bool detector=false,logical=false;
                for (auto t:op.target_data) {detector|=t.is_relative_detector_id();logical|=t.is_observable_id();}
                require(detector || !logical,"Undetected single-fault logical error");
            });
        }
        const auto two=physical_moments(generate(oct,d,2,0,false,false,0,true));
        std::vector<size_t> starts,ends;
        for (size_t m=0;m<two.size();m++) for (const auto& [gate,qs]:two[m]) if (qs[0]==n)
        {
            if (gate==int(stim::GateType::R)) starts.push_back(m);
            if (gate==int(stim::GateType::M)) ends.push_back(m);
        }
        require(starts.size()==2 && ends.size()==2,"Missing complete native cycles");
        for (size_t cycle=0;cycle<2;cycle++)
            require(ends[cycle]-starts[cycle]+1==(oct ? (d==3 ? 23u:24u):15u),"Incorrect native recurring-cycle depth");
    }
    if (argc==2)
    {
        for (uint32_t d:{3u,7u,11u,15u}) for (bool native:{false,true})
            for (bool x:{false,true}) for (int phase=0;phase<2;phase++)
            {
                const auto path=std::string(argv[1])+"/circuits/"+(native ? "native-cz/":"cnot/")+"planar_488_d"+
                    std::to_string(d)+"_r"+std::to_string(d)+"_"+(x ? "X":"Z")+"_phase"+std::to_string(phase)+".stim";
                std::ifstream input(path);
                require(bool(input),"Missing independent reference circuit");
                const std::string text((std::istreambuf_iterator<char>(input)),{});
                const auto actual=physical_moments(generate(true,d,d,0,x,true,phase,native));
                const auto expected=physical_moments(stim::Circuit(text));
                if (actual!=expected) throw std::runtime_error("Physical reference mismatch: "+path);
            }
        std::cout<<"All 32 published 4.8.8 physical gate schedules match exactly, including native moment timings\n";
    }
    for (uint32_t d:{5u,9u,13u})
    {
        bool rejected=false;
        try {(void)cc488_superdense(d,1,SI1000NoisePolicy(0),false);} catch (const std::invalid_argument&) {rejected=true;}
        require(rejected,"Unsupported 4.8.8 distance accepted");
    }
    std::cout<<"Superdense: "<<cases<<" noiseless cases, both bases/gate sets/phases, metadata and noisy DEMs passed\n";
}
