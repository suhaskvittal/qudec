// Author: OpenAI Codex. Shared generator policy and extension regressions.
#include "circuit_generator.h"
#include <stim/simulators/frame_simulator_util.h>
#include <stim/util_top/circuit_to_dem.h>
#include <iostream>
#include <set>
#include <stdexcept>

namespace
{
void require(bool value,const char* message)
{ if (!value) throw std::runtime_error(message); }

// A user policy can inspect/replace operation implementations without changing
// the geometry generator or registering a new error-model enum in it.
class RecordingPolicy : public CircuitNoisePolicy
{
public:
    explicit RecordingPolicy(const CircuitNoisePolicy& base) : base(base) {}
    CheckImplementation check_implementation() const override { return base.check_implementation(); }
    void validate(uint32_t rounds) const override { base.validate(rounds); }
    void append(stim::Circuit& c,const CircuitOperation& op) const override
    { seen.insert(op.kind);base.append(c,op); }
    const CircuitNoisePolicy& base;
    mutable std::set<CircuitOperationKind> seen;
};
void adapters()
{
    const ComponentErrorDistribution noise{.3,42};
    for (bool octagons : {false,true}) for (bool x : {false,true}) for (bool opposite : {false,true})
    {
        auto generate=[&](const CircuitNoisePolicy& policy,uint32_t rounds)
        {
            return octagons ? cc488_circuit(5,rounds,policy,x,opposite,noise)
                            : cc666_circuit(5,rounds,policy,x,opposite,noise);
        };
        const SI1000NoisePolicy si(.001);
        const PhenomenologicalNoisePolicy pheno(.01,.02);
        const CodeCapacityNoisePolicy cap(.01);
        require(generate(si,3)==(octagons ? cc488_si1000(5,3,.001,x,opposite,noise)
            : cc666_si1000(5,3,.001,x,opposite,noise)),"SI1000 adapter mismatch");
        require(generate(pheno,3)==(octagons ? cc488_phenomenological(5,3,.01,.02,x,opposite,noise)
            : cc666_phenomenological(5,3,.01,.02,x,opposite,noise)),"Pheno adapter mismatch");
        require(generate(cap,1)==(octagons ? cc488_code_capacity(5,.01,x,opposite,noise)
            : cc666_code_capacity(5,.01,x,opposite,noise)),"Capacity adapter mismatch");
        RecordingPolicy custom(si);
        require(generate(custom,3)==generate(si,3),"Custom gate policy bypassed");
        require(custom.seen.contains(CircuitOperationKind::Reset) &&
            custom.seen.contains(CircuitOperationKind::Unitary2) &&
            custom.seen.contains(CircuitOperationKind::Measure) &&
            custom.seen.contains(CircuitOperationKind::MeasureReset) &&
            custom.seen.contains(CircuitOperationKind::MeasurementIdle),"Gate hooks missing");
        RecordingPolicy product(pheno);
        require(generate(product,3)==generate(pheno,3),"Custom product policy bypassed");
        require(product.seen.contains(CircuitOperationKind::RoundStart) &&
            product.seen.contains(CircuitOperationKind::ProductMeasurement),"Product hooks missing");
        bool rejected=false;
        try { generate(cap,3); } catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"Capacity accepted multiple rounds");
    }
    for (bool x : {false,true}) for (bool opposite : {false,true})
    {
        SI1000NoisePolicy si(.001);RecordingPolicy custom(si);
        require(sc_circuit(3,3,custom,x,opposite,noise)==sc_si1000(3,3,.001,x,opposite,noise),"Surface policy mismatch");
        require(toric_circuit(3,3,custom,x,opposite,noise)==toric_si1000(3,3,.001,x,opposite,noise),"Toric policy mismatch");
        require(custom.seen.contains(CircuitOperationKind::MeasureReset),"Surface/toric hook missing");
    }
}
void surface_ideal()
{
    for (uint32_t d : {3,5}) for (bool x : {false,true}) for (bool opposite : {false,true})
        for (bool capacity : {false,true})
        {
            const uint32_t rounds = capacity ? 1 : d;
            const PhenomenologicalNoisePolicy clean(0,0), noisy(.01,.02);
            const CodeCapacityNoisePolicy cap_clean(0), cap_noisy(.01);
            const CircuitNoisePolicy& zero = capacity ? static_cast<const CircuitNoisePolicy&>(cap_clean) : clean;
            const CircuitNoisePolicy& noise = capacity ? static_cast<const CircuitNoisePolicy&>(cap_noisy) : noisy;
            const auto circuit = sc_circuit(d,rounds,zero,x,opposite);
            auto rng = std::mt19937_64(42);
            const auto detections = stim::sample_batch_detection_events<64>(circuit,128,rng);
            require(!detections.first.data.not_zero() && !detections.second.data.not_zero(),"Ideal surface circuit is not deterministic");
            const auto dem = stim::circuit_to_dem(sc_circuit(d,rounds,noise,x,opposite),{.decompose_errors=true});
            require(dem.count_errors()>0 && dem.count_observables()==1,"Ideal surface DEM missing faults or observable");
            if (capacity) require(dem.count_errors()<=d*d,"Surface capacity has more channels than data faults");
        }

}
}
int main()
{
    adapters();
    surface_ideal();
    std::cout << "Shared operation policies: code generators, legacy adapters and custom hooks passed\n";
}
