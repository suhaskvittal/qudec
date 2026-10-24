// Author: OpenAI Codex. Noise policies implement scheduled operations.
#include "circuit_noise_policy.h"
#include <cmath>
#include <stdexcept>
#include <string>

void SI1000NoisePolicy::validate(uint32_t rounds) const
{
    if (!rounds || !std::isfinite(p) || p < 0 || p > .2)
        throw std::invalid_argument("SI1000: need rounds >= 1 and finite 0 <= p <= 0.2");
}
void SI1000NoisePolicy::append(stim::Circuit& c,const CircuitOperation& op) const
{
    auto noise=[&](std::string_view gate,double rate)
    {
        if (!op.targets.empty() && (rate>0 || op.emit_zero_noise))
            c.safe_append_ua(gate,op.targets,rate);
    };
    auto flip=[&](double rate) { if (rate>0) noise(op.basis=='X' ? "Z_ERROR":"X_ERROR",rate); };
    auto gate=[&](std::string_view name) { c.safe_append_u(name,op.targets); };
    switch(op.kind)
    {
        case CircuitOperationKind::Reset:
            gate(op.basis=='X' ? "RX":"R");flip(2*p);break;
        case CircuitOperationKind::Unitary1: gate(op.gate);noise("DEPOLARIZE1",p/10);break;
        case CircuitOperationKind::Unitary2: gate(op.gate);noise("DEPOLARIZE2",p);break;
        case CircuitOperationKind::Measure:
            flip(5*p);gate(op.basis=='X' ? "MX":"M");break;
        case CircuitOperationKind::MeasureReset:
            flip(5*p);gate(op.basis=='X' ? "MRX":"MR");flip(2*p);break;
        case CircuitOperationKind::GateIdle: noise("DEPOLARIZE1",p/10);break;
        case CircuitOperationKind::MeasurementIdle: noise("DEPOLARIZE1",2*p);break;
        case CircuitOperationKind::PostMeasurement: noise("DEPOLARIZE1",p);break;
        case CircuitOperationKind::RoundStart: break;
        case CircuitOperationKind::ProductMeasurement:
            throw std::invalid_argument("SI1000 needs an explicit extraction gate schedule");
    }
}
void PhenomenologicalNoisePolicy::validate(uint32_t rounds) const
{
    if (!rounds || !std::isfinite(data_p) || !std::isfinite(measurement_p) ||
        data_p<0 || data_p>=.5 || measurement_p<0 || measurement_p>=.5)
        throw std::invalid_argument("ideal-check noise: need rounds >= 1 and finite 0 <= probabilities < 0.5");
}
void PhenomenologicalNoisePolicy::append(stim::Circuit& c,const CircuitOperation& op) const
{
    switch(op.kind)
    {
        case CircuitOperationKind::Reset: c.safe_append_u(op.basis=='X' ? "RX":"R",op.targets);break;
        case CircuitOperationKind::Unitary1:
        case CircuitOperationKind::Unitary2: c.safe_append_u(op.gate,op.targets);break;
        case CircuitOperationKind::Measure: c.safe_append_u(op.basis=='X' ? "MX":"M",op.targets);break;
        case CircuitOperationKind::MeasureReset: c.safe_append_u(op.basis=='X' ? "MRX":"MR",op.targets);break;
        case CircuitOperationKind::RoundStart:
            if (data_p>0)
            {
                c.safe_append_ua("X_ERROR",op.targets,data_p,"qudec_data");
                c.safe_append_ua("Z_ERROR",op.targets,data_p,"qudec_data");
            }
            break;
        case CircuitOperationKind::ProductMeasurement:
            c.safe_append_ua("MPP",op.targets,measurement_p,"qudec_measurement");break;
        case CircuitOperationKind::GateIdle:
        case CircuitOperationKind::MeasurementIdle:
        case CircuitOperationKind::PostMeasurement: break;
    }
}
void CodeCapacityNoisePolicy::validate(uint32_t rounds) const
{
    PhenomenologicalNoisePolicy::validate(rounds);
    if (rounds!=1) throw std::invalid_argument("code capacity requires exactly one extraction round");
}
