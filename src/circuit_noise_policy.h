// Author: OpenAI Codex. Operation implementations shared by circuit generators.
#ifndef CIRCUIT_NOISE_POLICY_H
#define CIRCUIT_NOISE_POLICY_H

#include <stim.h>
#include <cstdint>
#include <string_view>
#include <vector>

// The generator owns geometry, scheduling and detector records. A policy owns
// the implementation of each scheduled operation, including its noise channels.
// Product extraction is ideal MPP; Bell extraction uses the flagged gate schedule.
enum class CheckImplementation { GateSequence, Product };
enum class CircuitOperationKind {
    Reset, Unitary1, Unitary2, Measure, MeasureReset,
    GateIdle, MeasurementIdle, PostMeasurement, RoundStart, ProductMeasurement
};
struct CircuitOperation {
    CircuitOperationKind kind;
    const std::vector<uint32_t>& targets;
    std::string_view gate{};
    char basis{'Z'};
    // Retain historical surface/toric zero-probability noise instructions.
    bool emit_zero_noise{false};
};

class CircuitNoisePolicy {
public:
    virtual ~CircuitNoisePolicy() = default;
    virtual CheckImplementation check_implementation() const = 0;
    virtual void validate(uint32_t rounds) const = 0;
    virtual void append(stim::Circuit& circuit, const CircuitOperation& operation) const = 0;
};

class SI1000NoisePolicy : public CircuitNoisePolicy {
public:
    explicit SI1000NoisePolicy(double physical_error_rate) : p(physical_error_rate) {}
    CheckImplementation check_implementation() const override { return CheckImplementation::GateSequence; }
    void validate(uint32_t rounds) const override;
    void append(stim::Circuit&, const CircuitOperation&) const override;
    double p;
};

class PhenomenologicalNoisePolicy : public CircuitNoisePolicy {
public:
    PhenomenologicalNoisePolicy(double data_error_rate, double measurement_error_rate)
        : data_p(data_error_rate), measurement_p(measurement_error_rate) {}
    CheckImplementation check_implementation() const override { return CheckImplementation::Product; }
    void validate(uint32_t rounds) const override;
    void append(stim::Circuit&, const CircuitOperation&) const override;
    double data_p;
    double measurement_p;
};

class CodeCapacityNoisePolicy : public PhenomenologicalNoisePolicy {
public:
    explicit CodeCapacityNoisePolicy(double p) : PhenomenologicalNoisePolicy(p, 0) {}
    void validate(uint32_t rounds) const override;
};
#endif
