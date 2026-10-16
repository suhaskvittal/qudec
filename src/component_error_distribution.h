// Implementation authored by OpenAI Codex.
#ifndef COMPONENT_ERROR_DISTRIBUTION_H
#define COMPONENT_ERROR_DISTRIBUTION_H

#include <stim.h>
#include <cstdint>

// Fixed device disorder, sampled once when constructing a circuit. A component
// is a noise-channel type, ordered qubit target(s), and tag.
// Repeated applications of the same component keep the same sampled rate.
struct ComponentErrorDistribution {
    double log_stddev = 0;
    uint64_t seed = 0;

    void validate() const;
    // Draw min(1, nominal_probability * exp(sigma*Z - sigma^2/2)).
    // `component_id` identifies the hardware component independently of time.
    double probability(double nominal_probability, uint64_t component_id) const;
};

// Supports the scalar Pauli/depolarizing channels used by SI1000. Preserves
// gates, annotations, tags and REPEAT structure; refuses unsupported noise.
// Depolarizing rates saturate at 3/4 (one qubit) or 15/16 (two qubits), the
// maximally mixed channels, because Stim DEM conversion rejects over-mixing.
// At zero spread this returns an exact circuit copy without splitting targets.
stim::Circuit apply_component_error_distribution(
    const stim::Circuit& circuit, const ComponentErrorDistribution& distribution);

#endif
