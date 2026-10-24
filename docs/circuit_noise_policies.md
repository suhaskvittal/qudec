<!-- Author: OpenAI Codex. -->
# Circuit operation and noise policies

Use `sc_circuit`, `toric_circuit`, `cc666_circuit`, or `cc488_circuit` from
`src/circuit_generator.h`. Each accepts distance, rounds, a `CircuitNoisePolicy`,
memory basis, the opposite-basis detector option, and optional fixed component
error-rate disorder. Geometry, schedules, record bookkeeping, detector metadata,
and observable definitions belong to the code generator.

```cpp
const PhenomenologicalNoisePolicy noise(0.01, 0.02);
auto circuit = cc666_circuit(7, 7, noise, false);
const CodeCapacityNoisePolicy capacity(0.01);
auto single_round = cc666_circuit(7, 1, capacity, false);
const SI1000NoisePolicy circuit_noise(0.0005);
auto bell_flagged = cc666_circuit(7, 7, circuit_noise, false);
```

`CircuitNoisePolicy::append` implements each scheduled `CircuitOperation` and its
noise: reset, one- and two-qubit gates, measurement/reset, measurement, gate and
measurement idles, post-measurement noise, start-of-round data errors, and product
measurement. A new model subclasses this interface; no new code generator is
needed. The interface is runtime polymorphic so implementation details stay out
of generator headers. Policies must preserve the operation's measurement count,
ordering, and ideal action; changing these requires a corresponding extraction
schedule because detector record offsets depend on them. Policies must not use
per-call randomness to assign hardware rates: a round is constructed once and
repeated. Use the existing `ComponentErrorDistribution` for fixed device disorder.

Color codes share a single geometry/round/detector implementation. The policy's
`check_implementation()` selects either the existing Bell-flagged gate schedule
or direct stabilizer-product measurements. SI1000 supplies the former; code
capacity and phenomenological policies supply the latter. These are extraction
implementations, not independent geometry generators. Colorless flag decoding,
single-round flag detectors, color/basis metadata, and optional interior
opposite-basis checks retain their existing behavior. Ideal-check models retain
`qudec_data`/`qudec_measurement` provenance tags and have no flags or ancillas.

Surface codes support both the gate schedule and product extraction, sharing
check supports and detector/readout bookkeeping. Code-capacity and
phenomenological policies select product extraction. Toric codes currently
require a gate-level (`GateSequence`) policy. Neither gate schedule adds flags.

`CodeCapacityNoisePolicy` requires exactly one round. Phenomenological extraction
accepts explicit round counts (qudec defaults to distance). SI1000 validates the
physical probability against its largest channel, measurement error `5p`.
Nonuniform component-rate rewriting is applied after the complete circuit is
constructed and retains tags and repeat blocks.

The old `*_si1000`, `*_code_capacity`, and `*_phenomenological` entry points remain
thin compatibility adapters. `qudec` selects a policy once and calls only the
model-independent code generator. Existing flags and default model selection are
unchanged. Ideal circuits now use repeat/coordinate-shift blocks like circuit-level
circuits; their expanded operations and detector error models remain equivalent.
