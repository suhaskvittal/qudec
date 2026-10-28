<!-- Author: OpenAI Codex. -->
# Code-capacity and phenomenological models

The default qudec noise model remains SI1000. Select the ideal-check models with
these mutually exclusive aliases:

```sh
qudec concat_mwpm 7 --code cc666 -cap -p 0.01
qudec concat_mwpm 7 --code cc666 --code-capacity -p 0.01
qudec concat_mwpm 7 --code cc666 -pheno -p 0.01
qudec concat_mwpm 7 --code cc666 --phenomenological -p 0.01 \
    --measurement-error-rate 0.005
```

The ideal selectors support `sc`, `cc666`, and `cc488`. Surface-code ideal
extraction uses the existing check geometry with MPP measurements.
The restriction variants require cc666 because their lifting geometry is cc666,
and require at least two syndrome rounds for their canonical time slice. They
reject code capacity and one-round phenomenological requests. ConcatMWPM accepts
either color-code family and supports both ideal models.

Code capacity uses one noisy data-error layer and exactly one perfect stabilizer
extraction round. It defaults to one round regardless of distance; an explicit
`--rounds` value other than one is rejected. Preparation, check measurements and
final data readout are perfect. Final readout supplies a redundant perfect
closure-detector layer, not a second extraction round.

Phenomenological noise uses independent X(p) and Z(p) errors on every data qubit
before each extraction round, ideal stabilizer gates, and independent check
result flips with probability q. The default is q=p and `rounds=distance`.
`--measurement-error-rate` sets q separately, and `-r`/`--rounds` overrides the
positive round count. Initialization and final data readout remain perfect.
The final readout closes last-round measurement faults, giving `rounds+1`
memory-basis detector layers for `rounds` noisy extractions. LER is for the
whole memory experiment. No circuit gate noise or flag measurements are added.

The circuits measure both X and Z stabilizers using MPP. By default only the
memory-basis check detectors are emitted. With
`--include-opposite-basis-detectors`, the other basis also has comparisons
between interior rounds; its initial nondeterministic measurement and final
readout have no opposite-basis detectors. Color-code detector coordinates preserve
`(x,y,time,color,basis,is_flag=0)` with R/G/B colors 0/1/2 and X/Z bases 0/1.

`--noise-log-stddev` and `--noise-seed` work for both models. Data channel rates
remain fixed per qubit/channel, and MPP readout rates remain fixed per check and
basis across rounds. Changing the nominal p/q rescales the same component draws.

Native data channels carry the reserved Stim tag `qudec_data`, while noisy MPP
check measurements carry `qudec_measurement`. Stim propagates these tags to DEM
errors, including after nonuniform noise rewriting. The circuit adapter retains
their spatial/temporal provenance for diagnostics; logical effects are never used
to infer physical fault locations. Untagged SI1000 circuits retain their previous
adapter convention.

## C++ generators

Declared in `src/circuit_generator.h` and implemented with the existing layout
functions in `src/color_code_circuits.cpp`:

```cpp
auto capacity = cc666_code_capacity(distance, p, memory_x);
auto pheno = cc666_phenomenological(distance, rounds, p, q, memory_x);
auto capacity488 = cc488_code_capacity(distance, p, memory_x);
auto pheno488 = cc488_phenomenological(distance, rounds, p, q, memory_x);
```

All four accept the optional opposite-basis detector flag and component error
distribution after the required arguments. Nominal data/readout probabilities
must be finite in `[0,0.5)`; distance must be odd and at least three, and pheno
rounds must be positive. Zero noise is supported.

Surface-code C++ callers use the shared policy interface:

```cpp
auto surface_capacity = sc_circuit(distance, 1, CodeCapacityNoisePolicy(p), memory_x);
auto surface_pheno = sc_circuit(distance, rounds, PhenomenologicalNoisePolicy(p, q), memory_x);
```

Surface-code detector coordinates remain `(x,y,time)`; distance is at least two.
