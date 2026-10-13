# Libra decoder and corrected rare event sampling

The normal `qudeclib` library includes `dec::Libra`, implemented in
`src/decoder/libra.cpp`. Its auxiliary assignment uses PyMatching
and exact subset dynamic programming for components with at most eight
selected CSS edges. The decoder retains correlated matching, a 60-member
ensemble, a 20 dB gap gate, degeneracy scoring with up to 30 positive cycles
per class, and the paper's unconstrained correlated fallback.

## Build

Use a C++23 compiler. All local build output can stay in the ignored artifact directory:

```sh
cmake -S . -B .agentic_artifacts/build -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=ON
cmake --build .agentic_artifacts/build -j16
ctest --test-dir .agentic_artifacts/build --output-on-failure
```

Use `-DENABLE_MPI=OFF` for a serial build. `BUILD_TESTING=OFF` omits the test
executables. The `libra/` CMake entry point forwards to the normal project;
its binaries are in the build's `qudec/` subdirectory. The `libra_decoder`
library alias refers to `qudeclib`, which owns the single Libra implementation.

The Libra implementation and integration were contributed by OpenAI Codex.

## Monte Carlo through the normal CLI

```sh
mpirun -np 64 .agentic_artifacts/build/qudec libra 9 -p 1e-3 -r 9 \
    -m monte_carlo --seed 20261001 --libra-seed 20261001 \
    --mc-max-samples 10000000 --mc-stop-at-errors 25 --print-progress
```

Selecting `libra` enables opposing basis detectors before generating the
SI1000 circuit. The CLI uses Z memory and defaults to rounds equal to distance.
Other decoder commands can request the same detectors with
`--include-opposite-basis-detectors`. Comparing decoders requires the same
circuit, detector model, and sampled shots; the shared-shot evaluator below
provides this directly.

Options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--libra-ensemble-size` | 60 | Number of perturbed matchers |
| `--libra-gap-threshold-db` | 20 | Run the ensemble below this initial gap |
| `--libra-seed` | `0x4c49425241` | Shared perturbed priors across MPI ranks |
| `--seed` | 0 | Sampling seed; MPI ranks have distinct streams |
| `--libra-profile` | off | Collect decoding phase times |
| `--libra-validate` | off | Reconstruct and check candidate syndromes |

CLI seeds are nonnegative decimal integers. Each MPI rank owns its own
mutable decoder. Decoder statistics are summed across ranks. Timings describe
complete histories, and preparation happens before sampling. An empty syndrome
returns the identity without computing complementary class scores.

## Corrected count-stratified rare event sampling

The normal `rare_event` method now works with every decoder accepted by
`qudec`, including Libra:

```sh
mpirun -np 64 .agentic_artifacts/build/qudec libra 9 -p 1e-3 -r 9 \
    -m rare_event --seed 20261001 --libra-seed 20261001 \
    --rare-samples-per-level 10000 --rare-max-errors-per-level 25 \
    --rare-max-level 128 --rare-tail-tolerance 1e-12 --print-progress
```

`ConditionalErrorSampler` in `src/rare_event.h` samples the heterogeneous
independent DEM Bernoulli distribution conditioned on exactly K active error
instructions. The public `generate_syndromes_with_k_errors()` uses this sampler
as well. The previous uniform-subset generator was incorrect for unequal
channel probabilities.

The estimator sums `P(K) * P(decoder failure | K)` from K=0 through the
exclusive `--rare-max-level` limit, clipped to the DEM's possible count range.
It never drops low-count levels or prunes levels based on an observed error
rate. A Chernoff bound on omitted count probability also bounds its omitted
logical-error contribution. If that bound exceeds `--rare-tail-tolerance`,
the command fails and asks for a larger count limit. The reported tail bound
does not describe Monte Carlo uncertainty within sampled levels.

Sample and error budgets are global across MPI ranks. Error stopping occurs
after synchronized batches and can overshoot the error target; the sample
limit is respected exactly. Legacy configurations that skip levels or start
above K=0 are rejected. `LibraRareEventSampler` remains a compatibility alias.
IEEE floating point behavior is retained for Libra scores and sampler bounds.

## Corrected splitting evaluator

For very small logical error rates, the existing corrected splitting method
remains available as a normal build target:

```sh
.agentic_artifacts/build/libra_splitting_evaluate \
    13 64 3000000 20261001 .agentic_artifacts/splitting_d13.txt \
    0.8 8 1000 512 4
```

Arguments are distance, worker threads (at most 64), tilted base-shot cap,
seed, output path, initial odds tilt, bridge count, burn-in, records per
worker per stage, and proposal stride. It uses SI1000 p=1e-3, Z memory,
rounds=distance, and opposing basis detectors. Run one process: this executable
uses threads internally. It enables the validated local-swap kernel,
systematic resampling, and fast failure-union callback used for the reported
distance 13–17 results. It evaluates Libra and ordinary PyMatching together.
Its sampling records are correlated, and its diagnostics should be retained.

The shared-shot Monte Carlo and latency evaluator is also built normally:

```sh
.agentic_artifacts/build/libra_evaluate mc \
    9 64 10000000 20261001 .agentic_artifacts/libra_mc_d9.txt --max-errors 20
```

## C++ API

```cpp
#include "decoder/libra.h"

dec::Libra::Config config;
config.seed = 20261001;
dec::Libra decoder(dem, config);
decoder.prepare_ensemble();
auto result = decoder.decode(syndrome, observed_logical_bits);
```

`decode()` ignores observed truth when predicting. `flipped_obs` contains the
correction. `matching_data.gap` is the unsigned initial complementary gap in
dB; `matching_data.probability` is an approximate, unnormalized class score,
not a calibrated confidence. `decode_detailed()` provides error-channel IDs,
double class weights and log scores, and branch indicators. `syndrome_of()`
reconstructs the returned assignment.

Supported models have one logical observable, at most two graphlike CSS
components per channel, boundary observable crossings, and an independent
single-edge explanation for each CSS edge. Unsupported models are rejected.
The current toric generator has two observables and is outside this support.

Historical measurements remain in `Libra-PyMatching.txt` and `Libra.txt`.
