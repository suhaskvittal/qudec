# Libra

This adds the matching synthesis decoder from Sections III–V of
[arXiv:2408.12135](https://arxiv.org/abs/2408.12135). The default configuration
uses 60 correlated matching ensemble members, a 20 dB complementary gap gate,
and degeneracy estimation from up to 30 small positive cycles.

## Build and test

The separate build entry point includes the existing project without changing
its build files. All generated output can stay in the ignored artifact directory.

```sh
cmake -S libra -B .agentic_artifacts/libra/build \
  -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=OFF
cmake --build .agentic_artifacts/libra/build \
  --target libra_evaluate libra_splitting_evaluate libra_splitting_evaluate_fast \
           libra_splitting_evaluate_local libra_splitting_evaluate_systematic \
           libra_test libra_rare_event_test \
           libra_splitting_test libra_local_splitting_test -j8
ctest --test-dir .agentic_artifacts/libra/build --output-on-failure
```

The project requires a compiler supporting C++23, as in its existing build.
See `Libra.txt` for the compiler used for the reported measurements.

## Decoder API

```cpp
#include "decoder/libra.h"

dec::Libra::Config config;
config.seed = 20261001;
dec::Libra decoder(dem, config);
decoder.prepare_ensemble(); // Optional: move lazy setup out of decode timing.
auto result = decoder.decode(syndrome, observed_logical_bits);
```

`decode` ignores the supplied observed logical bits. They are only available to
the evaluator for checking failures. Use `decode_detailed` for the selected
hyperedge IDs, class weights, approximate log masses, and gap. `syndrome_of`
reconstructs the detector and observable parity of those hyperedges.
The log masses and `matching_data.probability` are scores without the common
Bernoulli normalization factor, rather than calibrated failure probabilities.

An instance owns mutable matching workspaces. Create one instance per worker.
Supported models have one observable, at most two graphlike CSS components per
physical channel, observable crossings on boundary edges, and independent
single component channels. Unsupported models are rejected. The zero syndrome
returns the identity correction. This is a heuristic decoder, not an exact
maximum likelihood solver.
For a complementary gap of at least 20 dB, the prediction is the original
unconstrained correlated matcher's prediction, as in the paper. The detailed
result identifies whether this fallback or the ensemble was used. The empty
syndrome shortcut leaves the complementary class scores uncomputed.

## Evaluations

Every evaluator uses Z memory, `rounds = distance`, SI1000 at `p = 1e-3`, and
opposing basis detectors. Both decoders see the same sampled syndromes and
observable truth. PyMatching is the existing ordinary, uncorrelated
`dec::PyMatching(dem, false)` baseline.

```sh
# Stop when both decoders have at least 20 failures, or at the shot cap.
.agentic_artifacts/libra/build/libra_evaluate mc \
  9 64 10000000 20261001 .agentic_artifacts/libra/example_mc.txt \
  --max-errors 20

# Timing excludes circuit generation, sampling, and decoder construction.
taskset -c 0 .agentic_artifacts/libra/build/libra_evaluate latency \
  9 1 65536 20261001 .agentic_artifacts/libra/example_latency.txt

# Failure conditioned splitting: distance, workers, base shot cap, seed,
# output, initial odds tilt, bridge levels, burn-in, records, stride.
.agentic_artifacts/libra/build/libra_splitting_evaluate_systematic \
  13 64 3000000 20261001 .agentic_artifacts/libra/example_splitting.txt \
  0.8 8 1000 512 4
```

The splitting evaluator stops its initial tilted Monte Carlo after both
decoders reach 25 failures. It then samples the union of decoder failures with
an exact Metropolis kernel and uses the known Bernoulli likelihood ratios to
return to the original distribution. Decoders retain the original `p = 1e-3`
prior throughout. Stage resampling uses ratio weighted reservoirs. Results
include acceptance, resampling diversity, and endpoint mixing diagnostics;
finite chains and populations still introduce uncertainty.
The fast target tests PyMatching first during Metropolis proposals and skips
Libra when union membership is already known. Every recorded state is decoded
with both decoders. The ordinary target provides the full proposal callback
for comparison.
The local target also mixes in probability-correct swaps between error
instructions sharing a detector. Those moves improve local fault mixing while
preserving the same conditional distribution and likelihood ratios.
The systematic target also resamples with one random CDF offset and evenly
spaced draws over the weighted reservoir slots, then shuffles worker assignment.
This reduces finite population drift. Sampling streams use distinct phase and
worker keys with a bijective SplitMix64 seed mapping.

`LibraRareEventSampler` provides an exact alternative to the original sampler's
fixed count subset generator. It samples heterogeneous Bernoulli channels
conditioned on their total count using a Poisson binomial convolution tree.
The `rare` evaluator estimates a single conditional level, which by itself is
not a full logical error rate. `rare_mix` combines levels with an odds tilted
proposal and reports the omitted count tail bound:

```sh
.agentic_artifacts/libra/build/libra_evaluate rare_mix \
  13 64 1000000 20261001 .agentic_artifacts/libra/example_count_mix.txt \
  --theta 0.8 --max-count 256 --max-errors 25
```

The original rare event sampler is unchanged. Its uniform sampling of fixed
size subsets is incorrect for unequal channel probabilities. `Libra.txt`
records that issue and the evaluation results.
