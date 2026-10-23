<!-- Author: OpenAI Codex, 22 October 2026. -->
# Correlated restriction decoder

`dec::CorrelatedRestrictionDecoder` implements one synchronous correlation
update for triangular cc666. It uses vendored PyMatching for the initial RG,
RB, and BG matchings, discounts parent hyperedge weights, and repeats all
three matchings once. The updated paths then undergo connected-component
construction, canonical flattening, and local lifting. It does not iterate
or fall back to the uncorrelated solution.

For a parent event with k >= 3 non-flag check detectors, let
W = log((1-p)/p). When updating a target graph, count r distinct detector
pairs selected in the OTHER initial restricted matchings. Set
W' = (1-r/C(k,2))*W. For a three-detector event, one supported pair gives
2W/3, and two give W/3. One- and two-detector events retain their weights.
Each update starts from the original event weight. The target graph's own
matching never supplies evidence. Parallel event contributions are combined
using independent-parity probabilities, as in the validated experiment.

## CLI names

Use `qudec correlated_restriction ...` for the fractional-weight single-pass
variant. Use `qudec restriction ...` for the uncorrelated restricted-lattice
matching baseline. They share the same flag-aware Bell-circuit model, connected
components, flattening, and lifting. The baseline is useful for paired
comparisons with the single-pass rule.

## Circuit-level use

```sh
qudec correlated_restriction 7 -p 0.0005 --mc-max-samples 10000 \
    --mc-stop-at-errors 100000 --seed 810007
```

This selects cc666 SI1000, memory Z, with distance-many extraction rounds by
default. The CLI preserves the original undecomposed DEM. The sampler's seed
and batching determine its samples; the command is not a replay of a previous
experimental batch. The decoder reports its logical prediction and timing
through the normal sampler interface and supports the existing rare-event
sampler. Each worker needs its own decoder instance.

The C++ DEM constructor takes `config_type` with distance, physical_p,
memory_x (default false), and use_flags (default true). It supports one logical
observable and requires detector coordinates
`(x,y,time,color,basis,is_flag)` from the color-code circuit generator and at
least two extraction rounds. Checks have colors 0/1/2; flags are colorless,
and flags from both sectors condition the memory-sector error classes.
Isolated flag-only DEM events determine per-flag measurement priors. The
initial weights include flagged/unflagged mixing and penalties for unrelated
observed flags. Correlation discounts those conditioned weights. The graph
must be graphlike after each color restriction. Spatial canonical detectors
come from the second extraction round, avoiding the initial time boundary.

## Code-capacity and phenomenological use

```cpp
#include "decoder/correlated_restriction.h"

dec::CorrelatedRestrictionDecoder capacity(7, 1, 0.03);
dec::CorrelatedRestrictionDecoder pheno(7, 7, 0.01);
```

These constructors reproduce the one-sector experiments: independent data
faults at probability p, and (for phenomenological noise) measurement faults
also at p, with perfect initial/final boundaries. Capacity has one detector
layer. Phenomenological noise has d noisy rounds and d+1 detector layers.
The correction is chosen across reference-color lifts using correction weight
for capacity and the existing locally optimized feasible-history cost for
phenomenological noise. That cost is not an exact logical posterior.

The earlier `dec::RestrictionDecoder` remains available and unchanged.
The integrated correlated implementation does not depend on files or
libraries under `.agentic_artifacts/`. Its main sources are
`src/decoder/correlated_restriction.{h,cpp}` and
`src/decoder/correlated_restriction/`.
