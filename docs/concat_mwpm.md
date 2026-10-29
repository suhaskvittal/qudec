<!-- Author: OpenAI Codex. 22 October 2026. -->

# Concatenated MWPM for color codes

`dec::ConcatMWPM` implements the two-stage, three-color matching decoder of
[Lee, Li, and Bartlett, Quantum 9, 1609 (2025)](https://doi.org/10.22331/q-2025-01-27-1609).
The matching stages use the vendored PyMatching C++ library.

Run a cc666 Z-memory SI1000 experiment with:

```sh
qudec concat_mwpm -d 7 -p 0.0005 -r 7 \
    --mc-max-samples 10000 --mc-stop-at-errors 10000 --seed 42
```

`concat_mwpm` selects `cc666` by default; `--code cc666` is also accepted.
Other decoder names retain their surface-code default. Both Monte Carlo
and rare-event samplers accept the new decoder. Component noise options
such as `--noise-log-stddev` and `--noise-seed` also apply to the circuit.

The CLI passes the original **undecomposed** detector error model to this
decoder. Decomposing its hyperedges before flag conditioning loses the
event relationships needed by the two-stage construction.

## C++ interface

```cpp
#include "decoder/concat_mwpm.h"
#include "circuit_generator.h"
#include <stim/util_top/circuit_to_dem.h>

auto circuit = cc666_si1000(7, 7, 0.0005, false);
auto dem = stim::circuit_to_dem(circuit, {});
dec::ConcatMWPM decoder(dem);
// result.flipped_obs contains the predicted logical frame.
auto result = decoder.decode(syndrome, actual_observable);
```

The true observable argument is unused by decoding. The sampler compares
the returned prediction with truth separately. For X memory, generate an
X-memory circuit and set `config_type::memory_x=true`. Flags are enabled
by default; `use_flags=false` marginalizes their information for controlled
comparisons. `flag_cache_capacity` bounds the number of cached flag contexts.
Each instance owns mutable matching workspaces and should have one worker.

The DEM must have one logical observable and detector coordinates
`(x,y,time,color,basis,is_flag)`. Check colors are `0/1/2` for red/green/blue;
bases are `0/1` for X/Z. Flags are colorless: their color field is ignored,
and flags from both bases affect every applicable matching projection.

The lower-level `dec::concat::Decoder`, in
`decoder/concat_mwpm/matching.h`, accepts a single-sector colored `Model`
directly. This supports code-capacity or phenomenological models without
the circuit metadata adapter. Both matching stages preserve time coordinates.

## Flag probabilities

For each flag detector, the adapter infers its readout probability from
original DEM instructions that flip exactly that flag and no other detector
or logical observable. Duplicate independent instructions combine by XOR
probability. This inference precedes check-basis marginalization, so a flag
accompanied by an opposite-basis check does not count as an isolated flip.
There is no hardcoded `5p` readout prior; an absent isolated channel has zero
probability. Correlated multiple-flag-only channels are not singleton flips.

Variants are grouped by identical memory-check support and logical frame.
With one associated observed flag, their local probability is
`(1-p_M)Y + p_M X`. Multiple associated flags use independent readout
likelihoods with the inferred probability of each flag. Flagged variants
require an associated observed flag. With no associated observed flags,
the class uses its unflagged variant.

Every observed flag outside a class's flag union also contributes its
`p_M` factor. Consequently, a wholly unflagged class is penalized by the
product of the readout probabilities of all observed flags. This fixes the
missing full-flag mismatch penalty discussed in
[Vittal et al., Equation 9](https://arxiv.org/html/2409.14283v1#S6.E9).
The local mixture and per-event concatenated matching weights are retained;
the equation's event-arity exponent is not added. This is a flag-reweighting
heuristic and does not perform exact joint inference of all flags and faults.

The conditioned probabilities feed all three restricted graphs and all
three monochromatic graphs. Each syndrome therefore receives six matching
calls. The least total second-stage log-odds weight selects the prediction.
Flags are not matched as additional check vertices. Each candidate is
checked against the complete supplied memory-check syndrome.

## Validation and limits

The artifact implementation was evaluated on code-capacity and
phenomenological models through distance 11. Under Bell-flagged SI1000
noise at `p=0.0005`, corrected flag weighting passed all single DEM events
at distances 7 and 9, 10,000 sampled connected two-event chains at distance
7, and 10,000 sampled connected three-event chains at distance 9. Chain
connectivity includes shared original flags and either-basis checks before
XOR cancellation. These sampled results do not establish universal bounds.

The algorithm can fail below half distance; a saved distance-7 code-capacity
weight-three counterexample remains uncorrectable with its deterministic
matching choices. The paper also describes asymptotic counterexamples.

The production adapter currently requires every retained event to be
graphlike in both stages for all three colors. It rejects unsupported
projections, parallel second-stage edges with different logical effects,
and single logical faults without memory-check support. No unsupported
events are silently dropped. The circuit-level validation covers cc666
distances 7 and 9; cc488 has not been validated for this adapter.

`concat_mwpm_test` covers the flag-prior inference and full-flag penalty,
colorless flags, both memory bases, graph-weight audits, cache eviction,
all single and double code-capacity errors at distance 5, and the saved
distance-7 two-hook regression. CLI tests exercise Monte Carlo sampling,
zero-noise decoding, and rare-event sampling through the executable.
Decoder statistics report full-history latency,
including graph construction when a new flag context is first observed.
