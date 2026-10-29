<!-- Author: OpenAI Codex. -->
# Decoder diagnostics

`dec::Result` has a default-empty `std::string debug_info`. Its contents are
owned by each decoder and intended for human inspection, not a stable parsing
schema. New decoders can populate it without changing sampler interfaces.

The color-code wrapper `ConcatMWPM` accepts `config_type.debug = true`.
Debugging defaults to false. Decoding does not depend on supplied truth.

```cpp
dec::ConcatMWPM decoder(dem, {.ensemble_size=20, .debug=true});
auto result = decoder.decode(syndrome, unused_truth);
std::cout << result.debug_info;
```

`qudec` exposes `--debug` for `concat_mwpm`. Unsupported decoders reject this
option. Monte Carlo and rare-event sampling print every nonempty returned
diagnostic string.

```sh
qudec concat_mwpm -d 7 -pheno --mc-max-samples 100 --debug
qudec concat_mwpm --code cc488 -d 5 --concat-ensemble 3 --debug
```

## Concatenated MWPM

The trace records actual matching calls rather than rerunning matching for
presentation. It reports:

- Each first-stage attempt on `!R` (GB), `!G` (BR) and `!B` (RG), with its
  ensemble member.
- Matched endpoints and ordered first-stage paths, with edge probabilities and
  natural-log weights `W=-log(p)`. `base_p` is the flag-conditioned projected
  probability before ensemble perturbation; `p` is the value used by that member.
- Virtual detector definitions and second-stage input. `Vj` represents projected
  first-stage edge `Pj`; these IDs are local to the ignored color. Boundary-only
  projected edges do not create fired virtual detectors, as in normal decoding.
- The second-stage matched endpoints and ordered paths, model event IDs,
  candidate correction weights and logical predictions, and the final winner.
- First-stage solutions coalesced before stage two, and infeasible candidates.

Real vertices use original DEM detector IDs with color, detector time and basis,
for example `D85[R,t=4,Z]`. These global IDs agree with the chain event listings
and the ASCII legend. `BoundarySink` denotes PyMatching's collapsed boundary.

`Cj` denotes an equivalence class keyed **only by the unflagged check-detector
support**, using the adapter's existing memory-basis projection and full detector
time IDs. Logical effects are retained on individual variants and do not split
class membership. Each variant stores `(flag IDs, probability, observable mask)`.
The class flag list is the union of its variants' flag sets. Identical variants
with the same logical effect are combined by independent XOR probability.
`CircuitDecoder::equivalence_classes()` exposes the constructed records without
running matching or selecting a logical representative.

Before projecting away opposite-basis checks, construction rejects two
positive-probability original DEM instructions with identical **full** detector
support (checks in both bases plus flags) but different logical frame masks. The
exception lists their original event IDs, conflicting masks and detector support;
such a pair implies circuit-fault distance at most two. Zero-probability channels
are excluded. Flag-only channels remain outside the check classes and supply
isolated-flag measurement priors.

For each observed flag set F, the adapter selects one actual member per class:

- If F has no overlap with the union of a class's flags, only unflagged members
  are eligible. If no unflagged representative exists, omit the class from the
  matching model. This also applies when F is empty.
- Otherwise, select the member whose expected flag set G minimizes
  `|G XOR F|`. For Z memory, F and member flag sets retain only flags from
  X-check extraction; for X memory, only flags from Z-check extraction. These
  detect the hook-error type that anticommutes with the checks being decoded.
  F includes all retained observed flags, not just those appearing in that class.
  Partial flag matches remain candidates. Unrelated observed flags do not enable
  a class whose members are all flagged and have no overlap with F.
- Equal-distance ties favor higher member probability. Exact probability ties
  retain the first member in the deterministic stored variant order.

The selected frame is retained. Precomputed flag marginals include **every**
original DEM fault that flips the flag, including data-correlated errors:
`q_f = (1 - product_e(1-2*p_e))/2`. These are separate from isolated flag-only
measurement priors. For selected member probability p and matching observed
flags H = G intersect F, the model uses `p_cond = p / product_(f in H) q_f`.
Unflagged and disjoint members stay unchanged: unrelated observed flags are
assumed independent, so their numerator and denominator factors cancel. Empty F
retains the unflagged graph. `use_flags=false` treats F as empty.

This is the requested rare-fault approximation, ignoring additional faults
cancelling a produced flag. Multiple matched flags use a product of marginals;
correlated flags can invalidate that approximation. An estimate materially above
one raises an explicit error instead of silently becoming a non-probability.
Floating-point endpoint overshoots within 1e-12 are rounded to one. The selected
raw probability and denominator are available in `inspect().weights`.

Both matching stages and final candidate scores use `-log(p_cond)`. To get that
cost from PyMatching's DEM log-odds conversion, the adapter supplies
`r=p_cond/(1+p_cond)`, yielding `log((1-r)/r)=-log(p_cond)`. Thus probabilities above
one half are supported directly, and probability one has zero cost. First-stage
parallel contributions from different classes are still combined by XOR as
before. Representative selection itself remains minimum symmetric difference.

`Ej` denotes an event in the effective memory-sector decoding model, **not an
original DEM instruction ID**. E IDs appear in second-stage paths and final
corrections. Debug output reports matching inputs, paths, edge probabilities and
weights, virtual detector definitions, candidate corrections and the winner.
It omits flag priors, observed/expected flag sets, incidence neighborhoods and
per-class representative tables. Class and representative data remain available
through `equivalence_classes()` and `inspect()` for programmatic inspection.

Printed matching weights use natural-log units. Projected parallel contributions
can be merged in graph construction, so a projected edge probability can differ
from an individual selected member's probability. Debugging does not change
representative selection or decoding.

Final candidate selection remains minimum representative-model correction weight, with the
existing tolerance and candidate order for ties. Debug output includes all
ensemble attempts and can be large. Graph construction and candidate generation
are unchanged. Trace generation adds path reconstruction, weight auditing and
formatting to measured decode latency, so debug latency is not a benchmark of
normal decoding.
