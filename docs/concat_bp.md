Author: OpenAI Codex

# Auxiliary-detector BP preprocessing for concatenated MWPM

Enable `--concat-bp` with `concat_mwpm` in either `qudec` or `qudec_chains`:

```
qudec concat_mwpm 7 --code cc666 -p 0.001 --concat-ensemble 60 --concat-bp \
    --include-opposite-basis-detectors
qudec_chains concat_mwpm -d 7 -k 2 -n 100 --concat-bp --debug
```

This optional mode replaces flag representative selection and flag
renormalization. Existing defaults are unchanged. Combining it with
`--concat-flag-renormalization` is an error. The library option is
`ConcatMWPM::config_type::belief_propagation`; direct `CircuitDecoder` users
can pass its final `belief_propagation` argument. BP requires `use_flags=true`.

The inference variables are merged original DEM faults, preserving their memory
check support, auxiliary signature and logical frame before marginalization.
Auxiliary detectors are **all flags plus available opposite-basis checks**. Both
quiet and fired outcomes condition the event probabilities. A fault with no
memory check support still participates in BP when it flips auxiliary detectors.
Flags remain colorless. Matching enforces only memory-basis check constraints;
these are not BP factors, avoiding conditioning twice on the same parity evidence.

Generate opposite-basis interior checks with
`--include-opposite-basis-detectors`. `--concat-bp` automatically uses them when
present in the DEM; it also retains both X- and Z-extraction flags. Without
opposite-basis checks it reduces to flag-only inference. The legacy non-BP flag
protocol retains its opposite-basis-only extraction flag filter. The context key
includes every auxiliary detector, so an opposite-check change with unchanged
flags also updates both matching graphs.

When each variable touches at most one auxiliary detector, the factor graph is
a set of independent stars. Each event's two posterior probabilities (quiet/fired
auxiliary detector) are precomputed once using sum-product messages. Per-shot inference selects the
appropriate probability. This gives exact *individual* marginals under the
independent DEM-event model; events become correlated after conditioning.
For a variable touching multiple auxiliary detectors, a generic sum-product
fallback runs at most 30 iterations with 0.5 old-message damping, stopping when the maximum marginal
change is below 1e-10. Such posteriors are approximate and convergence is not
guaranteed. Unlike a BP hard-decision decoder, satisfying the auxiliary syndrome is
not an early-stop criterion.

For each unflagged check support, posterior probabilities are summed separately
by logical frame. The frame with largest mass is retained (smaller frame mask
wins exact ties); that mass, clipped to [1e-12, 1-1e-9], is the hyperedge
probability. This is an approximation, not an exact odd-parity or joint class
probability. Competing logical frames are discarded before matching.
Hyperedges remain undecomposed at the input to concatenated matching.

The posterior-weighted model updates **both matching stages**, including
projected first-stage edges and second-stage event edges. The existing projection
aggregation is reused, although posterior event marginals are not independent.
Ensembling, boundary-gap, gap gating and grouped MPI use the same posterior
model. Final candidate selection minimizes the sum of -log(posterior model
probability) over correction events, using unperturbed probabilities. No p^m
penalty or unexplained-flag cost is added. Debug output identifies posterior
mode and includes the normal first/second-stage traces.

The optimization precomputes inference only. The existing full model/graph update
path is still used, so this feature does not claim to resolve the graph-update
latency bottleneck. Tests cover exact small flag posteriors, a multi-flag tree,
both-basis flag retention and legacy filtering, opposite-check-only context
invalidation, noiseless decoding, original-signature posterior model agreement,
independent star parity formulas, both-stage matching agreement, and CLI validation.

## Candidate-frame confidence gap

Contributor: OpenAI Codex.

With `--concat-bp`, decoding also reports the best alternative-frame margin from
candidates actually evaluated, without extra matching calls. For selected weight
W and the smallest unperturbed weight W_alt having a different observable frame,
the unsigned gap is `(W_alt-W)*10/ln(10)` dB (numerical ties clamp to zero).
This is a sampled-candidate likelihood-ratio proxy, not a calibrated posterior or
an exhaustive complementary matching. It does not change candidate selection.

`Result::candidate_gap_db` and `Result::has_candidate_gap` expose this value.
When no opposing-frame candidate is found, availability is false; no infinity or
artificial large confidence is inserted into the histograms. `print_stats` emits
available/unavailable counts and three distributions:

- `CONCAT_BP_UNSIGNED_GAP_DB`: all available magnitudes;
- `CONCAT_BP_SIGNED_GAP_DB`: positive for correct predictions, negative for
  incorrect predictions, using the supplied true observable only for statistics;
- `CONCAT_BP_UNSIGNED_GAP_DB_ERRORS`: magnitudes on logical failures.

Signed bins have width2 dB on [-256,256); unsigned bins width1 dB on [0,256).
Underflow/overflow counts retain samples outside these ranges. Exact ties are
zero regardless of correctness; the error histogram distinguishes their outcomes.
No-alternative shots are excluded from all three distributions. If truth is
absent, only the unsigned distribution is updated. Existing sampler callers
supply truth. `matching_data.gap` follows the signed convention when available.

The existing first-stage cluster gap and `--concat-gap-gating` are unchanged and
remain distinct from this final-correction gap. When gating skips ensembles,
only candidates actually evaluated contribute. Grouped MPI computes the margin
from gathered candidates and records statistics on the group leader; histogram
and availability reductions are included in `mpi_accumulate`.
