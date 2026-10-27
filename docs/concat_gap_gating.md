<!-- Author: OpenAI Codex. -->

# Cluster-gap-gated ConcatMWPM ensembling

```sh
qudec concat_mwpm 7 -cap -p .01 --concat-ensemble 20 --concat-gap-gating
qudec concat_mwpm 7 -cap -p .05 --concat-ensemble 20 \
    --concat-gap-gating --concat-gap-gating-threshold-db 12
```

`--concat-gap-gating` is disabled by default. When enabled, its threshold
defaults to **20 dB**. `--concat-gap-gating-threshold-db` sets a finite,
nonnegative threshold; a nondefault threshold requires gating to be enabled. The ensemble size
remains controlled by `--concat-ensemble`; with K=1, gating does no extra work
and the returned result does not claim a measured cluster gap.

The decoder first runs the three ordinary, unperturbed color candidates. It
selects the usual minimum-full-correction-weight K=1 output and evaluates
**that candidate's color's first-stage cluster gap**. If
`gap_db + 1e-9 < threshold_db`, it runs the remaining perturbed members.
Otherwise it returns the baseline correction. The tolerance preserves strict
less-than behavior for numerical ties, including a mathematically 20 dB vacuum
gap represented slightly below 20. Member-zero results and virtual-syndrome
coalescing are reused; expansion does not repeat those matching problems.

Cluster radii are collected through PyMatching's pre-extraction callback during
the three original first-stage matchings. Confidence estimation adds a weighted
boundary-path calculation, without another MWPM call. The metric uses the
existing cluster-ball coverage construction and converts natural-log weights
to dB with `10/ln(10)`. Empty matchings have their actual vacuum gap. Positive
physical corner edges between boundaries are retained, including the existing
supported-only rule under flag renormalization.

`dec::ConcatMWPM::config_type` exposes `gap_gating` and
`gap_gating_threshold_db`. When a gap is evaluated, `dec::Result` carries it in
`cluster_gap_db` with `has_cluster_gap=true`. It is the baseline winner's gap,
even if an expanded ensemble eventually returns another color's correction.
The decoder's statistics report `gap_evaluations` and `gap_expansions`; normal
MPI statistics accumulation includes these counters. `--debug` retains the
actual first/second-stage matching paths and adds the baseline gap, trigger
decision and final winner. There are no matching replays for debug output.

This integration supports ordinary boundaries and one MPI rank per decoder.
Independent-shot MPI sharding still works, but combining gating with
`--concat-parallel` or `--concat-boundary-gap` is explicitly rejected. Existing
parallel and boundary-gap modes work as before when gating is disabled.

The gap is a boundary-sector confidence proxy, rather than a calibrated
probability of final logical failure. A corner boundary-to-boundary edge can
change parity without changing the coalesced virtual syndrome and can cap the
metric. Consequently a useful threshold depends on the noise model and rate:
the prototype's 20 dB gate saves substantial matching work at code-capacity
p=1%, but triggers every shot at p=5%. See the
[paired prototype study](../.agentic_artifacts/concat_gap_gating/REPORT.md)
for failure counts, timings and its separate 12 dB holdout.

Flag-selected event probabilities and logical frames retain their existing
semantics. All configured ensemble graphs are still updated when the flag
pattern changes, even when gating skips their matchings. This integration
therefore does not claim a circuit-level speedup from deferred flag updates.

Validation includes never-expand equivalence to K=1, always-expand equivalence
to full ensembling, selected-threshold decisions, exact vacuum ties, matching
call/coalescing counts, changing SI1000 flag patterns, debug traces, invalid
configuration checks and existing serial/parallel regressions.
