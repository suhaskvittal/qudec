# Cluster soft output

Contributor: OpenAI Codex.

`qudec pymatching DISTANCE --cluster-gap` and `qudec uf DISTANCE --cluster-gap`
compute a full, nonnegative cluster gap for every decoded syndrome. Add
`--include-opposite-basis-detectors` when using that surface-code detector model.
The default output contains a compact summary; `--cluster-gap-histogram` also
prints the nonempty histogram bins.

C++ configuration: set `enable_cluster_gap=true` in `PyMatching::config_type` or
`UnionFind::config_type`. Each result contains `has_cluster_gap=true` and
`cluster_gap_db`. With the feature disabled, `has_cluster_gap=false` and no
cluster searches run. `print_cluster_gap_histogram` controls printed bins.

The implementation follows the weighted cluster geometry described in
[Meister, Pattison, and Preskill, *Efficient soft-output decoders for the surface
code*](https://arxiv.org/html/2405.07433v2). After treating the covered portions of
edges as free, the gap is the shortest remaining weighted path between the two
inequivalent physical boundaries. Partial coverage is preserved. Full values
are returned; they are not capped at a routing threshold.

## Decoder-specific clusters

PyMatching uses its ordinary uncorrelated matching solve. A callback before
Blossom extraction retains each fired detector's radius, including contributions
from every containing Blossom. Multisource weighted propagation reconstructs
the union of metric balls. This reconstruction permits propagation through each
collapsed physical boundary and therefore can conservatively over-cover the
cluster geometry. Its output can be smaller than the ideal paper quantity.
Opposing-basis detectors are excluded from the memory confidence graph, while
ordinary matching still processes all supplied detectors and retains endpoint
reporting when requested.

Union-find uses its own grown erasure intervals: the uncovered edge length is
`edge.length - growth_[edge]`, including partially grown frontier edges. Peeling
and endpoint reporting are unchanged. The current weighted UF growth schedule
is the existing decoder's simultaneous growth, not the paper's illustrative
smallest-cluster schedule. UF's cluster output is not claimed to bound an MWPM
complementary gap.

## Units and restrictions

The output uses decibels: natural-log weight is multiplied by `10/ln(10)`.
PyMatching uses its matching graph's normalization constant. UF uses its doubled
integer edge lengths and configured `weight_scale`. `MatchingData::gap` retains
its existing exact complementary-gap meaning and is not repurposed.

Initial support is one planar memory logical observable whose logical labels
occur only on single-detector boundary components. Non-graphlike components,
undetectable logical-only channels, missing inequivalent boundaries, and negative
log-odds weights are rejected. `--cluster-gap` cannot be combined with `--gap` or
`--correlated`, and currently applies only to standalone `pymatching` and `uf`.

A cluster value is soft information, not a calibrated conditional failure
probability. The MWPM cluster construction has a lower-bound relationship to
complementary weights under the paper's graph/dual assumptions. Tests against
exact gaps support the implementation, but do not prove a circuit-level SI1000
failure certificate or account for correlations discarded by independent graph
marginalization.

## Timing

Existing `decode_latency_ns` and printed decoder latency include the full
ordinary solve, cluster extraction/search, and requested endpoint reporting.
They describe the entire syndrome history, not one round. MPI summaries reduce
all ranks. A decoder instance owns mutable scratch buffers and must not be used
concurrently by multiple threads.

Earlier artifact measurements of approximately 11–12 microseconds were capped
threshold gates, and measurements of approximately 41–43 microseconds were
uncapped memory-component diagnostics. Neither is a production full-detector
latency measurement. Use the integrated decoder's timing for comparisons.

## Monte Carlo distribution study

Contributor: OpenAI Codex.

Use `scripts/run_cluster_gap_study.py` to collect a distribution for either
standalone decoder:

```sh
python3 scripts/run_cluster_gap_study.py --decoder pymatching --distance 9 --shots 1000000 --cores 64
python3 scripts/run_cluster_gap_study.py --decoder uf --distance 9 --shots 1000000 --cores 64
```

The script accepts `--build-dir` or `--qudec`, plus `--rounds`, `--p`, `--seed`,
`--max-errors`, `--uf-weight-scale`, `--mpi-launcher`, and `--output-root`.
Use `--cores 1` for serial operation; `--dry-run` prints the command. Each run
gets a unique `out/cluster_gap_<decoder>_.../` folder containing the sampling
log, command, `summary.json`, and `cluster_gap_bins.csv`. The summary includes
observed shot count, mean, standard deviation, minimum and maximum when
available. New runs also record a `CLUSTER_GAP_DB_ERRORS` histogram: logical
errors in each of the same 1 dB gap bins.

Fit the conditional logical error probability above a gap threshold with:

```sh
python3 scripts/fit_cluster_gap_threshold.py out/cluster_gap_pymatching_d9_...
```

This writes `cluster_gap_thresholds.csv` and `cluster_gap_threshold_fit.json`
into the run folder. For each bin lower edge `G`, it computes
`P(logical error | cluster gap >= G)` from the cumulative error and shot counts.
It fits `log10(P) = intercept + slope_per_db * G` by ordinary least squares
across thresholds with at least five errors and 100 shots; change these limits
with `--min-errors` and `--min-shots`. Zero-error thresholds remain in the CSV
but cannot enter a logarithmic fit. Cumulative threshold estimates share shots,
so the reported `r_squared` describes fit shape, not independent observations.
Existing run folders without per-gap logical error counts cannot support this
fit and must be collected again with the updated binary and runner.

These are empirical Monte Carlo probabilities under the selected SI1000 error
rate. The runner deliberately has no rare-event option: per-level conditional
histogram counts from rare-event sampling require weighting by the corresponding
`P(K)` masses before they represent the physical distribution.
