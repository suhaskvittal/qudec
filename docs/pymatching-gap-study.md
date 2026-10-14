# PyMatching gap studies

Contributor: OpenAI Codex.

## Complementary gap statistics

`scripts/run_gap_study.py` collects complementary-gap statistics from standalone
uncorrelated PyMatching. It does not run the two-level evaluator or fit a routing
threshold.

```sh
python3 scripts/run_gap_study.py --distance 9 --shots 1000000 --cores 64
```

Specify `--qudec PATH` or a configured `--build-dir DIR` when needed. Options
include `--rounds`, `--p`, `--seed`, `--max-errors`, `--mpi-launcher`, and
`--output-root`. Use `--cores 1` for a serial run. `--dry-run` prints the
command without creating an output folder. Each run creates a unique
`out/complementary_gap_.../` folder with the full sampling log, commands,
`summary.json`, and `gap_bins.csv`. The summary reports observed logical error
rate and signed/unsigned gap counts, mean, standard deviation, minimum and
maximum when available. Empty histograms use null statistics.

The log includes `SIGNED_GAP`, `UNSIGNED_GAP`, and `UNSIGNED_GAP_ERRORS`; the
latter reports error counts in each gap bin and supports conditional-error
analyses. The physical logical error rate is the observed error count divided
by observed shots. This script always uses Monte Carlo sampling.

## Two-level Libra evaluator

Build `qudec_eval_gap` and run it directly when evaluating a first decoder plus
Libra:

```sh
mpirun -np 64 PATH/TO/qudec_eval_gap 9 --mc-max-samples 1000000 \
  --gap-metric complementary --first-decoder pymatching --gap-threshold-db 20
```

The default metric remains exact PyMatching complementary gap. Its sign encodes
whether PyMatching's correction agrees with truth, so routing uses its absolute
value and equality stays with the first decoder. Ground truth contributes only
to error measurements.

Cluster-gap routing uses decoder soft output independent of truth:

```sh
PATH/TO/qudec_eval_gap 9 --mc-max-samples 1000 \
  --gap-metric cluster --first-decoder pymatching --gap-threshold-db 20
PATH/TO/qudec_eval_gap 9 --mc-max-samples 1000 \
  --gap-metric cluster --first-decoder uf --gap-threshold-db 20
```

`--first-decoder uf` is currently supported with `--gap-metric cluster` only.
`--uf-weight-scale` selects UF's log-weight quantization, default 1024. Both
decoders use SI1000 Z memory, p=0.001, and rounds equal to distance by default;
the evaluator includes opposite-basis detectors. The evaluator continues to
run Libra on low-gap shots and reports first-stage and two-level outcomes.

The complementary metric is retained as a separate option and is not replaced
by cluster gap. Cluster gap uses `Result.cluster_gap_db`, while complementary
routing uses `abs(Result.matching_data.gap)`. Correlated PyMatching is not
enabled in this evaluator.

The evaluator stores a unique run directory under `out/` by default, or the new
directory given by `--output-dir`. It writes `summary.txt`, per-level
`error_count_levels.csv`, and per-rank low-gap syndrome and result files. The
low-gap records include the selected route metric, first-decoder latency and
frame, Libra latency and frame, and truth. The reported two-level LER counts
first-decoder failures on high-gap shots plus Libra failures on low-gap shots.
Its requested combined mean latency uses first-decoder time for high-gap shots
and Libra time for low-gap shots; the sequential mean also includes the first
decoder's gap calculation on those low-gap shots. All latency values cover the
full multi-round history.

Rare-event evaluation computes probability-weighted error and latency estimates
by active-error count. Its raw per-shot histograms remain conditional samples,
not the physical gap distribution. The cluster-gap statistics script rejects
rare-event mode and therefore reports only Monte Carlo empirical distributions;
do not interpret unweighted rare-event counts as physical probabilities. For
two-level rare-event runs, physical errors and means use
`sum_k P(K=k) * conditional_mean_k` without renormalizing the retained mass.
The omitted count-probability bound limits the missing logical-error
contribution but does not bound omitted latency.
