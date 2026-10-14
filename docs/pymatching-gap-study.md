# PyMatching complementary-gap study

Contributor: OpenAI Codex.

Build the `qudec_eval_gap` CMake target. For example:

```sh
mpirun -np 64 PATH/TO/qudec_eval_gap 9 --mc-max-samples 1000000 --gap-threshold-db 20
```

The defaults are SI1000 Z memory, p=0.001, rounds equal to distance, and
opposing-basis detectors enabled. The first decoder is uncorrelated PyMatching
with complementary-gap estimation enabled. `--gap-threshold-db` routes a shot
to Libra exactly when the **absolute** PyMatching gap is below the threshold;
equality stays with PyMatching. Ground truth is used only for error statistics.

## Automatically choose and evaluate a threshold

```sh
python3 scripts/run_gap_study.py --distance 9 --shots 1000000 --target-error-rate 1e-8
```

This replaces the former separate collection and fitting scripts. The new
workflow collects Monte Carlo calibration data from `qudec pymatching --gap`,
fits the conditional error rate versus the absolute gap, computes the threshold,
and runs `qudec_eval_gap` on separate rank seed streams. It defaults to 64 ranks,
p=0.001, rounds equal to distance, up to 100 million calibration shots or 10,000
calibration errors, and a target conditional error rate of 6.9e-11 (the former
script's default). Every run has its own `out/gap_study_.../` directory.

The rule fits `log10(P(error | gap)) = slope * gap + intercept`, weighted by
each bin's error count, then solves for the gap at `--target-error-rate`.
`--min-errors` defaults to 5 observed errors per fitted bin; at least two bins
and a decreasing fitted slope are required. Actual histogram bounds are used,
including the current two-dB PyMatching bins. Missing zero bins are supported,
and lumped overflow bins are excluded from fitting. Insufficient calibration
data produce an explicit failure rather than an invented threshold.

The fitted condition concerns a per-gap error probability; it does not
guarantee the overall decoder error rate. A threshold outside the observed fit range is marked as
extrapolated; extremely small target rates will typically require extrapolation.
The output includes `threshold.json`, `fit_points.csv`, calibration/evaluation
logs, a manifest, and `study_summary.txt`. The evaluation subdirectory contains
the executable's distributions, error counts, and retained syndromes.

Use `--calibration-shots`, `--calibration-max-errors`, `--cores`, `--p`, and
`--rounds` to change the experiment. `--calibration-log FILE` reuses a completed
Monte Carlo PyMatching gap log; the user must ensure its distance, rounds, error
model and rate match this experiment. `--method rare_event` changes the second
evaluation stage's sampler; calibration always uses Monte Carlo, because raw
rare-event histogram counts do not represent physical conditional probabilities.
`--dry-run` prints commands without creating a run or sampling.

The script discovers built executables, or accepts `--qudec PATH` and
`--evaluator PATH`. With `--build-dir DIR`, it first builds both targets in that
already-configured CMake directory. NumPy and an Open MPI launcher are required
for the default configuration; use `--cores 1` for serial runs.

To provide a threshold directly and skip calibration and fitting:

```sh
python3 scripts/run_gap_study.py --distance 9 --shots 1000000 --gap-threshold-db 20
```

Direct-threshold runs need only `qudec_eval_gap`. The final
`two_level_logical_error_rate` is printed to the console and included in the
saved summaries. It counts high-gap PyMatching failures plus low-gap Libra
failures, divided by total shots for Monte Carlo or weighted by count-level
probability for rare-event evaluation.

## Two stages

1. Run PyMatching with gap estimation on every sampled syndrome. Split its
   per-shot decoding times into high-gap and low-gap distributions. Retain
   copies of all low-gap syndromes and their true logical frames.
2. After first-stage sampling finishes, run the retained original syndromes
   through Libra, recording its latency and frame. There is no PDZ processing.
   Libra retains its own confidence/ensemble policy; configure it independently
   with `--libra-gap-threshold-db` and `--libra-ensemble-size`.

`Result::decode_latency_ns` reports the existing PyMatching and Libra internal
timers, without a second external timer. PyMatching timings include both
complementary-class matching calls and its gap calculation. Constructors,
Libra ensemble preparation, retained-syndrome copies, sampling, and file writes
are outside the reported decoding times.

The requested combined mean is:

```text
(sum(PyMatching latency on high-gap shots) + sum(Libra latency on low-gap shots))
/ total sampled shots
```

This treats Libra latency as the replacement cost for low-gap shots. The report
also prints the sequential mean, which adds PyMatching gap-estimation latency
on those shots, since the executable itself actually runs both decoders there.
Both metrics cover a complete multi-round syndrome history, not one round.

## Sampling and output

Monte Carlo uses the existing sampler and completes one million shots by
default. `--mc-stop-at-errors` optionally stops based on first-stage PyMatching
failures. Use `--method rare_event` for the corrected heterogeneous conditional
sampler, with the same samples-per-level, errors-per-level, count limit, and
tail-tolerance controls. Early stopping uses the running **PyMatching** LER;
the report states the omitted probability bound, which also bounds the missing
two-level error contribution but does not bound omitted latency. Set
`--rare-relative-tail-tolerance 0` to disable this early stop.

Rare-event raw histograms describe the conditional shots actually drawn. Its
physical mean and LER are computed as `sum_k P(K=k) * conditional_mean_k` using
global per-level shot counts, without renormalizing retained probability. This
is distinct from the unweighted sample mean. MPI reductions combine counts,
timing totals, and histograms across workers before computing means.

Each run creates a new directory under `out/`, or the new directory supplied
by `--output-dir`. Existing run directories are rejected. It contains:

- `summary.txt`: configuration, counts, accuracy, mean latencies, and all three
  latency histograms (nonempty bins, plus populated overflow bins).
- `error_count_levels.csv`: per-level probability, counts, error counts, and
  integer nanosecond timing totals. Monte Carlo has one synthetic level of mass 1.
- `low_gap_syndromes_rank_N.b8`: retained detector syndromes in Stim b8 order;
  each record has `ceil(num_detectors / 8)` bytes with zero padding bits.
- `low_gap_records_rank_N.csv`: matching row indices, count level, PyMatching
  gap, both latencies, initial and final frames, and the true frame.

All retained low-gap syndromes are held in memory until stage two. If the
threshold forwards most shots, memory use grows with the shot budget and code
size. The latency histogram defaults to one-microsecond bins up to 20,000 us;
`--latency-bin-us` and `--latency-max-us` adjust that range.
Histograms are saved in `summary.txt`; console bins are hidden by default.
Use `--print-histograms` to display them. For example,
`--latency-bin-us 100 --latency-max-us 2000` gives 20 bins and an overflow
bucket for latencies at least 2,000 us. Overflow samples remain included in
all reported mean latencies.
