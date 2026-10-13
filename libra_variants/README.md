# Libra matching experiments

This entry point preserves the existing Libra implementations and adds
`dec::LibraVariant`, with two independent configuration options:

- `large_gap_fallback`: run the unconstrained matcher when the initial
  complementary gap is at least 20 dB. Disabling it returns the cheaper
  constrained representative on these shots.
- `correlated_matching`: enable PyMatching's correlated two-pass matching for
  all initial, ensemble, and fallback candidates. Disabling it uses the
  ordinary single-pass edge-output API and graphs without correlation rules.

Both options default to `true`, reproducing the preserved PyMatching-backed
Libra implementation. Hyperedge assignment, original-prior scoring, ensemble
perturbations, cycle synthesis, and degeneracy estimation remain identical.
The gap and therefore the ensemble trigger can change when matching changes.

Build with:

```sh
cmake -S libra_variants -B .agentic_artifacts/libra_variants/build \
  -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=OFF
cmake --build .agentic_artifacts/libra_variants/build \
  --target libra_variant_evaluate libra_variant_test -j 16
ctest --test-dir .agentic_artifacts/libra_variants/build -R '^libra_variant_' \
  --output-on-failure
```

The paired Monte Carlo evaluator accepts:

```text
libra_variant_evaluate fallback|correlation DIST WORKERS MAX_SHOTS SAMPLE_SEED PREFIX ERROR_TARGET KEEP_FALLBACK [ENSEMBLE_SEED]
```

`fallback` compares correlated Libra with and without the extra unconstrained
call. `correlation` compares correlated and uncorrelated matching, both with
the specified fallback policy (`KEEP_FALLBACK` is 0 or 1). Ordinary PyMatching
also decodes every shot. The SI1000 physical error rate is 0.001, rounds equal
distance, Z memory, and opposing basis detectors are enabled.

`ENSEMBLE_SEED` defaults to `SAMPLE_SEED`. Specify it separately to validate an
unchanged decoder on independent samples, retaining the same perturbed priors.

Each worker has a fixed sampling stream. All decoders receive exactly the
same syndrome and truth on each shot; truth is used only to score predictions.
The run stops when both Libra variants reach the error target or the shot cap.
Parallel scheduling can change the endpoint and worker stream lengths, so
the output records worker and local shot indices for each failure.

The evaluator refuses existing output filenames and creates:

- `PREFIX.txt`: aggregate results and mean decode times.
- `PREFIX_hw.csv`: exact fired-detector Hamming weight bins, including every
  bin's shot count, failures, disagreements, ensemble triggers, and timing.
- `PREFIX_fault_hw.csv`: analogous bins by sampled DEM channel count.
- `PREFIX_joint_hw.csv`: joint detector and DEM channel count bins.
- `PREFIX_failures.csv`: all shots where any decoder failed, including
  syndrome/channel IDs, predictions, gaps, and ensemble triggers for replay.
- `PREFIX_thresholds.csv`: cumulative results for dispatching the alternative
  when detector Hamming weight is at most the threshold, otherwise using
  correlated Libra. Timing excludes the small threshold dispatch overhead.

DEM channel count is a simulation diagnostic, not an observed quantity and
not an exact count of physical gate faults. Fired-detector count is available
to an adaptive decoder. Thresholds selected from these data require separate
samples before treating them as validated operating points.

See `Libra-Ablations.txt` for the comparison and threshold analysis results.
