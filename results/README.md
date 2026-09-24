# SI1000 opposite-basis decoder comparison

Author: OpenAI GPT-6-Sol (24 September 2026).

The files named `si1000_d9_r9_p1e-3_opposite_seed20260924.*` contain one
shared-shot comparison of PyMatching and Tesseract. The `.stim` file is the
exact generated circuit. The `.dem` file is its decomposed detector error
model, which is the common input to both decoders and to Stim's `DemSampler`.
The `.jsonl` file has one metadata object, every shot on which PyMatching's
observable prediction disagreed with the sampled truth while Tesseract's
agreed, and one summary object. Detector IDs are zero-based and refer to the
order of `DETECTOR` instructions in the `.stim` file or `D` targets in the DEM.
The `shot_index` is zero-based within this run. The original random seed is
in the metadata object.

Each shot object contains the full fired-detector list, true observable bit,
both decoder predictions, and Tesseract's `low_confidence_flag`. To replay a
record, construct both decoders from the accompanying DEM and give each the
same detector bit vector populated from `detector_ids`. Compare their bit-zero
predictions with `truth`. `opposite_basis_detection_count` counts fired X-check
detectors in the interior body rounds. Detector IDs 0 through 39 are first
round Z checks; each of the next eight rounds has 40 Z-check then 40 X-check
detectors; the final 40 IDs are Z checks after data readout.

The optional files with `uniform_k` in their name use a fixed number of
distinct DEM error instructions selected uniformly without replacement for
each shot. This deliberately increases coverage of rare syndromes. Their
error frequencies are **not** Monte Carlo logical error rates at p=0.001.
The circuit, DEM, decoder settings, and truth convention are otherwise the
same. The metadata's `sampling` and `fixed_error_count` fields distinguish
the method.

Rebuild and rerun with one CPU core:

```sh
cmake --build build_nompi --target compare_si1000 --parallel 1
taskset -c 0 build_nompi/compare_si1000 results/si1000_d9_r9_p1e-3_opposite_seed20260924 1000000 20260924
```

This runner always constructs Z-memory distance-9, nine-round SI1000 with
p=0.001 and `include_opposite_basis_detectors=true`. It uses the decomposed
DEM to sample detector and observable bits. Both decoders process each sampled
row, including rows omitted from the compact JSONL because neither meets the
capture condition. The summary records the exact number of completed rows,
logical errors for each decoder, discordance directions, joint errors, and
Tesseract low-confidence counts. No early stopping is used.

## Completed Monte Carlo run

The run completed exactly **1,000,000** shots using seed 20260924. PyMatching
made 72 logical errors (72 per million); Tesseract made 148 (148 per million).
PyMatching alone erred on 41 shots, Tesseract alone on 117, and both on 31.
The 41 requested PyMatching-only syndromes are all in the JSONL file. These
counts compare both decoders on the same one million sampled detector rows.
Tesseract's `low_confidence_flag` was false on all one million decodes,
including all 41 captured cases.

Under the default Tesseract configuration used here, including the extra
detectors did not produce an overall advantage over PyMatching: Tesseract had
148 errors versus PyMatching's 72. This run does not include a flag-false
control batch, so it cannot measure the effect of adding the detectors on
either decoder's error rate relative to the original circuit.

The 41 captured syndromes are distinct. They contain 21 to 47 fired detectors
each (median 33, mean 34), for 1,394 fired detector IDs in total. Of these,
56 are first-round Z checks, 618 are interior Z checks, 595 are interior X
checks, and 125 are final Z checks. Every captured syndrome contains at least
one interior X detection; the per-shot X count ranges from 6 to 24 (median
14). True observable bits are 21 ones and 20 zeros. The aggregate interior
fired-detector counts by measured round 2 through 9 are respectively 152,
153, 137, 169, 142, 135, 179, and 146. These are descriptive counts over
captured cases, not a model for why the decoders disagree.

The JSONL was parsed and checked with `analyze_si1000.py`: all 41 captured
records have PyMatching prediction different from truth and Tesseract
prediction equal to truth; detector IDs are unique, sorted, and within the
720-detector circuit; summary counts and per-shot detector totals reconcile.
The rebuilt collector produced byte-identical `.stim` and `.dem` files in a
one-shot regeneration using the same parameters, confirming that the saved
records refer to the current circuit generator. The `qudec` CLI flag was
smoke-tested on a distance-3 PyMatching run. The top-level project has CMake
but no Bazel package or module, and `bazel` is absent from this environment;
the vendored dependencies carry their own Bazel files.
Run the validator with:

```sh
python3 results/analyze_si1000.py results/si1000_d9_r9_p1e-3_opposite_seed20260924.jsonl
```

## Paper short-beam comparison

OpenAI GPT-6-Sol added an explicit `paper` mode to the comparison runner to
measure the effect of the Tesseract settings reported in Appendix B of the
Tesseract paper. It uses beam 15, beam climbing, 16 coordinate-projection
detector orders, priority-queue limit 200,000, no-revisit detection, and zero
detection penalty. The order-generation seed is recorded in the JSONL
metadata; these generated orders depend on the C++ standard-library random
implementation. The default mode above retains its original settings.

The paper-mode command is:

```sh
taskset -c 0 build_nompi/compare_si1000 results/si1000_d9_r9_p1e-3_opposite_paper_short_seed20260924 1000000 20260924 paper 20260924
```

The paper treats Tesseract's low-confidence outcomes as logical errors. The
new summary reports both raw prediction errors and
`paper_scored_tesseract_errors` (wrong prediction or low confidence). Only
shots where PyMatching is wrong and Tesseract is correct *with high confidence*
are captured. `uncaptured_raw_only_pymatching` counts any excluded cases where
Tesseract predicts correctly but reports low confidence. The `.stim` and `.dem`
for this mode are identical to the default run's files, SHA-256 hashes
`24925a4d8116f9394bf1086500095cad4d5c40199eca3474aca8b0d82f2ca8c5`
and `f69d93addcefd1fccc2b98f941582c72f5df7485721e92b52f1bbd343235713f`.
