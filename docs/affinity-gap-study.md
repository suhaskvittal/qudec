# Predicting surface-code gaps before decoding

Contributor: OpenAI Codex.

This study asks whether syndrome and fixed detector-error-model (DEM) features
predict an **unsigned complementary gap** or a **PyMatching cluster gap**. Labels
are computed offline; feature extraction and inference do not call a decoder.
The initial experiment uses rotated surface-code Z memory, distance 7, seven
rounds, SI1000 noise at p=0.001, and both detector bases. Cross-distance
generalization is a later experiment.

The [million-shot follow-up](gap-confidence-million-study.md) tests data-size
scaling and the fixed 40 dB screening decision. The results below are the
initial 50,000-shot pilot, rather than the largest completed experiment.

Gap is a measure of logical ambiguity. It is not itself decoder running time.
In the complementary target, the two logical classes have minimum matching
weights W0 and W1, and the label is `abs(W0-W1) * 10/ln(10)`. The signed gap in
the repository encodes agreement with ground truth and must not be a prediction
target or input. The cluster target uses this repository's uncorrelated
PyMatching cluster reconstruction; see [cluster-gap.md](cluster-gap.md) for its
boundary propagation and correlation limitations. The two targets are retained
separately. Context for cluster soft outputs is
[Meister, Pattison, and Preskill](https://arxiv.org/abs/2405.07433).

## What affinity represents

The repository has two different representations:

* `pp::measure_affinity` projects the undecomposed DEM's hyperedges into pair
  adjacencies, adds their probabilities, and finds shortest negative-log paths.
  Its output is the product of projected probabilities along a best path.
  It is a connectivity heuristic, not a sum over error chains or an exact
  pair probability. A shared boundary is absorbing, preventing a shortcut
  through that vertex.
* `pp::AffinityTable` stores unweighted hop distance in bytes. Values 0 through
  30 are exact; 31 represents saturation or disconnection. It discards edge
  probabilities and also has a single shared boundary.

Close pairs suggest that defects can be explained locally. Several similarly
close partners suggest competing explanations. Neither observation alone
identifies whether competing explanations cross a logical boundary. A pair
matrix with a shared boundary can lose the distinction between the memory
component and the other basis, as well as the distance to the other physical
boundary. This is a property of the representation; increasing model capacity
cannot recover information absent from it.

## Pipeline and controls

`export_affinity_dataset` samples the raw circuit DEM, preserving correlated
fault supports, and runs three offline uncorrelated PyMatching instances on the
decomposed DEM: ordinary decoding, complementary gap, and cluster gap. It writes
the syndrome, labels, fixed lookup tables, circuit/DEMs, detector coordinates,
and unmerged graph components. Output directories must not exist. Metadata is
written last so its presence marks completed output.

The Python script requires separate training, validation, and test seeds and
checks circuit, DEM, coordinate, graph-component, and lookup-table identity
between splits. The pilot uses 50,000 / 10,000 / 20,000 shots with seeds
7001 / 7002 / 7003. Repeated syndromes are valid IID samples; test metrics are
also reported on syndromes absent from training. Features are created solely
from syndrome and DEM. The truth, decoded observable, timings, and both gap
labels are excluded from input features.

The initial CNN inputs pack active detectors in ascending detector-ID order
after persistent boundary anchors and pad to 48 by 48. This is an **active-event
affinity matrix**, with axes representing event rank, not spatial lattice
coordinates. Translation invariance is therefore only a model hypothesis.
There is no truncation: preparation fails if a syndrome exceeds capacity.
Changing event order changes the convolution's neighborhoods; the current study
does not establish permutation or code-distance invariance.

Input variants are:

| Variant | Pair features | Boundary information |
|---|---|---|
| `cnn_hop` | Hop distance | Shared absorbing boundary |
| `cnn_weighted` | Negative log raw-DEM affinity | Shared absorbing boundary |
| `cnn_boundary` | Negative log raw-DEM affinity | Two labeled physical boundaries |
| `cnn_context` | Same as `cnn_boundary` | Both endpoints' boundary costs broadcast to every pair pixel |

Each matrix has a path-strength channel, connected/presence mask, and node-role
diagonal. Shared-boundary values are symmetrized from detector-to-boundary
values. Weighted distances are divided by the fixed `d*ln(1/p)` scale; hops by
`d`. Path strengths are clipped to `[0,1]`; these are engineering
transformations, not calibrated probabilities. No per-syndrome normalization
removes absolute distance. The context representation adds four channels.

Two-boundary distances are precomputed from the decomposed graph's natural-log
odds. Independent repeated edges are merged using XOR probabilities. Directed
shortest paths start at a boundary but cannot traverse a boundary to reenter a
component. Interior pair features still come from raw-DEM affinity, so boundary
augmentation also introduces information from the matching graph. It is not an
isolated test of boundary labels on otherwise identical graph metrics.

Baselines are a training mean, a nonlinear syndrome-weight predictor, ridge
regression and extra trees on affinity summaries, extra trees on summaries plus
basis counts (`trees_basis`) or full boundary features (`trees_boundary`), and
extra trees on the complete syndrome bit vector. The
summary baseline includes effective-partner count `(sum a)^2/sum(a^2)`, nearest
and second-nearest log distances, pair distances, and shared-boundary distances.
The basis-count control adds memory/opposite-component weights without the
distinct boundary distances, separating that effect from boundary geometry.
The boundary baseline includes memory-component weight and distinct
boundary distances. Trees are research baselines, not an FPGA implementation
proposal. The full-syndrome tree baseline uses square-root feature subsampling;
the low-dimensional tree baselines consider all features. `mlp_boundary` is a
small neural control with two 64-unit ReLU layers on the standardized full
summary vector. It tests whether a convolution is needed for these features.

The small CNN uses two convolutions, ReLU, average pooling, and two dense layers.
Both targets are standardized using training statistics; training minimizes
their mean squared errors jointly. Validation loss selects the epoch. Prediction
is projected to nonnegative values for every model. Test labels do not select
epochs or model parameters. The plotted CNN is selected by mean validation MAE
across both targets.

Reported metrics include MAE, RMSE, R-squared, Spearman correlation, novel
syndrome metrics, and bootstrap MAE intervals conditional on a fixed trained
model. Low-gap labels use the validation tenth percentile. The routing
prediction threshold targets 95% validation recall, then is frozen for
test recall, precision, and routed fraction. This quantile is a pilot operating
point, not a logical-failure guarantee. Retrospective test routing curves are
descriptive. Training-seed repeats are required before treating small model
differences as established.

## Pilot findings and fresh holdout

The pilot's independently seeded 20,000-shot test set was inspected during
development of the context CNN and MLP controls. Consequently its results are
exploratory. After all models were frozen, an additional 20,000-shot holdout
with seed 7010 was generated. `.agentic_artifacts/affinity_gap/holdout_lock.json`
records the model hashes and the validation-selected CNN and MLP before those
shots were generated. The holdout evaluator applies these existing models
without training or choosing new parameters.

The fresh holdout confirms the exploratory comparison:

| Model/input | Complementary MAE (dB) | Cluster MAE (dB) |
|---|---:|---:|
| Syndrome weight | 14.53 | 14.71 |
| Affinity summary trees | 13.29 | 13.40 |
| Add basis counts | 11.67 | 11.70 |
| Add distinct boundary distances | 9.44 | 9.43 |
| CNN: weighted shared boundary | 14.61 | 14.77 |
| CNN: separate boundary anchors | 12.72 | 12.72 |
| CNN: broadcast boundary context | 12.08 | 11.99 |
| Summary MLP, validation-selected seed 101 | 9.59 | 9.49 |

The three MLP training seeds averaged 9.64 / 9.53 dB MAE for complementary /
cluster gaps. Its roughly 5,900 parameters are fewer than the CNN's roughly
39,000. These are initial architectures, not a claim that CNNs cannot work
for this problem. There were 19,613 novel syndromes relative to training in
the 20,000-shot final holdout.

Average error conceals substantial tail error: the seed-101 MLP's fresh-holdout
99th-percentile absolute complementary-gap error was 39.43 dB and its largest
observed absolute error was 91.84 dB. Neither these measurements nor a
statistical calibration would establish a guarantee for every syndrome.

Affinity concentration is partly a proxy for syndrome size. On the pilot test
set, mean effective-partner count had Spearman correlation -0.402 with
complementary gap across all shots but only -0.098 among the 1,685 shots of
weight eight (the most frequent nonzero weight). Nearest-versus-second-nearest
margin had corresponding correlations 0.360 and 0.106. This supports comparing
against syndrome weight and retaining boundary/topology information.

There is also a concrete obstruction to exact prediction with the current
shared-boundary active-event representation. At d=7, single fired memory-basis
detectors 324 and 12 produce **identical float32 weighted CNN tensors**, while
their complementary/cluster labels are 62.0063 and 76.9575 dB. Distinct-boundary
tensors distinguish them. Any deterministic predictor of the identical input
must incur at least half their 14.9511 dB separation on one of them: a uniform
absolute-error guarantee below 7.4756 dB is impossible for that encoding.
The analogous hop-only collision between memory detectors 324 and 9 has labels
62.0063 and 127.0228 dB. This concerns the tested encoding; it does not imply
that a fixed-detector-ID matrix retaining position has the same limitation.
The example and comparison figure are saved in
`.agentic_artifacts/affinity_gap/memory_only_collision.json` and
`.agentic_artifacts/affinity_gap/boundary_collision.png`.

## Reproduce the datasets and training

All generated output, dependencies, caches, logs, and build files should stay
inside `.agentic_artifacts/`. The following uses the initial artifact paths:

```sh
mkdir -p .agentic_artifacts/affinity_gap/tmp
export TMPDIR="$PWD/.agentic_artifacts/affinity_gap/tmp"
export PIP_CACHE_DIR="$PWD/.agentic_artifacts/affinity_gap/pip_cache"
export MPLCONFIGDIR="$PWD/.agentic_artifacts/affinity_gap/mpl"
export PYTHONDONTWRITEBYTECODE=1
python3 -m pip install --target .agentic_artifacts/affinity_gap/python \
    --index-url https://download.pytorch.org/whl/cpu torch
python3 -m pip install --target .agentic_artifacts/affinity_gap/python \
    numpy scipy scikit-learn matplotlib
export PYTHONPATH="$PWD/.agentic_artifacts/affinity_gap/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

cmake -S . -B .agentic_artifacts/affinity_gap/build \
    -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=OFF -DBUILD_TESTING=OFF
cmake --build .agentic_artifacts/affinity_gap/build \
    --target export_affinity_dataset -j 16

.agentic_artifacts/affinity_gap/build/export_affinity_dataset \
    --distance 7 --shots 50000 --seed 7001 \
    --output-dir .agentic_artifacts/affinity_gap/train_d7
.agentic_artifacts/affinity_gap/build/export_affinity_dataset \
    --distance 7 --shots 10000 --seed 7002 \
    --output-dir .agentic_artifacts/affinity_gap/validation_d7
.agentic_artifacts/affinity_gap/build/export_affinity_dataset \
    --distance 7 --shots 20000 --seed 7003 \
    --output-dir .agentic_artifacts/affinity_gap/test_d7

python3 scripts/affinity_gap_study.py prepare \
    --train .agentic_artifacts/affinity_gap/train_d7 \
    --validation .agentic_artifacts/affinity_gap/validation_d7 \
    --test .agentic_artifacts/affinity_gap/test_d7 \
    --output .agentic_artifacts/affinity_gap/prepared_d7
python3 scripts/affinity_gap_study.py train \
    --prepared .agentic_artifacts/affinity_gap/prepared_d7 \
    --model cnn_context --seed 101 --epochs 40 --threads 1 \
    --output .agentic_artifacts/affinity_gap/cnn_context_s101
python3 scripts/plot_affinity_gap_study.py \
    --prepared .agentic_artifacts/affinity_gap/prepared_d7 \
    --runs .agentic_artifacts/affinity_gap/cnn_context_s101 \
    --output .agentic_artifacts/affinity_gap/report
```

Repeat `train` with the other model names and pass all completed run directories
to `--runs`. A single thread avoids reliance on OpenMP scaling. Controlled
PyTorch CPU threads may be requested using `--threads` after checking the local
installation. The current pipeline is serial per dataset exporter; independent
seeds/model runs can execute as separate processes. MPI is not required.

The pilot used Python 3.12, CPU PyTorch 2.14.1, NumPy 2.5.4, SciPy 1.18.1,
scikit-learn 1.9.1, and Matplotlib 3.11.2. Four PyTorch threads worked in this
installation. Exact installed versions are recorded in
`.agentic_artifacts/affinity_gap/environment.json`.

To apply frozen models to a fresh exported dataset, use the holdout evaluator.
It verifies the circuit/graph hashes and the new seed, builds identical
features, loads existing checkpoints, and reuses validation operating points.

```sh
python3 scripts/evaluate_affinity_gap_holdout.py \
    --dataset .agentic_artifacts/affinity_gap/holdout_d7 \
    --reference .agentic_artifacts/affinity_gap/prepared_d7 \
    --runs .agentic_artifacts/affinity_gap/cnn_context_s101 \
           .agentic_artifacts/affinity_gap/mlp_boundary_s101 \
    --output .agentic_artifacts/affinity_gap/final_holdout_v2
python3 scripts/plot_affinity_gap_study.py \
    --prepared .agentic_artifacts/affinity_gap/final_holdout_v2/prepared \
    --runs .agentic_artifacts/affinity_gap/final_holdout_v2/cnn_context_s101 \
           .agentic_artifacts/affinity_gap/final_holdout_v2/mlp_boundary_s101 \
    --output .agentic_artifacts/affinity_gap/final_report
```

The completed study's `final_report` includes all frozen baselines and three
MLP seeds as well as the CNNs. Its `metrics.csv` and the individual run
`metrics.json` files report the fresh holdout results; the `report` directory
contains the earlier pilot plots. Generated datasets and reports are ignored
by Git, while the exporter, study, evaluator, plotting script, and this protocol
are lasting project changes.

The exporter binary files are row major, with `labels.f64` and `affinity.f64`
encoded as little-endian float64. `syndromes.u8` is one byte per detector per
shot. `metadata.json` describes shapes and label columns. Prepared `.npy`
matrices are memory mapped during training. Checkpoints retain target scaling
and protocol metadata alongside model parameters.

## Hardware and further experiments

The first question is whether the representation predicts the gap well enough.
Precomputed pair and boundary tables, event compaction, fixed masks, small
convolutions, and ReLU/pooling/dense operations provide a path to an FPGA design.
The experiment does not validate fixed-point accuracy, resources, throughput,
or timing. Quantization, finite table precision, event overflow behavior, and
operator export should follow a working predictive design.

Further scientific controls include matched-metric boundary ablations, explicit
coordinates/basis channels, event-order permutations, a graph/set model with
boundary anchors, more training data and seeds, and physical-noise and
cross-distance shifts. Low-gap accuracy and actual failure enrichment need
larger held-out samples than the present pilot.
