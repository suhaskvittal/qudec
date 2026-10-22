# Upper-triangle affinity inputs and small neural networks

Contributor: OpenAI Codex.

This follow-up to [the first affinity study](affinity-gap-study.md) tests stable
placement, fixed final boundary slots, and a fully connected network on the
upper triangle. It also tests the smaller control of a network on syndrome
bits at fixed detector positions. Numerical detector IDs are not input values.

The [million-shot follow-up](gap-confidence-million-study.md) repeats these
inputs at larger scale and compares gap regression with direct 40 dB binary
classification. The results below use 50,000 training shots.

## Encodings

Both matrix encodings use 48 slots. The two physical boundary classes always
occupy indices 46 and 47, including for even-weight and empty syndromes. The
remaining slots hold active detectors; inactive nodes are masked.

* **Packed triangle:** active detectors occupy sequential slots in detector-ID
  order. Their positions can shift when earlier events are added or removed.
* **Hashed triangle:** a detector's coordinates `(x,y,t)` and memory-component
  membership map to two distinct four-slot buckets. Eleven buckets occupy
  slots 0–43, and slots 44–45 form a shared collision stash. Thus each detector
  has at most ten allowed positions. Deterministic augmenting paths relocate
  occupants only within their allowed candidate sets. Impossible placement
  raises an error; detectors are never merged, overwritten, or dropped.

The hash is computed offline into a candidate lookup table; inference need not
implement BLAKE2b. The present layout is an experimental design, not an FPGA
resource or latency result. On the 50,000 training syndromes, all assignments
succeeded and 33 shots used the stash.

Only the upper triangle **including the diagonal** is retained. There are
`48*49/2 = 1176` values per channel, compared with 2304 in a full square. The
two channels are normalized path strength and a connected/presence mask,
giving **2352 input values**. Active-node and boundary diagonals remain present.
Physical pair affinities are symmetrized before log scaling; boundary entries
contain separate shortest log-odds distances. The normalization is the fixed
`d*ln(1/p)` scale from the first study. Storing the upper triangle preserves
this symmetric representation, rather than averaging collisions or dropping
lower-triangle information without symmetry.

The **raw-syndrome control** has 336 binary inputs, one stable position per
physical detector in this fixed code. It receives no affinity matrix, boundary
distance vector, numerical ID value, or decoded feature. Code geometry is
implicit in the fixed correspondence between a physical detector and an input
position. The **summary control** is the existing 24-feature affinity/boundary
network from the first study.

## Model and evaluation

All four controls use two 64-unit ReLU hidden layers and two linear gap outputs,
with training-only input/target scaling and nonnegative output projection.
Triangle feature standard deviations have a floor of 0.05 to avoid amplifying
rare or previously unused hash slots. The other controls retain their original
training scaling. Both targets use standardized MSE, and validation loss selects
the epoch. Newly trained models use up to 60 epochs with patience ten. Each
family has three training seeds: 101, 202, and 303. The summary checkpoints are
the original three runs, without retraining.

The code/noise configuration remains d=7, seven rounds, SI1000 p=0.001, Z memory,
and opposite-basis detectors included. Training and validation contain 50,000
and 10,000 shots. All twelve checkpoints were frozen before generating a new
20,000-shot holdout with seed 7011. Checkpoint hashes are recorded in
`.agentic_artifacts/affinity_gap/triangle_study/holdout_lock.json` and verified
before inference. This holdout contains 19,596 syndromes absent from training.
Features are created only from the syndrome and fixed DEM.

## Fresh holdout results

Mean MAE across the three training seeds:

| Network input | Input values | Parameters | Complementary MAE (dB) | Cluster MAE (dB) |
|---|---:|---:|---:|---:|
| Packed upper triangle | 2352 | 154882 | 12.22 | 12.13 |
| Hashed upper triangle | 2352 | 154882 | 13.28 | 13.31 |
| Raw syndrome bits | 336 | 25858 | 9.94 | 9.79 |
| Affinity/boundary summaries | 24 | 5890 | 9.53 | 9.42 |

The simple raw-syndrome network is a useful candidate: it approaches the summary
network's accuracy with simpler input preparation. This hash scheme performed
worse than packed placement. This does not establish that every stable layout
or hash design is worse; the candidate buckets, collision policy, data volume,
network capacity, and regularization remain experimental choices. The raw-bit
control supports retaining a stable physical correspondence between detector
and input position.

Upper-triangle storage nearly halves the square-matrix input size. These runs
do not isolate its effect on accuracy against an otherwise identical full-square
MLP. They also do not establish cross-distance generalization, FPGA feasibility,
or an all-syndrome error guarantee. For the validation-selected raw-bit network
(seed 202), the fresh holdout's 99th-percentile absolute complementary-gap error
is 40.03 dB, and the maximum observed is 92.50 dB. Tail behavior remains a major
limitation despite useful average prediction.

## Reproduce

Use the dependencies and environment settings in the first study. Generate
triangle inputs, then train each model/seed combination:

```sh
python3 scripts/affinity_triangle_study.py prepare \
    --reference .agentic_artifacts/affinity_gap/prepared_d7 \
    --output .agentic_artifacts/affinity_gap/triangle_study/prepared \
    --size 48 --hash-seed 17
python3 scripts/affinity_gap_study.py train \
    --prepared .agentic_artifacts/affinity_gap/triangle_study/prepared \
    --model mlp_hashed_triangle --seed 101 --threads 4 \
    --epochs 60 --patience 10 \
    --output .agentic_artifacts/affinity_gap/triangle_study/mlp_hashed_triangle_s101
```

The other new model names are `mlp_packed_triangle` and `mlp_syndrome`.
To build new test features without rewriting training data:

```sh
python3 scripts/affinity_triangle_study.py prepare \
    --reference .agentic_artifacts/affinity_gap/prepared_d7 \
    --reuse-training .agentic_artifacts/affinity_gap/triangle_study/prepared \
    --test .agentic_artifacts/affinity_gap/triangle_study/holdout_d7 \
    --output .agentic_artifacts/affinity_gap/triangle_study/final_prepared
python3 scripts/affinity_triangle_study.py evaluate \
    --prepared .agentic_artifacts/affinity_gap/triangle_study/final_prepared \
    --runs .agentic_artifacts/affinity_gap/triangle_study/mlp_hashed_triangle_s101 \
    --output .agentic_artifacts/affinity_gap/triangle_study/final_results
```

Evaluation verifies checkpoint layout, topology, normalization, and training /
validation metadata. Old summary checkpoints skip only triangle-layout checks.
The completed `final_results` contains all twelve runs. Reports, CSV metrics,
and plots are under `.agentic_artifacts/affinity_gap/triangle_study/report/`.

The encoding audit compared triangle entries against an independently permuted
dense matrix, checked the full prepared arrays for finite values in [0,1], and
verified that hash assignments stay in each detector's candidate set. All
dataset preparation and final inference runs completed successfully without
truncation or collision merging.
