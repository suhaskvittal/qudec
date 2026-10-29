# Million-shot study of pre-decoding gap confidence

The `export_affinity_dataset` commands below document a historical study. The
standalone exporter has been removed from the current build.

Contributor: OpenAI Codex.

This follows the [initial affinity study](affinity-gap-study.md) and the
[upper-triangle study](affinity-triangle-study.md). The task is to predict
whether an unsigned complementary gap or a PyMatching cluster gap is at most
40 dB, using only the syndrome and precomputed code/noise information.
Equality belongs to the low-confidence class. Each target has its own output
and decision threshold. Gap confidence is a proxy for logical ambiguity;
these experiments do not predict decoder running time directly.

The [distance-9 comparison](gap-confidence-distance-study.md) tests the same
training-shot budget at a higher distance and records the stopping decision.

## Data and experiment

The configuration is rotated surface-code Z memory, distance 7, seven rounds,
SI1000 noise at p=0.001, including opposite-basis detectors. The fixed circuit
has 336 detectors. Training has 1,000,000 independent shots (seed 8101), and
validation has 100,000 shots (seed 8102). A separate 200,000-shot test uses
seed 8103 after checkpoints and thresholds are frozen. All three use the
same circuit and detector error model (DEM).

Sampling uses Stim's sampler on the **undecomposed circuit DEM**, preserving
the support of correlated faults. This avoids simulating the circuit for each
shot. Direct DEM sampling is still Monte Carlo sampling from the physical
distribution; it is not enumeration or deliberate enrichment of rare faults.
The ordinary, complementary, and cluster decoders run offline on the
decomposed DEM with uncorrelated PyMatching to label each shot. Neither feature
preparation nor network inference calls a decoder.

The learning curves use nested prefixes of the same training dataset, with
all feature and target statistics recomputed from that prefix alone:

| Training shots | Complementary gaps <=40 dB | Cluster gaps <=40 dB |
|---:|---:|---:|
| 50,000 | 719 | 811 |
| 200,000 | 2,894 | 3,331 |
| 1,000,000 | 14,990 | 17,157 |

Validation contains 1,483 complementary lows and 1,702 cluster lows. The small
number of low-gap examples in the pilot makes a larger dataset consequential.
Independent splits can contain identical syndromes; repeated physical outcomes
are valid samples from the same noise distribution. This is an IID experiment,
not a split that withholds every previously observed syndrome.
The training set contains 970,829 distinct syndromes. Both gap labels are
identical across every repeated-syndrome group; the audit is recorded in
`duplicate_label_audit.json`.

Two additional controls rank shots by total Hamming weight or memory-component
Hamming weight. Their score direction and validation thresholds are frozen
before the test is generated. Inclusive ties can make attained recall exceed
the nominal target. These controls require only counters at inference.

## Inputs and training

All models are fully connected networks with two 64-unit ReLU hidden layers
and two outputs:

| Input | Values | Parameters |
|---|---:|---:|
| Raw syndrome bits at fixed physical detector positions | 336 | 25,858 |
| Affinity and boundary summaries | 24 | 5,890 |
| Packed affinity upper triangle | 2,352 | 154,882 |
| Hashed affinity upper triangle | 2,352 | 154,882 |

Numeric detector IDs are not input values. Raw bits retain a fixed mapping
from a physical detector to an input position. Both triangles retain the
previous 48-slot design with the two distinct boundary classes fixed at slots
46 and 47, including for empty and even-weight syndromes. They contain path
strength and presence/connectivity channels, including diagonal entries.
The hashed representation permits each detector only its precomputed candidate
slots and never merges or drops events. The full training preparation completed
without overflow; the largest observed syndrome weight was 40.

Triangle arrays use float16 **storage**, then float32 computation during
training and inference. Labels remain float64. This storage choice is common
to all learning-curve sizes. It is not a fixed-point FPGA evaluation. Inputs
use training-prefix moments; triangle standard deviations have a 0.05 floor.

Two training objectives are compared:

* **Regression:** jointly predict both gaps with standardized mean squared
  error, then use negative predicted gap as the low-confidence score.
* **Classification:** jointly predict the two `gap <= 40` labels with weighted
  binary cross entropy. Each output's positive weight is the training-prefix
  count of highs divided by its count of lows. Sigmoid outputs are ranking
  scores, not calibrated posterior probabilities.

AdamW uses learning rate 0.001, weight decay 0.0001, batches of 1,024, four CPU
threads per job, and validation-loss early stopping with patience five. Epoch
caps are 60 for 50k, 40 for 200k, and 25 for 1M. Raw-bit and packed-triangle
models have all three training sizes. Summary and hashed-triangle controls use
1M shots. Raw-bit 1M models also have initialization seeds 202 and 303 in
addition to the primary seed 101. Other configurations use seed 101.
Some runs select their final allowed epoch; the study measures performance
under these training budgets and does not establish convergence or the best
possible result for an architecture.

## How to read the screening metrics

The positive class is **low confidence**. Average precision measures ranking
of this rare class; its uninformed baseline is the low-gap prevalence. A
classifier that always predicts high confidence would already be roughly
98.5% accurate for the complementary target, while missing every low-gap shot.
Overall accuracy is therefore insufficient.

For each trained model, validation alone selects score thresholds targeting
95%, 99%, and 99.9% recall among low-gap shots. These thresholds are saved with
the checkpoint and reused unchanged on test. Reports give the actual test
recall, fraction flagged low, and false-high count. Flagging a shot low could
route it to an expensive decoder, so the flagged fraction describes a possible
workload tradeoff; it is not a measured decoder speedup.

A false-high is a true gap <=40 dB accepted as high confidence. Its rates among
all shots, among low-gap shots, and among accepted shots answer different
questions and are reported separately. Default thresholds are also retained:
0.5 on classification scores and 40 dB on regressed gaps.

The 99.9% operating point is an empirical quantile of approximately 1,500–1,700
validation lows. It cannot establish a universal bound, and may flag nearly
all shots. Neither a good average score nor zero observed misses certifies
all syndromes. Triangle capacity is also finite: a deployment must explicitly
handle overflow rather than assume the observed training maximum is universal.

## Fresh test results

All twenty checkpoints and their validation thresholds were hashed before
generating the new test set. The lock is `million_study/holdout_lock.json`.
The selection rule was maximum mean validation average precision across the
two targets. It selected `regression_syndrome_n1000000_s101`; within direct
classification it selected `classification_syndrome_n1000000_s303`.
The test contains **3,029 complementary lows and 3,445 cluster lows** among
200,000 shots. No test result was used to change a model or threshold.
Of those shots, 192,846 (96.423%) have syndromes absent from training, including
every low-gap test shot. On those novel syndromes alone, selected-model MAE is
6.39/6.19 dB; `novel_syndrome_audit.json` records this descriptive check.

### Learning curve

Raw-bit regression, primary initialization seed 101, evaluated on the same
fresh test set:

| Training shots | Complementary MAE (dB) | Cluster MAE (dB) | Complementary AP | Cluster AP |
|---:|---:|---:|---:|---:|
| 50,000 | 10.35 | 10.23 | 0.171 | 0.211 |
| 200,000 | 7.56 | 7.36 | 0.265 | 0.317 |
| 1,000,000 | 6.23 | 6.04 | 0.345 | 0.412 |

Here AP is average precision for the `gap <=40 dB` class using negative
predicted gap as the score. The improvement supports increasing the pilot
dataset size. At 1M, mean raw-bit MAE over three initializations is 6.11/5.92
dB (complementary/cluster), so the improvement is not confined to one seed.

### Inputs at one million training shots

Regression MAE on fresh test:

| Input | Seeds | Complementary MAE (dB) | Cluster MAE (dB) |
|---|---:|---:|---:|
| Raw syndrome bits | 3 | 6.11 | 5.92 |
| Affinity and boundary summaries | 1 | 8.96 | 8.83 |
| Packed upper triangle | 1 | 10.18 | 10.02 |
| Hashed upper triangle | 1 | 11.79 | 11.74 |

This network/training setup favors raw bits with stable physical placement.
The tested hash layout does not improve prediction. This comparison does not
rule out better affinity architectures, layouts, or optimization schedules.

### The 40 dB screening decision

At the operating point targeting **99% validation recall**, with every cutoff
applied unchanged to test:

| Score | Test AP, complementary / cluster | Test recall, complementary / cluster | Flagged fraction, complementary / cluster | False highs, complementary / cluster |
|---|---:|---:|---:|---:|
| Selected raw-bit regressor | 0.345 / 0.412 | 98.78% / 98.90% | 40.58% / 40.68% | 37 / 38 |
| Selected direct binary network | 0.175 / 0.218 | 98.94% / 98.90% | 62.88% / 62.06% | 32 / 38 |
| Memory-component weight | 0.074 / 0.089 | 99.04% / 99.10% | 79.54% / 79.54% | 29 / 31 |
| Total weight | 0.075 / 0.088 | 99.54% / 98.98% | 87.56% / 80.99% | 14 / 35 |

These are matched **validation targets**, not exactly matched attained test
recalls; discrete weight-score ties especially change attained recall.
For this recipe, regression retains useful gap information and yields a
better rare-class ranking than direct binary training. A binary output remains
useful as a routing decision, but discarding real-valued gap supervision did
not help these experiments. This is not a general claim that regression must
always outperform classification.

For the selected regressor, the validation thresholds correspond to predicted
gaps of **90.16 dB / 89.27 dB**, not 40 dB. Thresholding its predictions at
40 dB directly catches only **19.12% / 26.12%** of true lows. Small overall
MAE does not ensure good decisions in the rare low-gap tail.

At the more conservative **99.9% validation target**, the selected regressor
flags 66.05% of shots and attains test recall **99.802% / 99.826%**, with
**six false highs for each target**. These are the same six syndromes, all
absent from the exact training sample. Their true complementary gaps range
from 21.42 to 38.95 dB, while predictions range from 102.29 to 111.20 dB.
Errors therefore include large mistakes, not only cases close to 40 dB.
`selected_error_analysis.json` records their rows, detector sets, predictions,
and training-occurrence counts. No all-syndrome error guarantee is established.

The generated `million_study/final_report_v4/` contains CSV metrics, learning curves,
operating-point plots, and a report for all twenty models. Its learning-curve
plot averages both targets and available seeds; the table above fixes seed
101 to make the data-size comparison explicit. Weight baselines are recorded
in `weight_baseline_test.json`. The training-only plot under
`training_distribution/` shows class prevalence by syndrome weight; estimates
at the largest weights are noisy because there are very few such samples.

## Reproduction and artifacts

The new preparation and training entry points are
[`prepare_gap_learning.py`](../scripts/prepare_gap_learning.py) and
[`train_gap_confidence.py`](../scripts/train_gap_confidence.py). Use the Python
dependencies and local environment setup from the initial study. Generated
datasets, logs, checkpoints, predictions, and reports stay under
`.agentic_artifacts/affinity_gap/million_study/`.

For example, after building `export_affinity_dataset`:

```sh
.agentic_artifacts/affinity_gap/build/export_affinity_dataset \
    --distance 7 --shots 1000000 --seed 8101 \
    --output-dir .agentic_artifacts/affinity_gap/million_study/train_raw
python3 scripts/prepare_gap_learning.py \
    --dataset .agentic_artifacts/affinity_gap/million_study/train_raw \
    --output .agentic_artifacts/affinity_gap/million_study/train_prepared \
    --workers 8
python3 scripts/train_gap_confidence.py train \
    --train .agentic_artifacts/affinity_gap/million_study/train_prepared \
    --validation .agentic_artifacts/affinity_gap/million_study/validation_prepared \
    --feature syndrome --objective classification --training-shots 1000000 \
    --seed 101 --threads 4 --epochs 25 --patience 5 --threshold-db 40 \
    --output .agentic_artifacts/affinity_gap/million_study/runs/classification_syndrome_n1000000_s101
```

Generate and prepare the validation split analogously with 100,000 shots and
seed 8102. All output paths must be new. To evaluate frozen checkpoints:

```sh
python3 scripts/train_gap_confidence.py evaluate \
    --data .agentic_artifacts/affinity_gap/million_study/test_prepared \
    --runs .agentic_artifacts/affinity_gap/million_study/runs/* \
    --threads 4 \
    --output .agentic_artifacts/affinity_gap/million_study/final_results
```

Protocol checks require disjoint sampling seeds, matching static circuit/DEM
and topology hashes, and consistent feature layout and normalization. Scalers,
class weights, and target statistics use training data only. The optimized
preparation matched the earlier feature encoder after storage conversion on
an independent sample, and both objectives completed end-to-end smoke runs on
separate small training, validation, and test datasets before the main jobs.

## FPGA considerations

The raw-bit MLP is a useful implementation candidate because its input needs
only fixed detector wiring and its operations are dot products and ReLUs.
Input and output affine scaling can be folded into weights and biases.
Classification can compare logits against transformed thresholds without a
sigmoid circuit. At 16 bits per parameter, raw-bit weights and biases occupy
about 51.7 kB, versus 11.8 kB for the summary model and 309.8 kB for either
triangle model. These are storage estimates, not synthesis or numerical-error
results. Quantization, FPGA resource use, and behavior across code distances
remain to be evaluated.

## A separate route to an all-syndrome certificate

The learned scores do not provide a universal guarantee. For the unsigned
complementary gap, there is a useful conditional bound worth implementing
separately. In the same nonnegative weighted graph used by the decoder, let
`L` be the minimum cost of a detector-free chain with nontrivial logical
parity, and let `U(s)` upper-bound the cost of any feasible correction for
syndrome `s`. Then

```text
complementary_gap(s) >= max(0, L - 2*U(s)) * decibels_per_weight
```

To see this, take a feasible correction E of cost at most U and any correction
F in the opposite logical class. Their symmetric difference has nontrivial
logical parity, so `L <= cost(E xor F) <= cost(E) + cost(F)`. The minimum
opposite-class cost is therefore at least `L-U`, whereas E's class has optimum
at most U. If `L>2U`, the classes are separated by at least `L-2U`; otherwise
zero is the valid lower bound.

A cheap U can sum precomputed costs of paths from individual defects to
physical boundaries: XORing the paths gives a feasible correction whose cost
cannot exceed their sum. No per-shot matching is needed. This may be a weak
bound, but can certify a subset of high-gap syndromes. Distinct logical
boundary classes are essential. A strict implementation must use the same
integer-quantized graph weights and dB scale as complementary decoding;
unrounded affinity costs do not automatically certify its labels. This
argument does not bound the decoder-specific cluster gap. The bound is a
theoretical next step and is not included in the measured neural results.
