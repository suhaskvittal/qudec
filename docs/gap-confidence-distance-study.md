# Comparing pre-decoding gap prediction at distances 7 and 9

Contributor: OpenAI Codex.

This experiment follows the [million-shot d=7 study](gap-confidence-million-study.md).
It asks whether raising the code distance improves prediction with the same
training-shot budget and network family. Models are trained separately for each
distance; this is not transfer of a d=7 checkpoint to d=9.

## Matched configuration

Both distances use rotated surface-code Z memory, SI1000 p=0.001, and both
detector bases. The number of rounds equals the distance. Thus d=7 has seven
rounds and 336 detectors, while d=9 has nine rounds and 720 detectors.

Each distance has 1,000,000 training shots, 100,000 validation shots, and a
fresh 200,000-shot test set. The d=7 sampling seeds are 8101/8102/8103; d=9
uses 9101/9102/9103. Stim samples the undecomposed DEM directly, preserving
correlated-fault supports. Uncorrelated PyMatching on the decomposed DEM
provides unsigned complementary and cluster-gap labels offline. No per-shot
decoding is used to prepare model inputs or make predictions.

The low-confidence class remains `gap <=40 dB`. The 60 dB threshold is included
only in dataset counts, not used to train additional classifiers in this
comparison. A larger distance changes the gap distribution and low-gap
prevalence; an increase in overall accuracy alone would be insufficient
evidence of a better estimator.

The following configurations are matched to the completed d=7 runs:

* Raw syndrome bits at fixed physical detector positions: regression and
  weighted binary classification, each with initialization seeds 101, 202,
  and 303.
* The 24 affinity/boundary summaries: regression and classification, each
  with seed 101.

Each model has two 64-unit ReLU hidden layers and two outputs. Raw-bit models
have 25,858 parameters at d=7 and 50,434 at d=9 because the input width grows.
Summary models have 5,890 parameters at both distances. The hidden width is
fixed; the parameter budget is not identical for raw inputs.

Training uses AdamW, learning rate 0.001, weight decay 0.0001, batches of 1,024,
up to 25 epochs, patience five, and four CPU threads per job. All input and
target statistics and class weights use training data only. Validation loss
selects the epoch. Mean validation average precision across both gap targets
selects a primary screening model and a direct classifier before test
generation. Validation also fixes the thresholds targeting 95%, 99%, and
99.9% low-gap recall.

The matrix models were weaker than these controls at d=7 and are not repeated
here. In particular, the previous 48-slot triangle representation cannot hold
every sampled d=9 syndrome; the initial d=9 validation sample reaches weight
63 and training reaches weight 68. Raw bits and summaries do not truncate
events or impose that matrix limit.

The 1M-shot d=9 training set contains 3,124 complementary lows (0.3124%) and
3,767 cluster lows (0.3767%) at 40 dB. At 60 dB those counts are 17,072
(1.7072%) and 20,546 (2.0546%). None of these labels equals either threshold
exactly. Thus the matched shot budget provides about one fifth as many low-gap
training examples as at d=7. Training contains 999,986 distinct syndromes;
both labels agree exactly within every repeated-syndrome group.

## Comparison criteria

The main operational comparison is actual test low-gap recall against the
fraction flagged for an expensive route, with thresholds fixed by each
distance's validation data. False-high counts, rates among true lows and among
accepted shots, and the lowest true gap mistakenly called high describe the
remaining risk. Default thresholds are included alongside the recall-targeted
thresholds.

Regression MAE and MAE divided by distance are reported separately. The latter
is a relative error measure, not a change to the absolute 40 dB decision.
ROC AUC, average precision, prevalence, and AP/prevalence are included because
the positive-class frequency changes with distance. Total-weight and
memory-component-weight controls use their own validation-frozen thresholds.

The d=9 validation split contains only 335 complementary lows and 396 cluster
lows at 40 dB. Its 99.9% quantile therefore selects the minimum positive score;
the observed validation recall is 100%, without guaranteeing future recall.
The held-out low counts and false-high counts must accompany any comparison
of this operating point. Empirical sampling cannot establish an all-syndrome
guarantee.

## Fresh-test results

All eight d=9 checkpoints and thresholds were hashed before generating test
seed 9103. The validation-selected primary is
`regression_syndrome_n1000000_s101`; the selected direct classifier is
`classification_summary_n1000000_s101`. The corresponding d=7 selections are
raw-bit regression seed 101 and raw-bit classification seed 303. Model choices
and thresholds were not changed after test evaluation.

The 200,000-shot d=9 test contains 618 complementary lows (0.309%) and 726
cluster lows (0.363%), versus 3,029 (1.5145%) and 3,445 (1.7225%) at d=7.
An all-high decision would therefore already be 99.691%/99.637% accurate at
d=9, while missing every low-gap shot. Overall accuracy is unsuitable for
deciding whether the estimator improved.

### Regression accuracy

Raw-bit regressors, mean over initialization seeds 101/202/303:

| Distance | Complementary MAE (dB) | Cluster MAE (dB) | Complementary MAE/d | Cluster MAE/d |
|---:|---:|---:|---:|---:|
| 7 | 6.11 | 5.92 | 0.873 | 0.845 |
| 9 | 10.38 | 10.10 | 1.153 | 1.122 |

Both absolute and distance-normalized errors increase. Mean ROC AUC falls
from 0.944/0.949 to 0.914/0.918. Mean AP falls from 0.342/0.407 to
0.0824/0.1045, while AP/prevalence rises from 22.6/23.6 to 26.7/28.8.
These different movements reinforce the need to inspect actual screening
tradeoffs rather than one aggregate score.

There are 199,998 d=9 test syndromes absent from training. Selected-model MAE
on those novel syndromes is 10.39/10.10 dB. For comparison, selected d=7 MAE
on novel syndromes was 6.39/6.19 dB; the larger d=7 overlap with training does
not account for this error increase. Audits are saved alongside the data.

### Screening at the 99% validation target

Each row uses its distance's validation-selected primary regressor and frozen
thresholds. The target is the same, but attained test recalls are not forced
to match:

| Distance | Target | Test low recall | Fraction flagged low | False highs / true lows |
|---:|---|---:|---:|---:|
| 7 | Complementary | 98.78% | 40.58% | 37 / 3,029 |
| 9 | Complementary | 97.90% | 49.33% | 13 / 618 |
| 7 | Cluster | 98.90% | 40.68% | 38 / 3,445 |
| 9 | Cluster | 97.38% | 47.73% | 19 / 726 |

The selected d=9 regressor flags a larger share and catches a smaller share
of true lows at this operating point. Fewer false highs in absolute count
reflect substantially fewer low-gap test shots; the missed fraction among
true lows is higher. The smaller d=9 low sample also makes its recall
estimate less precise.

At the 99.9% validation target, selected d=9 regression flags 70.50%/70.68%
and misses 3/4 lows, yielding 99.515%/99.449% test recall. The selected d=7
regressor flags 66.05%/66.05%, misses 6/6 lows, and attains
99.802%/99.826% recall.

The selected d=9 **direct classifier** uses summaries. At the 99% validation
target it attains 99.515%/99.587% recall, flags 73.24%/72.93%, and misses
3/3 lows. At its 99.9% validation target it flags 83.35%/83.65% and misses
one low for each target. This stricter screen retains substantial workload
and still makes false-high decisions. The comparison report also includes
all seeds and the matched summary/raw model families, avoiding a conclusion
based only on the two selected runs.

### Lowest gap mistakenly called high

For each validation-selected regressor:

| Validation operating point | d=7 complementary / cluster | d=9 complementary / cluster |
|---|---:|---:|
| Default predicted-gap threshold 40 dB | 0.032 / 0.032 dB | 0.274 / 0.785 dB |
| 95% recall target | 2.309 / 2.309 dB | 2.263 / 2.263 dB |
| 99% recall target | 12.101 / 12.101 dB | 20.617 / 14.099 dB |
| 99.9% recall target | 21.422 / 21.422 dB | 24.790 / 24.790 dB |

The d=9 selected direct classifier's lowest missed gap is 29.881 dB for both
targets at the 99% and 99.9% points, and 0.785 dB at the default and 95%
points. Some observed minima are higher at d=9, but d=9 has far fewer low-gap
test cases; these finite-sample extrema cannot establish a better worst-case
guarantee. Low-gap false highs persist.

### Decision

Increasing distance did not improve this design under the fixed 1M-shot,
64-unit-hidden-layer training setup. Regression error increases, ranking
degrades, and the primary frozen screening policy has worse observed recall
at a larger flagged fraction. Following the requested stopping criterion,
this study ends at d=9; no d=11 sweep or further model search is started.

This conclusion concerns the tested design and budget. It does not establish
that every pre-decoding estimator is ineffective. The rarer positive class,
larger input, and epoch cap may all matter; several regressors select their
last allowed epoch. No universal error guarantee is supplied by these models.

All per-run metrics, routing counts, minima, and plots are in
`.agentic_artifacts/affinity_gap/distance9_study/comparison_report/`, including
`report.md`, `aggregate_metrics.csv`, `routing_comparison.csv`,
`selected_misses.csv`, `mae_distance_comparison.png`, and
`ap_distance_comparison.png`.

## Reproduce

Use the exporter, preparation script, training script, and Python environment
from the d=7 study. All generated d=9 files are under
`.agentic_artifacts/affinity_gap/distance9_study/`.

```sh
.agentic_artifacts/affinity_gap/build/export_affinity_dataset \
    --distance 9 --rounds 9 --p 0.001 --shots 1000000 --seed 9101 \
    --output-dir .agentic_artifacts/affinity_gap/distance9_study/train_raw
python3 scripts/prepare_gap_learning.py \
    --dataset .agentic_artifacts/affinity_gap/distance9_study/train_raw \
    --output .agentic_artifacts/affinity_gap/distance9_study/train_prepared \
    --features syndrome summary --workers 8
python3 scripts/train_gap_confidence.py train \
    --train .agentic_artifacts/affinity_gap/distance9_study/train_prepared \
    --validation .agentic_artifacts/affinity_gap/distance9_study/validation_prepared \
    --feature syndrome --objective regression --training-shots 1000000 \
    --seed 101 --threads 4 --epochs 25 --patience 5 --threshold-db 40 \
    --output .agentic_artifacts/affinity_gap/distance9_study/runs/regression_syndrome_n1000000_s101
```

Generate and prepare validation with 100,000 shots and seed 9102. After all
eight models and their thresholds are frozen, generate and prepare test with
200,000 shots and seed 9103. Output directories must be new. Evaluate with
`train_gap_confidence.py evaluate`, as in the d=7 study. Static graph hashes,
normalization, geometry, and independent source seeds are checked within each
distance. The exploratory `compare_distances.py` in the d=9 artifact directory
generates the comparison report from frozen evaluation outputs.
