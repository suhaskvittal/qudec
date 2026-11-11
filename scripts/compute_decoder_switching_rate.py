#!/usr/bin/env python3
"""Estimate UF-to-L2 and PyMatching-to-L2 switching from cluster-gap histograms.

Fit a separate binomial logistic model for P(L1 logical error | gap in dB) using
the d=5,7,9,11 shot and error histograms for each L1. Shift each fitted curve
to match a supplied d=31 per-round L1 error rate multiplied by 93 rounds
(3d) before finding the cutoff. For a target residual L1 error rate, choose
the smallest gap cutoff that satisfies

    sum_{gap > cutoff} P_d31(gap) P(L1 error | gap) <= target.

Both L1 estimates use the d=31 PyMatching gap frequencies. This treats L2
corrections as error-free; its own error contribution must be accounted for
separately. Printed histogram bins are one dB wide, so cutoffs are bin edges.
By default, a Gaussian approximation inferred from the d=31 histogram moments
resolves its >=256 dB overflow. Use --overflow-model bound to leave it unresolved.
"""

import argparse
from collections import defaultdict
import math
from pathlib import Path

from run_cluster_gap_study import parse_histograms


REPOSITORY = Path(__file__).resolve().parents[1]
TRAIN_DISTANCES = (5, 7, 9, 11)
SHOTS = "CLUSTER_GAP_DB"
ERRORS = "CLUSTER_GAP_DB_ERRORS"


def probability(z):
    if z >= 0:
        return 1 / (1 + math.exp(-z))
    exp_z = math.exp(z)
    return exp_z / (1 + exp_z)


def histogram_rows(path):
    histograms = parse_histograms(path)
    shots, errors = histograms[SHOTS], histograms[ERRORS]
    shot_bins = shots["bins"]
    error_bins = errors["bins"]
    if set(error_bins) - set(shot_bins):
        raise ValueError(f"{path}: error bins have no matching shot bins")
    if any(error_bins.get(edge, 0) > count for edge, count in shot_bins.items()):
        raise ValueError(f"{path}: errors exceed shots in a gap bin")
    if any(lo < 0 or hi <= lo for lo, hi in shot_bins):
        raise ValueError(f"{path}: invalid gap bin edges")
    if any(direction != ">=" for direction, _, count in shots["extremes"] if count):
        raise ValueError(f"{path}: gap underflow is unsupported")
    if errors["extremes"]:
        raise ValueError(f"{path}: errors outside the fitted gap bins are unsupported")
    overflow = [(edge, count) for direction, edge, count in shots["extremes"] if count]
    if len(overflow) > 1:
        raise ValueError(f"{path}: multiple overflow bins")
    rows = sorted((lo, hi, n, error_bins.get((lo, hi), 0))
                  for (lo, hi), n in shot_bins.items() if n)
    if any(right[0] < left[1] for left, right in zip(rows, rows[1:])):
        raise ValueError(f"{path}: overlapping gap bins")
    return rows, overflow, shots["count"], errors["count"]


def normal_interval(lo, hi, mean, sd):
    root_two = math.sqrt(2)
    if lo >= mean:
        return (math.erfc((lo - mean) / (sd * root_two))
                - math.erfc((hi - mean) / (sd * root_two))) / 2
    return (math.erfc((mean - hi) / (sd * root_two))
            - math.erfc((mean - lo) / (sd * root_two))) / 2


def expand_overflow_normal(path, rows, overflow, total_shots):
    """Estimate overflow bin frequencies from its first two moments."""
    if not overflow:
        return rows, None
    edge, overflow_shots = overflow[0]
    summary = parse_histograms(path)[SHOTS]
    mean, sd, maximum = summary["mean_db"], summary["std_db"], summary["max_db"]
    if any(value is None or not math.isfinite(value) for value in (mean, sd, maximum)):
        raise ValueError(f"{path}: missing d=31 gap moments")
    finite_sum = sum(n * (lo + hi) / 2 for lo, hi, n, _ in rows)
    finite_sq_sum = sum(n * ((lo + hi) / 2)**2 for lo, hi, n, _ in rows)
    overflow_mean = (total_shots * mean - finite_sum) / overflow_shots
    overflow_sq = (total_shots * (mean * mean + sd * sd) - finite_sq_sum) / overflow_shots
    overflow_sd = math.sqrt(max(0, overflow_sq - overflow_mean * overflow_mean))
    if overflow_sd <= 0 or maximum <= edge:
        raise ValueError(f"{path}: invalid overflow moments")
    upper = math.ceil(maximum)
    bins = [(float(lo), float(lo + 1)) for lo in range(int(edge), upper)]
    weights = [normal_interval(lo, hi, overflow_mean, overflow_sd) for lo, hi in bins]
    total_weight = sum(weights)
    if not total_weight:
        raise ValueError(f"{path}: overflow Gaussian has zero probability in range")
    expanded = [(lo, hi, overflow_shots * weight / total_weight, 0)
                for (lo, hi), weight in zip(bins, weights)]
    return rows + expanded, (overflow_mean, overflow_sd, edge, upper)


def log_likelihood(rows, intercept, slope):
    result = 0.0
    for gap, shots, errors in rows:
        z = intercept + slope * gap
        result += errors * (-math.log1p(math.exp(-z)) if z >= 0 else z - math.log1p(math.exp(z)))
        result += (shots - errors) * (-z - math.log1p(math.exp(-z)) if z >= 0 else -math.log1p(math.exp(z)))
    return result


def fit_error_curve(rows):
    """Two-parameter binomial maximum likelihood fit via damped Newton steps."""
    intercept, slope = 0.0, -0.2
    likelihood = log_likelihood(rows, intercept, slope)
    for _ in range(100):
        gradient_0 = gradient_1 = h00 = h01 = h11 = 0.0
        for gap, shots, errors in rows:
            p = probability(intercept + slope * gap)
            residual = errors - shots * p
            weight = shots * p * (1 - p)
            gradient_0 += residual
            gradient_1 += residual * gap
            h00 += weight
            h01 += weight * gap
            h11 += weight * gap * gap
        determinant = h00 * h11 - h01 * h01
        if determinant <= 0:
            raise ValueError("Cannot fit gap error curve: singular information matrix")
        step_0 = (h11 * gradient_0 - h01 * gradient_1) / determinant
        step_1 = (h00 * gradient_1 - h01 * gradient_0) / determinant
        fraction = 1.0
        while fraction >= 2**-30:
            trial_0 = intercept + fraction * step_0
            trial_1 = slope + fraction * step_1
            trial_likelihood = log_likelihood(rows, trial_0, trial_1)
            if trial_likelihood >= likelihood - 1e-9:
                break
            fraction *= 0.5
        if fraction < 2**-30:
            raise ValueError("Gap error curve fit failed to converge")
        intercept, slope, likelihood = trial_0, trial_1, trial_likelihood
        if max(abs(fraction * step_0), abs(fraction * step_1)) < 1e-10:
            break
    else:
        raise ValueError("Gap error curve fit exceeded iteration limit")
    if slope >= 0:
        raise ValueError("Fitted error probability does not decrease with gap")
    return intercept, slope


def modeled_error_rate(rows, overflow, total_shots, intercept, slope):
    risk = sum(n * probability(intercept + slope * lo) for lo, _, n, _ in rows)
    risk += sum(n * probability(intercept + slope * edge) for edge, n in overflow)
    return risk / total_shots


def calibrate_intercept(rows, overflow, total_shots, intercept, slope, baseline):
    """Preserve the fitted gap slope while matching the supplied d=31 LER."""
    low, high = intercept - 1000, intercept + 1000
    for _ in range(200):
        middle = (low + high) / 2
        if modeled_error_rate(rows, overflow, total_shots, middle, slope) < baseline:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def switching_threshold(rows, overflow, total_shots, intercept, slope, target):
    """Return cutoff, switched shots and conservative remaining L1 error rate."""
    # The lower bin edge upper-bounds risk throughout each finite bin because
    # the fitted error probability decreases with gap.
    residual = modeled_error_rate(rows, overflow, total_shots, intercept, slope)
    if residual <= target:
        return None, 0, residual
    switched = 0
    for lo, hi, shots, _ in rows:
        switched += shots
        residual -= shots * probability(intercept + slope * lo) / total_shots
        if residual <= target:
            return hi, switched, max(residual, 0.0)
    if overflow:
        edge, shots = overflow[0]
        # All overflow shots have gap >= edge; their internal gap distribution
        # is unknown, so only an all-L2 estimate is available at this point.
        return math.inf, switched + shots, 0.0
    return rows[-1][1], switched, 0.0


def positive_rate(value):
    try:
        rate = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a finite rate in (0, 1)") from error
    if not math.isfinite(rate) or not 0 < rate < 1:
        raise argparse.ArgumentTypeError("must be a finite rate in (0, 1)")
    return rate


def positive_int(value):
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a positive integer") from error
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-error-rate", type=positive_rate, required=True,
                        help="maximum modeled residual L1 logical errors per d=31 shot")
    parser.add_argument("--uf-ler-per-round", "--uf-ler", dest="uf_ler_per_round",
                        type=positive_rate, default=2.2e-14,
                        help="projected d=31 UF LER per round (default: 2.2e-14)")
    parser.add_argument("--pymatching-ler-per-round", "--pymatching-ler",
                        dest="pymatching_ler_per_round", type=positive_rate, default=8.0e-15,
                        help="projected d=31 PyMatching LER per round (default: 8.0e-15)")
    parser.add_argument("--ler-rounds", type=positive_int, default=93,
                        help="rounds used to convert per-round LER to per-shot LER (default: 93 = 3*31)")
    parser.add_argument("--overflow-model", choices=("normal", "bound"), default="normal",
                        help="estimate d=31 >=256 dB gaps with a Gaussian or retain only bounds")
    parser.add_argument("--input-root", type=Path, default=REPOSITORY / "out",
                        help="folder containing union_find_gap and pymatching_gap")
    args = parser.parse_args()
    if max(args.uf_ler_per_round, args.pymatching_ler_per_round) * args.ler_rounds >= 1:
        parser.error("per-round LER times --ler-rounds must be below 1")
    root = args.input_root.expanduser().resolve()
    try:
        d31_path = root / "pymatching_gap" / "d31.out"
        rows31, overflow31, total31, _ = histogram_rows(d31_path)
        if not total31:
            raise ValueError("d31 histogram contains no shots")
        modeled_overflow = None
        if args.overflow_model == "normal":
            rows31, modeled_overflow = expand_overflow_normal(
                d31_path, rows31, overflow31, total31)
            overflow31_for_fit = []
        else:
            overflow31_for_fit = overflow31
        results = []
        for label, folder, per_round_ler in (("UF", "union_find_gap", args.uf_ler_per_round),
                                             ("PyMatching", "pymatching_gap", args.pymatching_ler_per_round)):
            baseline = per_round_ler * args.ler_rounds
            pooled = defaultdict(lambda: [0, 0])
            observed = []
            training = []
            for distance in TRAIN_DISTANCES:
                rows, overflow, shots, errors = histogram_rows(root / folder / f"d{distance}.out")
                if overflow:
                    raise ValueError(f"{folder}/d{distance}: training histogram has an overflow bin")
                for lo, hi, count, bad in rows:
                    pooled[(lo + hi) / 2][0] += count
                    pooled[(lo + hi) / 2][1] += bad
                    if bad:
                        observed.append(hi)
                training.append((distance, shots, errors))
            if not observed:
                raise ValueError(f"{label}: no training errors to fit")
            fit_rows = sorted((gap, pair[0], pair[1]) for gap, pair in pooled.items())
            intercept, slope = fit_error_curve(fit_rows)
            raw_rate = modeled_error_rate(rows31, overflow31_for_fit, total31, intercept, slope)
            calibrated_intercept = calibrate_intercept(
                rows31, overflow31_for_fit, total31, intercept, slope, baseline)
            cutoff, switched, residual = switching_threshold(
                rows31, overflow31_for_fit, total31, calibrated_intercept, slope,
                args.target_error_rate)
            results.append((label, training, max(observed), intercept, slope,
                            raw_rate, per_round_ler, baseline, calibrated_intercept,
                            cutoff, switched, residual))
    except (OSError, ValueError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")

    print(f"Target residual L1 error rate: {args.target_error_rate:.8g}")
    print(f"LER conversion: per-round rate * {args.ler_rounds} rounds")
    print("L2 decoding errors are excluded from this estimate.")
    print(f"d=31 PyMatching gap histogram: {total31:,} shots (used for both L1 decoders)")
    if overflow31:
        edge, count = overflow31[0]
        print(f"d=31 >= {edge:g} dB overflow: {count:,} shots")
        if modeled_overflow:
            mean, sd, _, upper = modeled_overflow
            print(f"Overflow estimate: Gaussian over [{edge:g}, {upper:g}) dB,"
                  f" inferred mean={mean:.5g}, sd={sd:.5g} dB")
        else:
            print(f"Overflow bound: all {count:,} shots evaluated at {edge:g} dB")
    for (label, training, highest_error_gap, intercept, slope, raw_rate, per_round_ler, baseline,
         calibrated_intercept, cutoff, switched, residual) in results:
        print(f"\n{label} as L1:")
        for distance, shots, errors in training:
            print(f"  Training d={distance}: {shots:,} shots, {errors:,} L1 errors")
        print(f"  Small-distance fit: P(error | gap g dB) = sigmoid({intercept:.8g} {slope:+.8g} * g)")
        print(f"  Largest gap bin with a training error: <{highest_error_gap:g} dB")
        print(f"  Uncalibrated d=31 LER estimate: {raw_rate:.8g}")
        print(f"  Supplied d=31 LER per round: {per_round_ler:.8g}")
        print(f"  Baseline d=31 LER per shot: {baseline:.8g}")
        print(f"  Calibrated d=31 curve: P(error | gap g dB) = sigmoid({calibrated_intercept:.8g} {slope:+.8g} * g)")
        if cutoff is None:
            print("  Switching cutoff: below all observed d=31 gaps (no switching)")
        elif math.isinf(cutoff):
            print(f"  Switching cutoff: unresolved within >= {overflow31[0][0]:g} dB overflow")
            print(f"  Switching rate: between {(switched - overflow31[0][1]) / total31:.8g} and 1;"
                  " all-L2 gives the zero-residual bound")
        else:
            print(f"  Switching cutoff: {cutoff:g} dB (switch bins with gap < {cutoff:g} dB)")
        if cutoff is None or not math.isinf(cutoff):
            print(f"  Estimated switching rate: {switched / total31:.8g}"
                  f" (~{switched:,.0f}/{total31:,} shots)")
        print(f"  Modeled residual L1 error rate: {residual:.8g}")
        print("  d=31 error risk extrapolates beyond training gaps with observed errors.")
    if modeled_overflow:
        print("Overflow switching rates depend on the Gaussian approximation;"
              " rerun the d=31 histogram with a wider range to measure them.")


if __name__ == "__main__":
    main()
