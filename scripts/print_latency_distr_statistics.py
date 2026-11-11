#!/usr/bin/env python3
"""Print fitted latency models from PyMatching and seeded histogram logs.

The lognormal parameters describe latency in microseconds as
exp(N(mu_ln_us, sigma_ln)). The two-normal model is a weighted mixture of
N(mean_us, sd_us) components. Fits use the geometric midpoint of each printed
logarithmic bin, so parameters are approximate. By default, correlated runs
exclude bins above twice the median as likely per-rank startup outliers.
"""

import argparse
import math
from pathlib import Path
import re

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
HEADERS = {
    "pymatching": re.compile(r"^PYMATCHING_DECODE_LATENCY_US\s*=+"),
    "correlated_pymatching": re.compile(r"^PYMATCHING_DECODE_LATENCY_US\s*=+"),
    "seeded_correlated_pymatching": re.compile(r"^SEEDED_SECOND_DECODE_LATENCY_US\s*=+"),
}
SUMMARY = re.compile(r"\bcount\s*=\s*(\d+)\s*$")
BIN = re.compile(r"^([0-9.eE+-]+)\s*<=\s*X\s*<\s*([0-9.eE+-]+):\s*(\d+)\s*$")
EXTREME = re.compile(r"^(?:<|>=)[0-9.eE+-]+:\s*(\d+)\s*$")


def read_histogram(path, header):
    lines = path.read_text().splitlines()
    start = next((i for i, line in enumerate(lines) if header.match(line)), None)
    if start is None:
        raise ValueError(f"{path}: latency histogram is missing or incomplete")
    summary = SUMMARY.search(lines[start + 1]) if start + 1 < len(lines) else None
    if summary is None:
        raise ValueError(f"{path}: histogram count is missing")
    expected = int(summary[1])
    bounds, counts = [], []
    for line in lines[start + 2:]:
        if not line.strip():
            break
        match = BIN.match(line)
        if match:
            lo, hi, count = float(match[1]), float(match[2]), int(match[3])
            if not (0 < lo < hi):
                raise ValueError(f"{path}: invalid latency bin {line!r}")
            bounds.append((lo, hi))
            counts.append(count)
        elif (match := EXTREME.match(line)) and int(match[1]):
            raise ValueError(f"{path}: nonzero underflow/overflow cannot be fitted from bin limits")
        else:
            raise ValueError(f"{path}: unrecognized histogram line {line!r}")
    if not bounds or sum(counts) != expected:
        raise ValueError(f"{path}: histogram bin counts do not match {expected} shots")
    return np.asarray(bounds), np.asarray(counts), expected


def normal_bin_mass(lo, hi, mean, sd):
    """Normal probability in each interval, stable on either side of the mean."""
    root_two = math.sqrt(2.0)
    result = np.empty(len(lo))
    for i, (a, b) in enumerate(zip(lo, hi)):
        if a >= mean:
            result[i] = 0.5 * (math.erfc((a - mean) / (sd * root_two))
                               - math.erfc((b - mean) / (sd * root_two)))
        else:
            result[i] = 0.5 * (math.erfc((mean - b) / (sd * root_two))
                               - math.erfc((mean - a) / (sd * root_two)))
    return result


def lognormal_log_likelihood(params, lo, hi, counts):
    mu, log_sigma = params
    sigma = math.exp(log_sigma)
    masses = normal_bin_mass(np.log(lo), np.log(hi), mu, sigma)
    return float(counts @ np.log(np.maximum(masses, 1e-300)))


def fit_lognormal(lo, hi, counts):
    midpoints = np.sqrt(lo * hi)
    log_values = np.log(np.repeat(midpoints, counts))
    start = np.array([log_values.mean(), math.log(log_values.std())])
    step = np.array([log_values.std() * .1, .1])
    simplex = [start, start + [step[0], 0], start + [0, step[1]]]
    scores = [-lognormal_log_likelihood(p, lo, hi, counts) for p in simplex]
    for _ in range(350):
        order = np.argsort(scores)
        simplex = [simplex[i] for i in order]
        scores = [scores[i] for i in order]
        if max(np.linalg.norm(simplex[i] - simplex[0]) for i in (1, 2)) < 1e-7:
            break
        center = (simplex[0] + simplex[1]) / 2
        reflected = center + (center - simplex[2])
        reflected_score = -lognormal_log_likelihood(reflected, lo, hi, counts)
        if reflected_score < scores[0]:
            expanded = center + 2 * (reflected - center)
            expanded_score = -lognormal_log_likelihood(expanded, lo, hi, counts)
            simplex[2], scores[2] = ((expanded, expanded_score) if expanded_score < reflected_score
                                     else (reflected, reflected_score))
        elif reflected_score < scores[1]:
            simplex[2], scores[2] = reflected, reflected_score
        else:
            contracted = center + .5 * (simplex[2] - center)
            contracted_score = -lognormal_log_likelihood(contracted, lo, hi, counts)
            if contracted_score < scores[2]:
                simplex[2], scores[2] = contracted, contracted_score
            else:
                for i in (1, 2):
                    simplex[i] = simplex[0] + .5 * (simplex[i] - simplex[0])
                    scores[i] = -lognormal_log_likelihood(simplex[i], lo, hi, counts)
    mu, log_sigma = simplex[int(np.argmin(scores))]
    return float(mu), float(math.exp(log_sigma))


def normal_density(x, mean, sd):
    return np.exp(-.5 * ((x - mean) / sd) ** 2) / (sd * math.sqrt(2 * math.pi))


def fit_two_normals(lo, hi, counts):
    """Weighted midpoint EM with several deterministic starts; choose bin likelihood."""
    x = np.sqrt(lo * hi)
    samples = np.repeat(x, counts)
    best = None
    for quantile in (.02, .05, .1, .2, .35, .999):
        if quantile == .999:
            mean1, mean2 = np.quantile(samples, [.5, .999])
            sd1 = samples.std() * .5
            sd2 = max(1.0, samples[samples > np.quantile(samples, .995)].std())
            weight = .998
        else:
            mean1, mean2 = np.quantile(samples, [quantile, .7])
            sd1 = sd2 = samples.std() * .5
            weight = .2
        for _ in range(500):
            a = weight * normal_density(x, mean1, sd1)
            b = (1 - weight) * normal_density(x, mean2, sd2)
            responsibility = a / np.maximum(a + b, 1e-300)
            n1 = counts @ responsibility
            n2 = counts.sum() - n1
            if min(n1, n2) < 5:
                break
            weight = n1 / counts.sum()
            mean1 = (counts * responsibility) @ x / n1
            mean2 = (counts * (1 - responsibility)) @ x / n2
            sd1 = max(1.0, math.sqrt((counts * responsibility) @ ((x - mean1) ** 2) / n1))
            sd2 = max(1.0, math.sqrt((counts * (1 - responsibility)) @ ((x - mean2) ** 2) / n2))
        masses = (weight * normal_bin_mass(lo, hi, mean1, sd1)
                  + (1 - weight) * normal_bin_mass(lo, hi, mean2, sd2))
        score = float(counts @ np.log(np.maximum(masses, 1e-300)))
        if best is None or score > best[0]:
            best = (score, weight, mean1, sd1, mean2, sd2)
    _, weight, mean1, sd1, mean2, sd2 = best
    components = sorted(((weight, mean1, sd1), (1 - weight, mean2, sd2)), key=lambda x: x[1])
    return components


def print_statistics(root, seeded_root, include_extremes):
    print("Latency units: microseconds (us)")
    print("Lognormal: T = exp(N(mu_us, sigma_ln^2))")
    print("Two-normal mixture: select a component by weight, then sample N(mean_us, sd_us^2)")
    found = False
    datasets = (
        ("pymatching", root / "pymatching_latency", "PyMatching"),
        ("correlated_pymatching", root / "correlated_pymatching_latency", "Correlated PyMatching"),
        ("seeded_correlated_pymatching", seeded_root / "seeded_correlated_pymatching_latency",
         "Seeded correlated PyMatching second decoder"),
    )
    for decoder, directory, label in datasets:
        for path in sorted(directory.glob("d*.out"), key=lambda p: int(p.stem[1:])):
            if not re.fullmatch(r"d\d+", path.stem):
                continue
            bounds, counts, total = read_histogram(path, HEADERS[decoder])
            lo, hi = bounds[:, 0], bounds[:, 1]
            excluded = 0
            if decoder in ("correlated_pymatching", "seeded_correlated_pymatching") and not include_extremes:
                median = np.median(np.repeat(np.sqrt(lo * hi), counts))
                keep = np.sqrt(lo * hi) <= 2 * median
                excluded = int(counts[~keep].sum())
                lo, hi, counts = lo[keep], hi[keep], counts[keep]
            if counts.sum() < 10:
                raise ValueError(f"{path}: too few shots remain to fit distributions")
            log_mu, log_sigma = fit_lognormal(lo, hi, counts)
            (w1, mean1, sd1), (w2, mean2, sd2) = fit_two_normals(lo, hi, counts)
            print(f"\n{label}, d={path.stem[1:]} ({int(counts.sum()):,} of {total:,} shots fitted)")
            if excluded:
                print(f"  Excluded {excluded} shot(s) above twice the median (startup outlier heuristic)")
            print(f"  Lognormal: mu_us={np.exp(log_mu):.8f}, sigma_ln={log_sigma:.8f}")
            print("  Two-normal mixture:")
            print(f"    Component 1: weight={w1:.6f}, mean_us={mean1:.2f}, sd_us={sd1:.2f}")
            print(f"    Component 2: weight={w2:.6f}, mean_us={mean2:.2f}, sd_us={sd2:.2f}")
            found = True
    if not found:
        raise ValueError(f"No latency histograms found under {root} or {seeded_root}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=REPOSITORY / "out" / "single_core",
                        help="directory containing pymatching_latency and correlated_pymatching_latency")
    parser.add_argument("--seeded-input-root", type=Path, default=REPOSITORY / "out",
                        help="directory containing seeded_correlated_pymatching_latency")
    parser.add_argument("--include-extremes", action="store_true",
                        help="fit correlated startup outliers too (default: exclude bins above 2x median)")
    args = parser.parse_args()
    try:
        print_statistics(args.input_root.expanduser().resolve(),
                         args.seeded_input_root.expanduser().resolve(), args.include_extremes)
    except (OSError, ValueError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")


if __name__ == "__main__":
    main()
