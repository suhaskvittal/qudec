#!/usr/bin/env python3
"""Fit the per-gap logical error rate and report the complementary gap threshold.

Each ``out/gap/d{d}.out`` stats file ends with two unsigned-gap histograms:

    UNSIGNED_GAP         -- gap distribution over ALL shots
    UNSIGNED_GAP_ERRORS  -- gap distribution over shots that produced a logical error

Both are binned by integer gap ``g`` (lines of the form ``g <= X < g+1:\t<count>``),
so the per-gap conditional logical error rate is

    P(error | g) = UNSIGNED_GAP_ERRORS[g] / UNSIGNED_GAP[g]

We expect this to fall off linearly in log space (~10^-g), so we fit

    log10(P(error | g)) = m*g + b

per code distance and pooled across all distances, then invert the fit to find the
complementary gap ``g_th`` at which the per-shot error rate reaches a target value.
"""

import re
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DISTANCES = [7, 9, 11]                                                        
TARGET_LOGICAL_ERROR_RATE = 23 * 1e-12 * 0.01
MIN_ERRORS = 5                          # drop gap bins with fewer errors than this from the fit
OUT_DIR = "out/gap"                     # resolved relative to the repo root

REPO_ROOT = Path(__file__).resolve().parent.parent
GAP_DIR = REPO_ROOT / OUT_DIR

_BIN_RE = re.compile(r"^(\d+)\s*<=\s*X\s*<\s*\d+:\s*(\d+)\s*$")
_TOTAL_HEADER = "UNSIGNED_GAP ="
_ERRORS_HEADER = "UNSIGNED_GAP_ERRORS ="


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
def _parse_section(lines):
    """Parse ``g <= X < g+1: count`` bins from a block of lines into {g: count}."""
    counts = {}
    for line in lines:
        m = _BIN_RE.match(line)
        if m:
            counts[int(m.group(1))] = int(m.group(2))
    return counts


def parse_unsigned_histograms(path):
    """Return (total_counts, error_counts) dicts keyed by integer gap.

    Raises ValueError if the file does not yet contain both unsigned-gap
    histograms (e.g. a run that is still in progress).
    """
    text = path.read_text().splitlines()

    total_start = errors_start = None
    for i, line in enumerate(text):
        if _TOTAL_HEADER in line:
            total_start = i
        elif _ERRORS_HEADER in line:
            errors_start = i

    if total_start is None or errors_start is None:
        raise ValueError(
            f"{path} is missing the unsigned-gap histograms "
            "(is the Monte Carlo run still in progress?)"
        )

    total_counts = _parse_section(text[total_start + 1:errors_start])
    error_counts = _parse_section(text[errors_start + 1:])
    return total_counts, error_counts


# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------
def fit_points(total, errors):
    """Build (x, y, w) fit arrays from the total / error histograms.

    x = gap bin center, y = log10 P(error | g), w = error count (Poisson weight).
    Bins with too few errors or no samples are dropped.
    """
    xs, ys, ws = [], [], []
    for g, n_err in sorted(errors.items()):
        n_tot = total.get(g, 0)
        if n_tot <= 0 or n_err < MIN_ERRORS:
            continue
        p = n_err / n_tot
        xs.append(g + 0.5)
        ys.append(np.log10(p))
        ws.append(n_err)
    return np.array(xs), np.array(ys), np.array(ws, dtype=float)


def weighted_linear_fit(x, y, w):
    """Weighted least-squares fit y = m*x + b. Returns (m, b, r2)."""
    # np.polyfit weights are applied to the residuals, so pass sqrt of the
    # inverse-variance weights (weight ~ error count).
    m, b = np.polyfit(x, y, 1, w=np.sqrt(w))

    pred = m * x + b
    wmean = np.average(y, weights=w)
    ss_res = np.sum(w * (y - pred) ** 2)
    ss_tot = np.sum(w * (y - wmean) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return m, b, r2


def complementary_gap(m, b):
    """Gap g_th at which the fit predicts P(error | g) == TARGET_LOGICAL_ERROR_RATE."""
    return (np.log10(TARGET_LOGICAL_ERROR_RATE) - b) / m


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def report(label, x, y, w):
    m, b, r2 = weighted_linear_fit(x, y, w)
    g_th = complementary_gap(m, b)
    print(f"[{label}]")
    print(f"  fit points   : {len(x)}")
    print(f"  log10(P)     = {m:.6f} * g + {b:.6f}")
    print(f"  P(error|g)   ~ 10^({m:.6f} * g + {b:.6f})")
    print(f"  weighted R^2 : {r2:.6f}")
    print(f"  g_th         : {g_th:.4f}  (P(error|g_th) = {TARGET_LOGICAL_ERROR_RATE:.3e})")
    print()


def main():
    print(f"TARGET_LOGICAL_ERROR_RATE = {TARGET_LOGICAL_ERROR_RATE:.3e}")
    print(f"MIN_ERRORS per bin        = {MIN_ERRORS}")
    print(f"DISTANCES                 = {DISTANCES}")
    print()

    all_x, all_y, all_w = [], [], []
    for d in DISTANCES:
        path = GAP_DIR / f"d{d}.out"
        total, errors = parse_unsigned_histograms(path)
        x, y, w = fit_points(total, errors)
        if len(x) < 2:
            print(f"[d = {d}] not enough usable bins to fit ({len(x)} points)\n")
            continue
        report(f"d = {d}", x, y, w)
        all_x.append(x)
        all_y.append(y)
        all_w.append(w)

    if all_x:
        x = np.concatenate(all_x)
        y = np.concatenate(all_y)
        w = np.concatenate(all_w)
        report("ALL DISTANCES", x, y, w)


if __name__ == "__main__":
    main()
