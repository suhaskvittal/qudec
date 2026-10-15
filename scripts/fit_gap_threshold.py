#!/usr/bin/env python3
"""Fit log10 P(logical error | gap >= threshold) from a gap-study run."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def read_bins(path, shots_name, errors_name):
    histograms = {shots_name: {}, errors_name: {}}
    with path.open(newline="") as source:
        reader = csv.DictReader(source)
        required = {"histogram", "bin_low_db", "bin_high_db", "extreme", "edge_db", "count"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path.name} has no per-gap error histogram; rerun with an updated collector")
        for row in reader:
            name = row["histogram"]
            if name == "SIGNED_GAP" and shots_name == "UNSIGNED_GAP":
                continue
            if name not in histograms:
                raise ValueError(f"Unknown histogram {name!r}")
            count = int(row["count"])
            if count < 0:
                raise ValueError("Histogram counts must be nonnegative")
            if row["extreme"]:
                if row["extreme"] != ">=":
                    raise ValueError("Gap underflow cannot be used for threshold fitting")
                key = float(row["edge_db"])
            else:
                key = float(row["bin_low_db"])
                high = float(row["bin_high_db"])
                if not math.isfinite(high) or high <= key:
                    raise ValueError("Invalid gap bin edges")
            if not math.isfinite(key) or key < 0 or key in histograms[name]:
                raise ValueError("Invalid or duplicate gap threshold")
            histograms[name][key] = count
    shots = histograms[shots_name]
    errors = histograms[errors_name]
    if not shots:
        raise ValueError(f"No {shots_name} bins found")
    if set(errors) - set(shots):
        raise ValueError("Error histogram contains a bin absent from the shot histogram")
    if any(errors.get(edge, 0) > count for edge, count in shots.items()):
        raise ValueError("Error count exceeds shot count in a gap bin")
    return shots, errors


def gap_format(folder):
    cluster = folder / "cluster_gap_bins.csv"
    complementary = folder / "gap_bins.csv"
    if cluster.is_file() and complementary.is_file():
        raise ValueError("Both gap CSV formats are present; input folder is ambiguous")
    if cluster.is_file():
        return "cluster", cluster, "CLUSTER_GAP_DB", "CLUSTER_GAP_DB_ERRORS"
    if complementary.is_file():
        return "complementary", complementary, "UNSIGNED_GAP", "UNSIGNED_GAP_ERRORS"
    raise ValueError("No cluster_gap_bins.csv or gap_bins.csv found in input folder")


def threshold_rows(shots, errors, min_shots, min_errors):
    surviving_shots = surviving_errors = 0
    rows = []
    for threshold in sorted(shots, reverse=True):
        surviving_shots += shots[threshold]
        surviving_errors += errors.get(threshold, 0)
        probability = surviving_errors / surviving_shots
        included = surviving_shots >= min_shots and surviving_errors >= min_errors
        rows.append({"threshold_db": threshold, "shots_ge_threshold": surviving_shots,
                     "errors_ge_threshold": surviving_errors, "conditional_error_probability": probability,
                     "log10_conditional_error_probability": math.log10(probability) if probability else None,
                     "used_in_fit": included})
    return list(reversed(rows))


def fit_line(rows):
    selected = [row for row in rows if row["used_in_fit"]]
    if len(selected) < 2:
        raise ValueError("Fewer than two usable thresholds; reduce --min-errors/--min-shots or collect more logical errors")
    xs = [row["threshold_db"] for row in selected]
    ys = [row["log10_conditional_error_probability"] for row in selected]
    x_mean = math.fsum(xs) / len(xs)
    y_mean = math.fsum(ys) / len(ys)
    xx = math.fsum((x - x_mean) ** 2 for x in xs)
    if xx == 0:
        raise ValueError("Usable thresholds have no gap range")
    slope = math.fsum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / xx
    intercept = y_mean - slope * x_mean
    total = math.fsum((y - y_mean) ** 2 for y in ys)
    residual = math.fsum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    return {"intercept": intercept, "slope_per_db": slope,
            "r_squared": 1 - residual / total if total else None,
            "thresholds_fitted": len(selected), "fit_threshold_min_db": min(xs),
            "fit_threshold_max_db": max(xs)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_folder", type=Path, help="cluster or complementary gap-study folder")
    parser.add_argument("--min-errors", type=positive, default=5,
                        help="minimum errors above a threshold for fitting (default: 5)")
    parser.add_argument("--min-shots", type=positive, default=100,
                        help="minimum shots above a threshold for fitting (default: 100)")
    args = parser.parse_args()
    folder = args.input_folder.expanduser().resolve()
    try:
        with (folder / "summary.json").open() as source:
            summary = json.load(source)
        gap_type, bins_path, shots_name, errors_name = gap_format(folder)
        shots, errors = read_bins(bins_path, shots_name, errors_name)
        if sum(shots.values()) != summary.get("shots_observed") or sum(errors.values()) != summary.get("logical_errors_observed"):
            raise ValueError("Histogram totals disagree with summary.json")
        rows = threshold_rows(shots, errors, args.min_shots, args.min_errors)
        fit = fit_line(rows)
        result = {"input_folder": str(folder), "gap_type": gap_type,
                  "quantity": f"P(logical error | {gap_type} gap >= G)",
                  "gap_units": "dB", "model": "log10(probability) = intercept + slope_per_db * G",
                  "min_errors": args.min_errors, "min_shots": args.min_shots,
                  "shots_observed": sum(shots.values()), "logical_errors_observed": sum(errors.values()), **fit}
        with (folder / "gap_thresholds.csv").open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        (folder / "gap_threshold_fit.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as error:
        parser.exit(1, f"fit_gap_threshold: {error}\n")


if __name__ == "__main__":
    sys.exit(main())
