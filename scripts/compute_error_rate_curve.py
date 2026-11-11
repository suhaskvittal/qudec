#!/usr/bin/env python3
"""Fit LER = A * (p/B)**floor((d+1)/2) to simulation outputs in a folder.

The .out format used by this repository stores distance and physical error rate
in filenames (d5_p1.00e-01.out), and the logical error rate in the file body.
Run metadata supplies the number of rounds for per-round normalization.
By default, runs with fewer than 10 observed logical errors are excluded.
CSV files with columns d, p, ler (or logical_error_rate), and errors are also accepted.
"""

import argparse
import csv
import json
import math
import re
from pathlib import Path


FILENAME = re.compile(r"d(?P<d>\d+)_p(?P<p>[0-9.eE+-]+)\.out$")
LER_LINE = re.compile(r"^Logical error rate:\s*([0-9.eE+-]+)\s*$", re.MULTILINE)
CYCLES_LINE = re.compile(r"\brounds=(\d+)\b")
RUN_METADATA = re.compile(r"^# BEGIN LER_RUN (\{.*\})$", re.MULTILINE)
PROGRESS_ERRORS = re.compile(r"^\[\d+ / \d+\].*?-- errors = (\d+)\b", re.MULTILINE)


def nonnegative_int(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return number


def rounds_from_output(contents, distance, p, path):
    metadata_match = RUN_METADATA.search(contents)
    if metadata_match:
        metadata = json.loads(metadata_match[1])
        if int(metadata["distance"]) != distance or not math.isclose(float(metadata["p"]), p):
            raise ValueError(f"{path}: filename and run metadata disagree")
        rounds = int(metadata["rounds"])
    else:
        cycles_match = CYCLES_LINE.search(contents)
        rounds = int(cycles_match[1]) if cycles_match else 1
    if rounds < 1:
        raise ValueError(f"{path}: rounds must be positive")
    return rounds


def read_points(folder, min_errors=10):
    points = []
    for path in sorted(folder.iterdir()):
        if not path.is_file():
            continue
        match = FILENAME.search(path.name)
        if match:
            contents = path.read_text()
            ler_match = LER_LINE.search(contents)
            if not ler_match:
                continue
            if min_errors:
                counts = PROGRESS_ERRORS.findall(contents)
                if not counts:
                    raise ValueError(f"{path}: no observed error count; use --min-errors 0 to fit without filtering")
                if int(counts[-1]) < min_errors:
                    continue
            distance, p = int(match["d"]), float(match["p"])
            rounds = rounds_from_output(contents, distance, p, path)
            points.append((distance, p, float(ler_match[1]) / rounds, path))
        elif path.suffix.lower() == ".csv":
            with path.open(newline="") as stream:
                reader = csv.DictReader(stream)
                if min_errors and (reader.fieldnames is None or "errors" not in reader.fieldnames):
                    raise ValueError(f"{path}: CSV needs an errors column for --min-errors {min_errors}")
                for row in reader:
                    if min_errors and int(row["errors"]) < min_errors:
                        continue
                    ler_key = "ler" if "ler" in row else "logical_error_rate"
                    points.append((int(row["d"]), float(row["p"]), float(row[ler_key]), path))
    valid = [point for point in points if point[0] > 0 and point[1] > 0 and point[2] > 0]
    if len(valid) < 2:
        raise ValueError(f"Need at least two positive data points with >= {min_errors} errors in {folder}")
    return valid


def fit(points):
    # Taking logs makes this linear in log(A) and log(B):
    # ln(LER) = ln(A) + k ln(p) - k ln(B), k=(d+1)/2.
    # Nonlinear least squares in log-rate space is solved by grid-free
    # Gauss-Newton iterations over ln(A), ln(B).
    x = [(d + 1) // 2 for d, _, _, _ in points]
    lp = [math.log(p) for _, p, _, _ in points]
    ly = [math.log(y) for _, _, y, _ in points]
    b = math.log(max(p for _, p, _, _ in points))
    a = sum(ly) / len(ly)
    for _ in range(100):
        residual = [y - (a + k * (p - b)) for k, p, y in zip(x, lp, ly)]
        j00 = len(x)
        j01 = -sum(x)
        j11 = sum(k * k for k in x)
        g0 = sum(residual)
        g1 = -sum(k * r for k, r in zip(x, residual))
        determinant = j00 * j11 - j01 * j01
        if determinant == 0:
            raise ValueError("Fit requires at least two distinct distances")
        da = (g0 * j11 - g1 * j01) / determinant
        db = (j00 * g1 - j01 * g0) / determinant
        a += da
        b += db
        if max(abs(da), abs(db)) < 1e-12:
            break
    return math.exp(a), math.exp(b)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, help="Folder containing .out or CSV data")
    parser.add_argument("--min-errors", type=nonnegative_int, default=10,
                        help="minimum observed logical errors per point (default: 10; 0 disables)")
    args = parser.parse_args()
    points = read_points(args.folder, args.min_errors)
    amplitude, threshold = fit(points)
    print(f"Points: {len(points)}")
    print(f"A = {amplitude:.10g}")
    print(f"B = {threshold:.10g}")
    print(f"LER = {amplitude:.10g} * (p/{threshold:.10g})**floor((d+1)/2)")

    for d in [15, 21, 25, 27, 29, 31, 33, 35]:
        for p in [0.001]:
            ler = amplitude * (p/threshold) ** ((d+1)//2)
            print(f'[d = {d}] projected ler @ p = {p} : {ler}')

if __name__ == "__main__":
    main()
