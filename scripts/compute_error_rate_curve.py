#!/usr/bin/env python3
"""Fit LER = A * (p/B)**floor((d+1)/2) to simulation outputs in a folder.

The .out format used by this repository stores distance and physical error rate
in filenames (d5_p1.00e-01.out), and the logical error rate in the file body.
CSV files with columns d, p, and ler (or logical_error_rate) are also accepted.
"""

import argparse
import csv
import math
import re
from pathlib import Path


FILENAME = re.compile(r"d(?P<d>\d+)_p(?P<p>[0-9.eE+-]+)\.out$")
LER_LINE = re.compile(r"^Logical error rate:\s*([0-9.eE+-]+)\s*$", re.MULTILINE)
CYCLES_LINE = re.compile(r"\brounds=(\d+)\b")


def read_points(folder):
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
            cycles_match = CYCLES_LINE.search(contents)
            cycles = int(cycles_match[1]) if cycles_match else 1
            points.append((int(match["d"]), float(match["p"]), float(ler_match[1]) / cycles, path))
        elif path.suffix.lower() == ".csv":
            with path.open(newline="") as stream:
                reader = csv.DictReader(stream)
                for row in reader:
                    ler_key = "ler" if "ler" in row else "logical_error_rate"
                    points.append((int(row["d"]), float(row["p"]), float(row[ler_key]), path))
    valid = [point for point in points if point[0] > 0 and point[1] > 0 and point[2] > 0]
    if len(valid) < 2:
        raise ValueError(f"Need at least two positive data points in {folder}")
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
    args = parser.parse_args()
    points = read_points(args.folder)
    amplitude, threshold = fit(points)
    print(f"Points: {len(points)}")
    print(f"A = {amplitude:.10g}")
    print(f"B = {threshold:.10g}")
    print(f"LER = {amplitude:.10g} * (p/{threshold:.10g})**floor((d+1)/2)")

    for d in [15, 21, 25]:
        for p in [0.01, 0.001]:
            ler = amplitude * (p/threshold) ** ((d+1)//2)
            print(f'[d = {d}] projected ler @ p = {p} : {ler}')

if __name__ == "__main__":
    main()
