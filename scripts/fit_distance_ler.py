#!/usr/bin/env python3
"""Fit log10(logical error rate) as a linear function of code distance.

Contributor: OpenAI Codex.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path


RATE_PATTERN = re.compile(
    r"^Logical error rate:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$",
    re.MULTILINE,
)
DISTANCE_PATTERN = re.compile(r"^d(\d+)$")


@dataclass(frozen=True)
class Measurement:
    distance: int
    rate: float


def read_measurements(folder: Path) -> tuple[list[Measurement], list[int]]:
    if not folder.is_dir():
        raise ValueError(f"Data folder does not exist: {folder}")

    measurements: list[Measurement] = []
    zero_rate_distances: list[int] = []
    for directory in folder.iterdir():
        match = DISTANCE_PATTERN.fullmatch(directory.name)
        if match is None or not directory.is_dir():
            continue
        distance = int(match.group(1))
        output = directory / "output.log"
        if not output.is_file():
            raise ValueError(f"Missing output.log for distance {distance}: {output}")
        rates = RATE_PATTERN.findall(output.read_text())
        if len(rates) != 1:
            raise ValueError(f"Expected exactly one logical error rate in {output}; found {len(rates)}")
        rate = float(rates[0])
        if not math.isfinite(rate) or not 0 <= rate <= 1:
            raise ValueError(f"Invalid logical error rate in {output}: {rates[0]}")
        if rate == 0:
            zero_rate_distances.append(distance)
        else:
            measurements.append(Measurement(distance, rate))

    measurements.sort(key=lambda item: item.distance)
    zero_rate_distances.sort()
    if not measurements and not zero_rate_distances:
        raise ValueError(f"No dN/output.log results found in {folder}")
    return measurements, zero_rate_distances


@dataclass(frozen=True)
class Fit:
    intercept: float
    slope: float
    r_squared: float

    def predict(self, distance: int) -> float:
        return 10 ** (self.intercept + self.slope * distance)


def fit_log_rate(measurements: list[Measurement]) -> Fit:
    if len(measurements) < 2:
        raise ValueError("At least two positive logical error rates are needed for a fit")
    xs = [item.distance for item in measurements]
    ys = [math.log10(item.rate) for item in measurements]
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    spread_x = sum((x - mean_x) ** 2 for x in xs)
    if spread_x == 0:
        raise ValueError("The fit needs results from at least two distinct distances")
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / spread_x
    intercept = mean_y - slope * mean_x
    residual = sum((y - intercept - slope * x) ** 2 for x, y in zip(xs, ys))
    total = sum((y - mean_y) ** 2 for y in ys)
    r_squared = 1 - residual / total if total else (1.0 if residual == 0 else float("nan"))
    return Fit(intercept, slope, r_squared)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, help="Folder containing dN/output.log results")
    parser.add_argument("--min-distance", type=int, help="Smallest distance included in the fit")
    parser.add_argument("--max-distance", type=int, help="Largest distance included in the fit")
    parser.add_argument("--predict", type=int, nargs="*", default=[17, 21, 25, 29], metavar="D",
                        help="Distances to extrapolate (default: 17 21 25 29)")
    args = parser.parse_args()

    try:
        measurements, zeros = read_measurements(args.folder)
        selected = [item for item in measurements
                    if (args.min_distance is None or item.distance >= args.min_distance)
                    and (args.max_distance is None or item.distance <= args.max_distance)]
        fit = fit_log_rate(selected)
    except ValueError as error:
        parser.error(str(error))

    print(f"Data folder: {args.folder}")
    print(f"Fit: log10(LER) = {fit.intercept:.8g} {fit.slope:+.8g} * d")
    print(f"Equivalent: LER = {10 ** fit.intercept:.8g} * 10^({fit.slope:.8g} * d)")
    print(f"R^2 in log10 space: {fit.r_squared:.6g}")
    print(f"LER multiplier for d + 2: {10 ** (2 * fit.slope):.8g}")
    print("distance  measured LER   fitted LER    fitted/measured")
    for item in selected:
        predicted = fit.predict(item.distance)
        print(f"{item.distance:8d}  {item.rate:12.6g}  {predicted:12.6g}  {predicted / item.rate:15.6g}")
    if zeros:
        print("Zero-rate distances omitted from fit: " + ", ".join(map(str, zeros)))
    if args.predict:
        print("predicted distance  predicted LER")
        for distance in args.predict:
            print(f"{distance:18d}  {fit.predict(distance):.8g}")


if __name__ == "__main__":
    main()
