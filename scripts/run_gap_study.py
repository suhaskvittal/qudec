#!/usr/bin/env python3
"""Contributor: OpenAI Codex. Calibrate a PyMatching gap threshold and evaluate it."""

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
BIN = re.compile(rf"^({NUMBER})\s*<=\s*X\s*<\s*({NUMBER}):\s*(\d+)\s*$")
EXTREME = re.compile(rf"^(<|>=)({NUMBER}):\s*(\d+)\s*$")
HEADER = re.compile(r"^(SIGNED_GAP|UNSIGNED_GAP|UNSIGNED_GAP_ERRORS)\s*=+")


def positive_integer(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def parse_histograms(path):
    """Respect actual bin bounds, omitted zero bins, and under/overflow counts."""
    histograms, section = {}, None
    for line in path.read_text().splitlines():
        heading = HEADER.match(line)
        if heading:
            name = heading.group(1)
            section = None if name == "SIGNED_GAP" else name
            if section:
                histograms[section] = {"bins": {}, "extremes": [], "count": None}
            continue
        if section is None:
            continue
        match = BIN.match(line)
        if match:
            lo, hi, count = float(match[1]), float(match[2]), int(match[3])
            if lo < 0 or hi <= lo:
                raise ValueError("Invalid unsigned-gap bin")
            histograms[section]["bins"][(lo, hi)] = count
            continue
        match = EXTREME.match(line)
        if match:
            histograms[section]["extremes"].append((match[1], float(match[2]), int(match[3])))
        match = re.search(r"\bcount\s*=\s*(\d+)", line)
        if match:
            histograms[section]["count"] = int(match[1])
    for name in ("UNSIGNED_GAP", "UNSIGNED_GAP_ERRORS"):
        if name not in histograms:
            raise ValueError(f"Missing {name} in {path}; calibration may not have completed")
        histogram = histograms[name]
        count = sum(histogram["bins"].values()) + sum(item[2] for item in histogram["extremes"])
        if histogram["count"] is not None and count != histogram["count"]:
            raise ValueError(f"Histogram count mismatch for {name}")
        histogram["count"] = count
    return histograms["UNSIGNED_GAP"], histograms["UNSIGNED_GAP_ERRORS"]


def fit_threshold(total, errors, target, min_errors):
    points = []
    for bounds, error_count in sorted(errors["bins"].items()):
        shots = total["bins"].get(bounds, 0)
        if error_count > shots:
            raise ValueError("Gap-bin errors exceed gap-bin shots")
        if shots and error_count >= min_errors:
            points.append({"gap_low_db": bounds[0], "gap_high_db": bounds[1],
                           "gap_center_db": sum(bounds) / 2, "shots": shots,
                           "errors": error_count, "conditional_error_rate": error_count / shots})
    if len(points) < 2:
        raise ValueError(f"Only {len(points)} usable gap bins; increase calibration shots/error budget "
                         "or reduce --min-errors. A threshold cannot be inferred from zero errors.")
    x = np.array([p["gap_center_db"] for p in points])
    y = np.log10([p["conditional_error_rate"] for p in points])
    weight = np.array([p["errors"] for p in points], dtype=float)
    design = np.column_stack([x, np.ones(len(x))])
    coefficients, _, rank, _ = np.linalg.lstsq(design * np.sqrt(weight)[:, None],
                                              y * np.sqrt(weight), rcond=None)
    slope, intercept = map(float, coefficients)
    if rank != 2 or not math.isfinite(slope) or slope >= 0:
        raise ValueError("Calibration did not establish a decreasing error rate versus gap; collect more data")
    threshold = max(0., (math.log10(target) - intercept) / slope)
    if not math.isfinite(threshold):
        raise ValueError("Fitted threshold is not finite")
    residual = float(np.sum(weight * (y - design @ coefficients) ** 2))
    variation = float(np.sum(weight * (y - np.average(y, weights=weight)) ** 2))
    # Bracket the forwarding fraction rather than pretending a partially cut
    # histogram bin, or a lumped overflow bin, is known exactly.
    definitely_below, possibly_below = 0, 0
    for (lo, hi), count in total["bins"].items():
        if hi <= threshold:
            definitely_below += count
        if lo < threshold:
            possibly_below += count
    for kind, edge, count in total["extremes"]:
        if kind == "<":
            if threshold >= edge: definitely_below += count
            if threshold > 0: possibly_below += count
        elif threshold > edge:
            possibly_below += count
    count = total["count"]
    result = {"threshold_db": threshold, "source": "fit", "target_conditional_error_rate": target,
              "slope_log10_per_db": slope, "intercept_log10": intercept,
              "weighted_r_squared": 1 - residual / variation if variation else None,
              "usable_bins": len(points), "fit_gap_range_db": [float(x.min()), float(x.max())],
              "extrapolated": bool(threshold < x.min() or threshold > x.max()),
              "calibration_shots": count, "calibration_errors": errors["count"],
              "forward_fraction_bounds": [definitely_below / count, possibly_below / count]}
    return result, points


def mpi_enabled(directory):
    cache = directory / "CMakeCache.txt"
    if cache.is_file():
        for line in cache.read_text().splitlines():
            if line.startswith("ENABLE_MPI:BOOL="):
                return line.split("=", 1)[1] == "ON"
    return None


def find_binaries(args):
    directories = ([args.build_dir.expanduser().resolve()] if args.build_dir else
                   [REPOSITORY / p for p in ("build_mpi", "build", "build_nompi",
                    ".agentic_artifacts/libra_integration/build", ".agentic_artifacts/libra_integration/serial")])
    usable = [p for p in directories if (p / "qudec_eval_gap").is_file()
              and (args.gap_threshold_db is not None or (p / "qudec").is_file())
              and not (args.cores > 1 and mpi_enabled(p) is False)]
    usable.sort(key=lambda p: (p / "qudec_eval_gap").stat().st_mtime if args.gap_threshold_db is not None
                else min((p / "qudec").stat().st_mtime, (p / "qudec_eval_gap").stat().st_mtime), reverse=True)
    selected = usable[0] if usable else None
    qudec = args.qudec.expanduser().resolve() if args.qudec else selected / "qudec" if selected else None
    evaluator = args.evaluator.expanduser().resolve() if args.evaluator else selected / "qudec_eval_gap" if selected else None
    for binary in ((evaluator,) if args.gap_threshold_db is not None else (qudec, evaluator)):
        if binary is None or not binary.is_file() or not os.access(binary, os.X_OK):
            raise ValueError("Build qudec and qudec_eval_gap, or provide --build-dir / --qudec / --evaluator")
        if args.cores > 1 and mpi_enabled(binary.parent) is False:
            raise ValueError(f"Not an MPI build: {binary}; use --cores 1 or an MPI build")
    return qudec, evaluator


def launch_prefix(args):
    if args.cores == 1:
        return []
    launcher = shutil.which(args.mpi_launcher)
    if launcher is None:
        raise ValueError(f"Cannot find Open MPI launcher: {args.mpi_launcher}")
    return [launcher, "-np", str(args.cores), "--bind-to", "core", "--map-by", "core"]


def run_logged(command, path):
    print(shlex.join(command), flush=True)
    with path.open("w") as log:
        with subprocess.Popen(command, cwd=REPOSITORY, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, bufsize=1) as process:
            for line in process.stdout:
                log.write(line)
                if (line.startswith("[") or "mean_latency_us=" in line or line.startswith("Stage ")
                        or line.startswith("two_level_logical_error_rate=")):
                    print(line, end="", flush=True)
            status = process.wait()
    if status:
        tail = "\n".join(path.read_text().splitlines()[-20:])
        raise RuntimeError(f"Command failed ({status}); see {path}\n{tail}")


def evaluation_command(args, evaluator, output, threshold):
    return launch_prefix(args) + [str(evaluator), str(args.distance), "-r", str(args.rounds or args.distance),
            "-p", str(args.p), "--method", args.method, "--seed", str(args.seed),
            "--gap-threshold-db", str(threshold), "--mc-max-samples", str(args.shots),
            "--mc-stop-at-errors", str(args.shots + 1),
            "--libra-ensemble-size", str(args.libra_ensemble_size),
            "--libra-gap-threshold-db", str(args.libra_gap_threshold_db),
            "--libra-seed", str(args.libra_seed),
            "--rare-samples-per-level", str(args.rare_samples_per_level),
            "--rare-max-errors-per-level", str(args.rare_max_errors_per_level),
            "--rare-max-level", str(args.rare_max_level),
            "--rare-relative-tail-tolerance", str(args.rare_relative_tail_tolerance),
            "--rare-tail-tolerance", str(args.rare_tail_tolerance),
            "--output-dir", str(output), "--print-progress"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance", type=positive_integer, required=True)
    parser.add_argument("--shots", type=positive_integer, default=1000000, help="evaluation MC budget (global)")
    parser.add_argument("--calibration-shots", type=positive_integer, default=100000000)
    parser.add_argument("--calibration-max-errors", type=positive_integer, default=10000)
    threshold_source = parser.add_mutually_exclusive_group()
    threshold_source.add_argument("--calibration-log", type=Path, help="reuse a completed PyMatching --gap log instead of sampling calibration")
    threshold_source.add_argument("--gap-threshold-db", type=float, help="use this threshold directly; skip calibration and fitting")
    parser.add_argument("--target-error-rate", type=float, default=6.9e-11, help="target conditional P(error | gap); default matches the former fit script")
    parser.add_argument("--min-errors", type=positive_integer, default=5, help="minimum observed errors per fitted gap bin")
    parser.add_argument("--cores", type=positive_integer, default=64)
    parser.add_argument("--p", type=float, default=.001)
    parser.add_argument("--rounds", type=positive_integer)
    parser.add_argument("--seed", type=int, default=0, help="evaluation seed")
    parser.add_argument("--calibration-seed", type=int, help="default: evaluation seed + cores, giving disjoint rank streams")
    parser.add_argument("--method", choices=("monte_carlo", "rare_event"), default="monte_carlo", help="evaluation sampler; calibration is always Monte Carlo")
    parser.add_argument("--qudec", type=Path)
    parser.add_argument("--evaluator", type=Path)
    parser.add_argument("--build-dir", type=Path, help="configured CMake directory; build both targets before running")
    parser.add_argument("--mpi-launcher", default="mpirun")
    parser.add_argument("--output-root", type=Path, default=REPOSITORY / "out")
    parser.add_argument("--libra-ensemble-size", type=positive_integer, default=60)
    parser.add_argument("--libra-gap-threshold-db", type=float, default=20.)
    parser.add_argument("--libra-seed", type=int, default=0x4C49425241)
    parser.add_argument("--rare-samples-per-level", type=positive_integer, default=10000)
    parser.add_argument("--rare-max-errors-per-level", type=positive_integer, default=25)
    parser.add_argument("--rare-max-level", type=positive_integer, default=128)
    parser.add_argument("--rare-relative-tail-tolerance", type=float, default=.1)
    parser.add_argument("--rare-tail-tolerance", type=float, default=1e-12)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.calibration_seed is None: args.calibration_seed = args.seed + args.cores
    if args.distance < 2 or max(args.distance, args.rounds or args.distance) > 2**32 - 1:
        parser.error("distance must be at least 2; distance and rounds must fit uint32")
    if max(args.shots, args.calibration_shots, args.calibration_max_errors) > 2**63 - 2:
        parser.error("shot/error budgets must fit int64 with room for shots + 1")
    if min(args.seed, args.calibration_seed, args.libra_seed) < 0 or max(args.seed, args.calibration_seed) > 2**63 - args.cores:
        parser.error("sampling seeds must be nonnegative int64 values with room for rank offsets")
    if args.libra_seed > 2**63 - 1:
        parser.error("Libra seed must fit nonnegative int64")
    if not math.isfinite(args.p) or not 0 < args.p < .5:
        parser.error("physical error rate must be between 0 and 0.5")
    if not math.isfinite(args.target_error_rate) or not 0 < args.target_error_rate < 1:
        parser.error("target error rate must be between 0 and 1")
    if args.gap_threshold_db is not None and (not math.isfinite(args.gap_threshold_db) or args.gap_threshold_db < 0):
        parser.error("gap threshold must be finite and nonnegative")
    if not math.isfinite(args.libra_gap_threshold_db) or args.libra_gap_threshold_db < 0:
        parser.error("Libra gap threshold must be finite and nonnegative")
    if not math.isfinite(args.rare_relative_tail_tolerance) or not 0 <= args.rare_relative_tail_tolerance <= 1:
        parser.error("relative tail tolerance must be in [0, 1]")
    if not math.isfinite(args.rare_tail_tolerance) or not 0 <= args.rare_tail_tolerance < 1:
        parser.error("tail tolerance must be in [0, 1)")
    if args.calibration_log and not args.calibration_log.is_file():
        parser.error("calibration log does not exist")
    if args.build_dir:
        build = args.build_dir.expanduser().resolve()
        if not (build / "CMakeCache.txt").is_file():
            parser.error("--build-dir must already be configured with CMake")
        targets = ["qudec_eval_gap"] if args.gap_threshold_db is not None else ["qudec", "qudec_eval_gap"]
        build_command = ["cmake", "--build", str(build), "--target"] + targets + ["-j", str(args.cores)]
        if args.dry_run:
            print(shlex.join(build_command))
        else:
            subprocess.run(build_command, cwd=REPOSITORY, check=True)
    try:
        qudec, evaluator = find_binaries(args)
        prefix = launch_prefix(args)
    except ValueError as error:
        parser.error(str(error))
    directory = args.output_root.expanduser().resolve() / (
        f"gap_study_d{args.distance}_seed{args.seed}_" + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f"))
    calibration = prefix + [str(qudec), "pymatching", str(args.distance), "-r", str(args.rounds or args.distance),
                   "-p", str(args.p), "--gap", "--include-opposite-basis-detectors", "--method", "monte_carlo",
                   "--seed", str(args.calibration_seed), "--mc-max-samples", str(args.calibration_shots),
                   "--mc-stop-at-errors", str(args.calibration_max_errors), "--print-progress"]
    print(f"Run directory: {directory}", flush=True)
    if args.dry_run:
        if args.gap_threshold_db is None:
            print("Calibration log: " + str(args.calibration_log) if args.calibration_log else shlex.join(calibration))
        threshold = args.gap_threshold_db if args.gap_threshold_db is not None else "COMPUTED_THRESHOLD_DB"
        print(shlex.join(evaluation_command(args, evaluator, directory / "evaluation", threshold)))
        return 0
    directory.mkdir(parents=True, exist_ok=False)
    manifest = {"contributor": "OpenAI Codex", "started_utc": datetime.now(timezone.utc).isoformat(),
                "configuration": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                "calibration_command": calibration if args.gap_threshold_db is None and not args.calibration_log else None,
                "status": "running", "working_directory": str(REPOSITORY)}
    start = time.monotonic()
    status = 0
    try:
        calibration_lines = []
        if args.gap_threshold_db is not None:
            fit = {"threshold_db": args.gap_threshold_db, "source": "provided"}
            print(f"Using supplied gap threshold: {args.gap_threshold_db:.6f} dB (calibration skipped)", flush=True)
        else:
            calibration_log = directory / "calibration.log"
            if args.calibration_log:
                shutil.copyfile(args.calibration_log, calibration_log)
            else:
                run_logged(calibration, calibration_log)
            total, errors = parse_histograms(calibration_log)
            fit, points = fit_threshold(total, errors, args.target_error_rate, args.min_errors)
            with (directory / "fit_points.csv").open("w", newline="") as target:
                writer = csv.DictWriter(target, fieldnames=list(points[0]))
                writer.writeheader()
                writer.writerows(points)
            print(f"Computed threshold: {fit['threshold_db']:.6f} dB; "
                  f"target P(error | gap)={args.target_error_rate:.6g}", flush=True)
            if fit["extrapolated"]:
                print("Threshold extrapolates beyond the fitted gap range; target rate is a model prediction.", flush=True)
            calibration_lines = [f"Target conditional error rate: {args.target_error_rate:.6g}",
                 f"Calibration: {fit['calibration_shots']:,} shots, {fit['calibration_errors']:,} errors, {fit['usable_bins']} fitted bins",
                 f"Fit: log10(P(error | gap)) = {fit['slope_log10_per_db']:.6g} * gap + {fit['intercept_log10']:.6g}",
                 f"Extrapolated threshold: {fit['extrapolated']}"]
        (directory / "threshold.json").write_text(json.dumps(fit, indent=2, allow_nan=False) + "\n")
        threshold = fit["threshold_db"]
        command = evaluation_command(args, evaluator, directory / "evaluation", format(threshold, ".17g"))
        manifest["evaluation_command"] = command
        calibration_provenance = "Calibration skipped: threshold supplied directly"
        if args.gap_threshold_db is None:
            calibration_provenance = shlex.join(calibration) if not args.calibration_log else f"Calibration input: {args.calibration_log}"
        (directory / "commands.txt").write_text(
            "Contributor: OpenAI Codex.\n\n" + calibration_provenance
            + "\n\n" + shlex.join(command) + "\n")
        run_logged(command, directory / "evaluation.log")
        summary = (directory / "evaluation" / "summary.txt").read_text()
        lines = [f"PyMatching gap threshold: {threshold:.6f} dB (source: {fit['source']})"] + calibration_lines + ["", summary]
        (directory / "study_summary.txt").write_text("\n".join(lines))
        manifest["status"] = "complete"
        print(f"Completed. Threshold, fit, decoder results, latencies, and low-gap syndromes: {directory}")
    except (OSError, ValueError, RuntimeError, KeyboardInterrupt) as error:
        manifest["status"] = "failed"
        manifest["error"] = str(error) or type(error).__name__
        print(f"Study failed: {error}\nFiles preserved in {directory}", file=sys.stderr)
        status = 1
    finally:
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["elapsed_seconds"] = time.monotonic() - start
        (directory / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return status


if __name__ == "__main__":
    sys.exit(main())
