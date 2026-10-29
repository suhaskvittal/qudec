#!/usr/bin/env python3
"""Contributor: OpenAI Codex. Collect Monte Carlo cluster-gap statistics."""

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

REPOSITORY = Path(__file__).resolve().parents[1]
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
BIN = re.compile(rf"^({NUMBER})\s*<=\s*X\s*<\s*({NUMBER}):\s*(\d+)\s*$")
EXTREME = re.compile(rf"^(<|>=)({NUMBER}):\s*(\d+)\s*$")
HEADER = re.compile(r"^(CLUSTER_GAP_DB|CLUSTER_GAP_DB_ERRORS)\s*=+")


def positive(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def parse_histograms(path):
    histograms, current = {}, None
    for line in path.read_text().splitlines():
        header = HEADER.match(line)
        if header:
            current = {"bins": {}, "extremes": [], "count": None, "mean_db": None,
                       "std_db": None, "min_db": None, "max_db": None}
            histograms[header[1]] = current
            continue
        if current is None:
            continue
        match = BIN.match(line)
        if match:
            current["bins"][(float(match[1]), float(match[2]))] = int(match[3])
            continue
        match = EXTREME.match(line)
        if match:
            current["extremes"].append((match[1], float(match[2]), int(match[3])))
            continue
        match = re.search(r"mean\s*=\s*([^,]+).*std\s*=\s*([^,]+).*min\s*=\s*([^,]+).*max\s*=\s*([^,]+).*count\s*=\s*(\d+)", line)
        if match:
            current.update(mean_db=float(match[1]), std_db=float(match[2]),
                           min_db=float(match[3]), max_db=float(match[4]), count=int(match[5]))
            continue
        match = re.fullmatch(r"count\s*=\s*0", line.strip())
        if match:
            current["count"] = 0
    for name in ("CLUSTER_GAP_DB", "CLUSTER_GAP_DB_ERRORS"):
        if name not in histograms:
            raise ValueError(f"Missing {name} histogram in {path}; rerun with a qudec build that records cluster-gap errors")
        hist = histograms[name]
        observed = sum(hist["bins"].values()) + sum(x[2] for x in hist["extremes"])
        if hist["count"] != observed:
            raise ValueError(f"Cluster histogram count mismatch for {name}: header {hist['count']}, bins {observed}")
    return histograms


def mpi_capable(binary):
    cache = binary.parent / "CMakeCache.txt"
    if cache.is_file():
        return any(line.upper() == "ENABLE_MPI:BOOL=ON"
                   for line in cache.read_text().splitlines())
    try:
        deps = subprocess.run(["ldd", str(binary)], text=True, capture_output=True, check=False).stdout
        return "libmpi" in deps.lower()
    except OSError:
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decoder", choices=("pymatching", "uf"), default="pymatching")
    parser.add_argument("--distance", type=positive, required=True)
    parser.add_argument("--shots", type=positive, default=1_000_000)
    parser.add_argument("--max-errors", type=positive, default=10_000)
    parser.add_argument("--cores", type=positive, default=64)
    parser.add_argument("--p", type=float, default=.001)
    parser.add_argument("--rounds", type=positive)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--uf-weight-scale", type=positive, default=1024)
    parser.add_argument("--qudec", type=Path)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--mpi-launcher", default="mpirun")
    parser.add_argument("--output-root", type=Path, default=REPOSITORY / "out")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.distance < 2 or max(args.distance, args.rounds or args.distance) > 2**32 - 1:
        parser.error("distance must be at least 2; distance and rounds must fit uint32")
    if args.seed < 0 or args.seed > 2**63 - args.cores:
        parser.error("seed must be nonnegative int64 with room for rank offsets")
    if not 0 < args.p < .5:
        parser.error("physical error rate must be between 0 and 0.5")
    binary = args.qudec.expanduser().resolve() if args.qudec else None
    if args.build_dir:
        build = args.build_dir.expanduser().resolve()
        if not (build / "CMakeCache.txt").is_file():
            parser.error("--build-dir must already be configured with CMake")
        build_cmd = ["cmake", "--build", str(build), "--target", "qudec", "-j", str(args.cores)]
        if args.dry_run:
            print(shlex.join(build_cmd))
        else:
            subprocess.run(build_cmd, cwd=REPOSITORY, check=True)
        if binary is None:
            binary = build / "qudec"
    if binary is None:
        candidates = [REPOSITORY / p / "qudec" for p in ("build_mpi", "build", "build_nompi",
                      ".agentic_artifacts/libra_integration/build", ".agentic_artifacts/libra_integration/serial")]
        candidates = [p for p in candidates if p.is_file() and (args.cores == 1 or mpi_capable(p))]
        if candidates:
            binary = max(candidates, key=lambda p: p.stat().st_mtime)
    if (binary is None or not binary.is_file() or not os.access(binary, os.X_OK)) and args.dry_run and args.build_dir:
        binary = args.build_dir.expanduser().resolve() / "qudec"
    if binary is None or (not args.dry_run and (not binary.is_file() or not os.access(binary, os.X_OK))):
        parser.error("Build qudec or provide --qudec / --build-dir")
    if args.cores > 1 and not mpi_capable(binary):
        parser.error(f"Not an MPI build: {binary}; use --cores 1 or an MPI build")
    prefix = []
    if args.cores > 1:
        launcher = shutil.which(args.mpi_launcher)
        if not launcher:
            parser.error(f"Cannot find MPI launcher: {args.mpi_launcher}")
        prefix = [launcher, "-np", str(args.cores), "--bind-to", "core", "--map-by", "core"]
    directory = args.output_root.expanduser().resolve() / (
        f"cluster_gap_{args.decoder}_d{args.distance}_" + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f"))
    command = prefix + [str(binary), args.decoder, "-d", str(args.distance), "-r", str(args.rounds or args.distance),
        "-p", str(args.p), "--cluster-gap", "--cluster-gap-histogram", "--include-opposite-basis-detectors",
        "--method", "monte_carlo", "--seed", str(args.seed), "--mc-max-samples", str(args.shots),
        "--mc-stop-at-errors", str(args.max_errors), "--print-progress"]
    if args.decoder == "uf":
        command += ["--uf-weight-scale", str(args.uf_weight_scale)]
    print(f"Run directory: {directory}")
    print(shlex.join(command))
    if args.dry_run:
        return 0
    directory.mkdir(parents=True, exist_ok=False)
    (directory / "commands.txt").write_text("Contributor: OpenAI Codex.\n\n" + shlex.join(command) + "\n")
    log = directory / "sampling.log"
    try:
        with log.open("w") as output:
            subprocess.run(command, cwd=REPOSITORY, stdout=output, stderr=subprocess.STDOUT, check=True)
        histograms = parse_histograms(log)
        hist, errors = histograms["CLUSTER_GAP_DB"], histograms["CLUSTER_GAP_DB_ERRORS"]
        logical_error_rate = errors["count"] / hist["count"] if hist["count"] else None
        summary = {"contributor": "OpenAI Codex", "decoder": args.decoder, "sampler": "monte_carlo",
                   "distance": args.distance, "rounds": args.rounds or args.distance, "p": args.p,
                   "seed": args.seed, "shots_requested": args.shots, "max_errors": args.max_errors,
                   "shots_observed": hist["count"], "logical_errors_observed": errors["count"],
                   "logical_error_rate": logical_error_rate,
                   "cluster_gap_db": {k: v for k, v in hist.items() if k not in ("bins", "extremes")}}
        (directory / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        with (directory / "cluster_gap_bins.csv").open("w", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(["histogram", "bin_low_db", "bin_high_db", "extreme", "edge_db", "count"])
            for name, value in histograms.items():
                for (lo, hi), count in sorted(value["bins"].items()):
                    writer.writerow([name, lo, hi, "", "", count])
                for direction, edge, count in value["extremes"]:
                    writer.writerow([name, "", "", direction, edge, count])
        print(json.dumps(summary, indent=2))
        print(f"Completed cluster-gap collection: {directory}")
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"Cluster-gap study failed: {error}; files preserved in {directory}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
