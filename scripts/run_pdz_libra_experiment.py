#!/usr/bin/env python3
"""Contributor: OpenAI Codex. Run and analyze a recorded PDZ + Libra experiment."""

import argparse
import csv
from datetime import datetime, timezone
import importlib.util
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time


REPOSITORY = Path(__file__).resolve().parents[1]
ANALYZER = REPOSITORY / "scripts" / "analyze_pdz_libra_records.py"


def positive_integer(value):
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return result


def mpi_setting(executable):
    """Check the adjacent CMake cache when available; None means unknown."""
    cache = executable.parent / "CMakeCache.txt"
    if cache.is_file():
        for line in cache.read_text().splitlines():
            if line.startswith("ENABLE_MPI:BOOL="):
                return line.split("=", 1)[1] == "ON"
    return None


def find_executable(specified, cores):
    if specified is not None:
        candidates = [specified.expanduser().resolve()]
    else:
        candidates = [REPOSITORY / path for path in (
            "build_mpi/qudec", "build/qudec", ".agentic_artifacts/libra_integration/build/qudec",
            "build_nompi/qudec", ".agentic_artifacts/libra_integration/serial/qudec",
        )]
    # Prefer the newest compatible build rather than accidentally selecting an
    # older executable from another build directory.
    candidates.sort(key=lambda path: path.stat().st_mtime if path.is_file() else -1, reverse=True)
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            if cores > 1 and mpi_setting(candidate) is False:
                if specified is not None:
                    raise ValueError("The selected qudec was built without MPI; use --cores 1 or an MPI build")
                continue
            return candidate.resolve()
    raise ValueError("Cannot find a suitable built qudec; specify --qudec PATH (MPI build for multiple cores)")


def make_command(args, executable, directory):
    command = [str(executable), "pdz_libra", str(args.distance),
               "-r", str(args.rounds or args.distance), "-p", str(args.p),
               "-m", args.method, "--seed", str(args.seed),
               "--pdz-aff-keep-threshold", str(args.pdz_aff_keep_threshold),
               "--pdz-removal-threshold", str(args.pdz_removal_threshold),
               "--libra-ensemble-size", str(args.libra_ensemble_size),
               "--libra-gap-threshold-db", str(args.libra_gap_threshold_db),
               "--libra-seed", str(args.libra_seed), "--record-syndromes",
               "--record-syndromes-dir", str(directory), "--print-progress"]
    if args.method == "monte_carlo":
        command.extend(["--mc-max-samples", str(args.shots),
                        "--mc-stop-at-errors", str(args.shots + 1)])
    else:
        command.extend(["--rare-samples-per-level", str(args.rare_samples_per_level or args.shots),
                        "--rare-max-errors-per-level", str(args.rare_max_errors_per_level),
                        "--rare-max-level", str(args.rare_max_level),
                        "--rare-relative-tail-tolerance", str(args.rare_relative_tail_tolerance),
                        "--rare-tail-tolerance", str(args.rare_tail_tolerance)])
    if args.pdz_consider_compl_match:
        command.append("--pdz-consider-compl-match")
    if args.cores > 1:
        launcher = shutil.which(args.mpi_launcher)
        if launcher is None:
            raise ValueError(f"Cannot find MPI launcher: {args.mpi_launcher}")
        command = [launcher, "-np", str(args.cores), "--bind-to", "core", "--map-by", "core"] + command
    return command


def run_logged(command, log_path, progress=False):
    """Use argument lists, preserve full logs, and show sampling progress."""
    with log_path.open("w") as log:
        with subprocess.Popen(command, cwd=REPOSITORY, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, bufsize=1) as process:
            for line in process.stdout:
                log.write(line)
                if progress and (line.startswith("[") or line.startswith("Logical error rate:")):
                    print(line, end="", flush=True)
            status = process.wait()
    if status:
        tail = "\n".join(log_path.read_text().splitlines()[-20:])
        raise RuntimeError(f"Command exited with status {status}; see {log_path}\n{tail}")


def format_summary(directory, summary):
    groups = [[0, 0, 0, 0] for _ in range(4)]
    empty = 0
    with (directory / "analysis_by_removal_fraction.csv").open(newline="") as source:
        for row in csv.DictReader(source):
            bucket = int(row["fraction_bucket"])
            if bucket < 0:
                empty += int(row["shots"])
                continue
            index = 0 if bucket <= 2 else 1 if bucket <= 4 else 2 if bucket <= 6 else 3
            for i, key in enumerate(("shots", "ensemble", "helped", "harmed")):
                groups[index][i] += int(row[key])
    bypassed = summary.get("bypassed_shots", 0)
    lines = [f"Completed {summary['shots']:,} shots across {len(summary['rank_shots'])} ranks.",
             f"Libra bypassed: {bypassed:,} shots.",
             f"PyMatching failures: {summary['pymatching_failures']}; pipeline failures: {summary['pipeline_failures']}",
             "", "PDZ removed     Shots       Ensemble    Errors fixed    Errors introduced"]
    if summary.get("sampler") == "rare_event":
        lines.insert(0, "Rare-event tables describe unweighted conditional samples, not physical frequencies.")
    for label, (shots, ensemble, helped, harmed) in zip(("0-20%", ">20-40%", ">40-60%", ">60%"), groups):
        lines.append(f"{label:12} {shots:11,d} {ensemble:12,d} {helped:14,d} {harmed:20,d}")
    lines.extend([f"Empty input syndromes excluded from the table: {empty:,}",
                  "Fixed/introduced compare final pipeline and initial PyMatching on ensemble shots.",
                  "", f"Data, logs, reports, and CSV tables: {directory}"])
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance", type=positive_integer, required=True)
    parser.add_argument("--shots", type=positive_integer, required=True,
                        help="global Monte Carlo shots, or rare-event sample cap per level")
    parser.add_argument("--method", choices=("auto", "monte_carlo", "rare_event"), default="auto",
                        help="default: Monte Carlo below d=13, rare event at d>=13")
    parser.add_argument("--rare-samples-per-level", type=positive_integer,
                        help="override the rare-event per-level cap (otherwise --shots)")
    parser.add_argument("--rare-max-errors-per-level", type=positive_integer, default=25)
    parser.add_argument("--rare-max-level", type=positive_integer, default=128)
    parser.add_argument("--rare-relative-tail-tolerance", type=float, default=0.1)
    parser.add_argument("--rare-tail-tolerance", type=float, default=1e-12)
    parser.add_argument("--cores", type=positive_integer, default=64, help="MPI ranks (default: 64; 1 for serial)")
    parser.add_argument("--rounds", type=positive_integer, help="default: distance")
    parser.add_argument("--p", type=float, default=0.001, help="physical error rate (default: 0.001)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--qudec", type=Path, help="built qudec executable; otherwise discover an existing build")
    parser.add_argument("--mpi-launcher", default="mpirun", help="Open MPI launcher (default: mpirun)")
    parser.add_argument("--output-root", type=Path, default=REPOSITORY / "out")
    parser.add_argument("--pdz-aff-keep-threshold", type=int, default=2)
    parser.add_argument("--pdz-removal-threshold", type=float, default=.5,
                        help="fraction removed to bypass Libra (default: 0.5; 1 requires empty residual)")
    parser.add_argument("--pdz-consider-compl-match", action="store_true")
    parser.add_argument("--libra-ensemble-size", type=positive_integer, default=60)
    parser.add_argument("--libra-gap-threshold-db", type=float, default=20.)
    parser.add_argument("--libra-seed", type=int, default=0x4C49425241)
    parser.add_argument("--dry-run", action="store_true", help="print commands without creating files or running sampling")
    args = parser.parse_args()
    if args.method == "auto":
        args.method = "rare_event" if args.distance >= 13 else "monte_carlo"
    if args.distance < 2 or args.distance > 2**32 - 1 or (args.rounds or args.distance) > 2**32 - 1:
        parser.error("distance must be at least 2; distance and rounds must fit uint32")
    if args.shots > 2**63 - 2:
        parser.error("shot count must be at most 2**63 - 2")
    for name in ("rare_samples_per_level", "rare_max_errors_per_level", "rare_max_level"):
        if getattr(args, name) is not None and getattr(args, name) > 2**63 - 1:
            parser.error(f"{name} must fit int64")
    if not math.isfinite(args.rare_relative_tail_tolerance) or args.rare_relative_tail_tolerance < 0:
        parser.error("rare relative tail tolerance must be finite and nonnegative")
    if not math.isfinite(args.rare_tail_tolerance) or not 0 <= args.rare_tail_tolerance < 1:
        parser.error("rare tail tolerance must be finite and in [0, 1)")
    if not math.isfinite(args.p) or not 0 < args.p < .5:
        parser.error("physical error rate must be finite and between 0 and 0.5")
    if min(args.seed, args.libra_seed) < 0 or max(args.seed, args.libra_seed) > 2**63 - 1 - args.cores:
        parser.error("seeds must be nonnegative int64 values with room for rank offsets")
    if not 0 <= args.pdz_aff_keep_threshold <= 31:
        parser.error("PDZ threshold must be between 0 and 31")
    if not math.isfinite(args.pdz_removal_threshold) or not 0 <= args.pdz_removal_threshold <= 1:
        parser.error("PDZ removal threshold must be in [0, 1]")
    if not math.isfinite(args.libra_gap_threshold_db) or args.libra_gap_threshold_db < 0:
        parser.error("Libra gap threshold must be finite and nonnegative")
    if importlib.util.find_spec("numpy") is None:
        parser.error("analysis requires NumPy in this Python environment")
    try:
        executable = find_executable(args.qudec, args.cores)
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
        directory = args.output_root.expanduser().resolve() / (
            f"pdz_libra_d{args.distance}_{args.method}_{args.shots}_shots_seed{args.seed}_{timestamp}")
        command = make_command(args, executable, directory)
    except ValueError as error:
        parser.error(str(error))
    analysis_command = [sys.executable, str(ANALYZER), str(directory)]
    if args.method == "rare_event":
        analysis_command.append("--allow-rare-event-samples")
    print(f"Sampler: {args.method}", flush=True)
    print(f"Run directory: {directory}", flush=True)
    print(shlex.join(command), flush=True)
    if args.dry_run:
        print(shlex.join(analysis_command))
        return

    directory.mkdir(parents=True, exist_ok=False)
    (directory / "command.txt").write_text(
        "Contributor: OpenAI Codex.\n\n" + shlex.join(command) + "\n\n" + shlex.join(analysis_command) + "\n")
    manifest = {"contributor": "OpenAI Codex", "started_utc": datetime.now(timezone.utc).isoformat(),
                "decoder_command": command, "analysis_command": analysis_command,
                "working_directory": str(REPOSITORY), "sampler": args.method, "status": "running"}
    manifest_path = directory / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    start = time.monotonic()
    try:
        print("Sampling (full output in run.log)...", flush=True)
        run_logged(command, directory / "run.log", progress=True)
        manifest["sampling_seconds"] = time.monotonic() - start
        print("Validating and analyzing recorded syndromes...", flush=True)
        run_logged(analysis_command, directory / "analysis.log")
        summary = json.loads((directory / "analysis_summary.json").read_text())
        if ((args.method == "monte_carlo" and summary["shots"] != args.shots)
                or len(summary["rank_shots"]) != args.cores):
            raise RuntimeError("Recorded shot/rank counts differ from the requested experiment")
        ler_lines = [line for line in (directory / "run.log").read_text().splitlines()
                     if line.startswith("Logical error rate:")]
        if not ler_lines:
            raise RuntimeError("qudec did not report the logical error rate")
        manifest["logical_error_rate"] = float(ler_lines[-1].split(":", 1)[1])
        report = format_summary(directory, summary)
        report = ler_lines[-1] + "\n" + report
        (directory / "experiment_summary.txt").write_text(report)
        manifest["status"] = "complete"
        print(report, end="")
    except (OSError, RuntimeError, ValueError, KeyboardInterrupt) as error:
        manifest["status"] = "failed"
        manifest["error"] = str(error) or type(error).__name__
        print(f"Experiment failed: {error}\nFiles preserved in {directory}", file=sys.stderr)
        return_code = 1
    else:
        return_code = 0
    finally:
        manifest["total_seconds"] = time.monotonic() - start
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return return_code


if __name__ == "__main__":
    sys.exit(main())
