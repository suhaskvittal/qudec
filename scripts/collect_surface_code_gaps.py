#!/usr/bin/env python3
"""Collect UF and uncorrelated PyMatching cluster-gap histograms with 64 MPI ranks."""

import argparse
from pathlib import Path
import subprocess

from _surface_code_collection import (
    REPOSITORY, launcher_command, mpi_binary, positive_int, run_to_file,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=positive_int, default=1_000_000,
                        help="global shots per distance and decoder (default: 1,000,000)")
    parser.add_argument("--p", type=float, default=0.001, help="physical error rate (default: 0.001)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--cores", type=positive_int, default=64)
    parser.add_argument("--qudec", type=Path, help="MPI-enabled executable (default: build/qudec)")
    parser.add_argument("--mpi-launcher", default="mpirun")
    parser.add_argument("--output-root", type=Path, default=REPOSITORY / "out")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not 0 < args.p < 0.5 or not 0 <= args.seed <= 2**63 - args.cores:
        parser.error("--p must be between 0 and 0.5, and --seed must leave room for MPI rank offsets")
    try:
        binary = mpi_binary(args)
        prefix = launcher_command(args)
        for decoder, folder in (("uf", "union_find_gap"), ("pymatching", "pymatching_gap")):
            for distance in [5, 7, 9, 11, 31]:
                if distance == 31 and decoder == 'uf':
                    continue
                command = prefix + [str(binary), decoder, "-d", str(distance), "-r", str(3 * distance),
                    "-p", str(args.p), "--code", "sc",
                    "--cluster-gap", "--method", "monte_carlo", "--seed", str(args.seed),
                    "--mc-max-samples", str(args.shots), "--mc-stop-at-errors", str(args.shots + 1)]
                output = args.output_root / folder / f"d{distance}.out"
                run_to_file(command, output, "CLUSTER_GAP_DB", args.shots, args.dry_run)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")


if __name__ == "__main__":
    main()
