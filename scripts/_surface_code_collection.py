"""Shared MPI runner for surface-code histogram collection scripts."""

import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess


REPOSITORY = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = REPOSITORY / "build" / "qudec"


def positive_int(value):
    number = int(value)
    if number < 1:
        raise ValueError("expected a positive integer")
    return number


def mpi_binary(args):
    """Locate and validate an already compiled MPI executable."""
    binary = (args.qudec or DEFAULT_BINARY).expanduser().resolve()
    if not args.dry_run and not (binary.is_file() and os.access(binary, os.X_OK)):
        raise ValueError(f"qudec executable not found: {binary}")
    cache = binary.parent / "CMakeCache.txt"
    if cache.is_file():
        if "ENABLE_MPI:BOOL=ON" not in (line.upper() for line in cache.read_text().splitlines()):
            raise ValueError(f"qudec is not an MPI build: {binary}")
    elif not args.dry_run:
        dependencies = subprocess.run(["ldd", str(binary)], capture_output=True,
                                  text=True, check=True).stdout
        if "libmpi" not in dependencies.lower():
            raise ValueError(f"qudec is not linked to MPI: {binary}")
    return binary


def launcher_command(args):
    launcher = shutil.which(args.mpi_launcher)
    if launcher is None:
        if args.dry_run:
            launcher = args.mpi_launcher
        else:
            raise ValueError(f"MPI launcher not found: {args.mpi_launcher}")
    return [launcher, "-np", str(args.cores)]


def run_to_file(command, output, expected_histogram, shots, dry_run):
    print(f"{shlex.join(command)} > {shlex.quote(str(output))}", flush=True)
    if dry_run:
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    with output.open("x") as stream:
        subprocess.run(command, cwd=REPOSITORY, stdout=stream, stderr=subprocess.STDOUT, check=True)
    report = output.read_text()
    header = re.escape(expected_histogram)
    match = re.search(rf"(?m)^{header}\s*=+\n(?:mean\s*=.*?count\s*=\s*(\d+)|count\s*=\s*0)\s*$", report)
    if not match:
        raise ValueError(f"{output}: expected histogram {expected_histogram} was not printed")
    observed = int(match.group(1) or 0)
    if observed != shots:
        raise ValueError(f"{output}: expected {shots} histogram samples, found {observed}")
