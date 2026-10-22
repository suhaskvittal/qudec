#!/usr/bin/env python3
"""Prepare large, fixed-DEM affinity-gap datasets for neural-network training.

Each invocation prepares one raw exporter directory into a standalone dataset
directory. Arrays are written as NumPy files; triangle features are float16 on
disk (the trainer should cast them to float32 when loading). Feature workers
read and write disjoint contiguous row ranges, and the protocol file is written
last as the completion marker.
"""

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
import csv
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from affinity_matrix_layout import LayoutOverflow, build_candidates, deterministic_placement, packed_placement


RAW_FILES = ("metadata.json", "circuit.stim", "raw.dem", "decomposed.dem",
             "graph_edges.csv", "coords.csv", "affinity.f64", "hops.u8",
             "syndromes.u8", "labels.f64")
FEATURE_CHOICES = ("syndrome", "summary", "packed_triangle", "hashed_triangle")
SUMMARY_NAMES = [
    "weight", "weight_parity", "mean_nearest_log_distance", "max_nearest_log_distance",
    "mean_nearest_second_margin", "min_nearest_second_margin", "mean_effective_partners",
    "max_effective_partners", "isolated_rows", "mean_pair_log_distance", "min_pair_log_distance",
    "max_pair_log_distance", "mean_shared_boundary_log_distance", "min_shared_boundary_log_distance",
    "max_shared_boundary_log_distance", "memory_weight", "opposite_weight",
    "min_boundary0_distance", "min_boundary1_distance", "mean_boundary0_distance",
    "mean_boundary1_distance", "min_boundary_asymmetry", "mean_boundary_asymmetry",
    "max_boundary_asymmetry",
]

_WORK = {}


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_coordinates(path, n):
    coords = np.full((n, 3), np.nan, dtype=np.float64)
    with (path / "coords.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            i = int(row["detector_id"])
            if not 0 <= i < n:
                raise ValueError(f"Invalid detector id in coords.csv: {i}")
            coords[i] = [float(row[k]) for k in ("x", "y", "t")]
    if not np.isfinite(coords).all():
        raise ValueError("coords.csv is missing finite coordinates for one or more detectors")
    return coords


def graph_geometry(path, n):
    """Return two physical-boundary distances and memory-boundary relevance."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import dijkstra

    probabilities = {}
    with (path / "graph_edges.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            u, v = int(row["detector_u"]), int(row["detector_v"])
            if v == -1:
                v = n + int(row["boundary_class"])
            key = tuple(sorted((u, v)))
            p = 1. / (1. + np.exp(float(row["weight_natural_log"])))
            old = probabilities.get(key, 0.)
            probabilities[key] = old + p - 2. * old * p
    rows, cols, weights = [], [], []
    # Boundary vertices have outgoing edges only, so a path cannot leave one
    # physical boundary and transit through another component.
    for (u, v), p in probabilities.items():
        weight = np.log((1. - p) / p)
        rows.append(v); cols.append(u); weights.append(weight)
        if v < n:
            rows.append(u); cols.append(v); weights.append(weight)
    graph = csr_matrix((weights, (rows, cols)), shape=(n + 2, n + 2))
    distances = dijkstra(graph, directed=True, indices=[n, n + 1])[:, :n]
    return distances, np.isfinite(distances[1])


def summary_features(ids, affinity, boundary, relevant, scale):
    """Match affinity_gap_study.summary_features without importing its ML stack."""
    k = len(ids)
    if not k:
        return np.zeros(len(SUMMARY_NAMES), dtype=np.float32)
    raw = affinity[np.ix_(ids, ids)]
    distances = np.where(raw > 0, -np.log(np.maximum(raw, 1e-300)) / scale, 2.)
    np.fill_diagonal(distances, 2.)
    nearest = np.sort(distances, axis=1)[:, :2]
    if k == 1:
        nearest = np.full((1, 2), 2.)
    margin = nearest[:, 1] - nearest[:, 0]
    sums = raw.sum(axis=1)
    squared_sums = (raw ** 2).sum(axis=1)
    effective = np.divide(sums ** 2, squared_sums, out=np.zeros(k), where=squared_sums > 0)
    pairs = distances[np.triu_indices(k, 1)]
    pairs = pairs[pairs < 2.]
    shared = -np.log(np.maximum(affinity[ids, -1], 1e-300)) / scale
    base = [k, k % 2, nearest[:, 0].mean(), nearest[:, 0].max(), margin.mean(), margin.min(),
            effective.mean(), effective.max(), np.count_nonzero(sums == 0),
            pairs.mean() if len(pairs) else 0., pairs.min() if len(pairs) else 0.,
            pairs.max() if len(pairs) else 0., shared.mean(), shared.min(), shared.max()]
    memory = ids[relevant[ids]]
    if len(memory):
        b = boundary[:, memory] / scale
        asym = np.abs(b[0] - b[1])
        extra = [len(memory), k - len(memory), b[0].min(), b[1].min(),
                 b[0].mean(), b[1].mean(), asym.min(), asym.mean(), asym.max()]
    else:
        extra = [0., k, 0., 0., 0., 0., 0., 0., 0.]
    return np.asarray(base + extra, dtype=np.float32)


def triangle_row(ids, slots, pair_strength, pair_connected, boundary, size):
    """Construct existing strength/mask upper-triangle layout without dense scatter loops."""
    side = size * (size + 1) // 2
    result = np.zeros(2 * side, dtype=np.float32)
    b0, b1 = size - 2, size - 1
    boundary_diag = np.asarray([b0 * size - b0 * (b0 - 1) // 2,
                                b1 * size - b1 * (b1 - 1) // 2], dtype=np.int64)
    result[boundary_diag] = 1.
    result[side + boundary_diag] = 1.

    k = len(ids)
    if not k:
        return result
    node_diag = slots * size - slots * (slots - 1) // 2
    result[node_diag] = 1.
    result[side + node_diag] = 1.

    ii, jj = np.triu_indices(k)
    si, sj = slots[ii], slots[jj]
    lo, hi = np.minimum(si, sj), np.maximum(si, sj)
    triangle_index = lo * size - lo * (lo - 1) // 2 + (hi - lo)
    connected = pair_connected[ids[ii], ids[jj]]
    values = pair_strength[ids[ii], ids[jj]]
    diag = ii == jj
    values[diag] = 1.
    connected[diag] = True
    result[triangle_index] = values
    result[side + triangle_index] = connected.astype(np.float32)

    for boundary_index, boundary_slot in enumerate((b0, b1)):
        distances = boundary[boundary_index, ids]
        good = np.isfinite(distances)
        strength = np.zeros(k, dtype=np.float32)
        strength[good] = np.clip(1. - distances[good] / _WORK["scale"], 0., 1.)
        lo = np.minimum(slots, boundary_slot)
        hi = np.maximum(slots, boundary_slot)
        idx = lo * size - lo * (lo - 1) // 2 + (hi - lo)
        result[idx] = strength
        result[side + idx] = good.astype(np.float32)
    return result


def init_worker(config):
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    _WORK.update(config)
    n = config["n"]
    _WORK["syndromes"] = np.memmap(config["source"] / "syndromes.u8", mode="r", dtype="u1",
                                    shape=(config["shots"], n))
    _WORK["affinity"] = np.memmap(config["source"] / "affinity.f64", mode="r", dtype="<f8",
                                   shape=(n + 1, n + 1))
    _WORK["output_arrays"] = {}
    for name in config["features"]:
        path = config["output"] / f"{name}.npy"
        _WORK["output_arrays"][name] = np.load(path, mmap_mode="r+")
    if "summary" in config["features"]:
        _WORK["summary_out"] = _WORK["output_arrays"]["summary"]
    if "packed_triangle" in config["features"] or "hashed_triangle" in config["features"]:
        affinity = np.asarray(_WORK["affinity"], dtype=np.float64)
        raw = (affinity[:-1, :-1] + affinity[:-1, :-1].T) * .5
        connected = raw > 0
        strength = np.zeros_like(raw)
        strength[connected] = np.clip(1. + np.log(np.maximum(raw[connected], 1e-300)) / config["scale"], 0., 1.)
        _WORK["pair_strength"] = strength
        _WORK["pair_connected"] = connected


def process_chunk(bounds):
    start, stop = bounds
    cfg = _WORK
    size, n = cfg["size"], cfg["n"]
    syndromes = cfg["syndromes"][start:stop]
    output_arrays = cfg["output_arrays"]
    max_weight = 0
    for offset, syndrome in enumerate(syndromes):
        row = start + offset
        ids = np.flatnonzero(syndrome)
        max_weight = max(max_weight, len(ids))
        if "summary" in output_arrays:
            output_arrays["summary"][row] = summary_features(
                ids, cfg["affinity"], cfg["boundary"], cfg["relevant"], cfg["scale"])
        if "packed_triangle" in output_arrays:
            try:
                slots = packed_placement(ids, size)
            except LayoutOverflow as exc:
                raise LayoutOverflow(f"row {row}, syndrome weight {len(ids)}: {exc}") from exc
            output_arrays["packed_triangle"][row] = triangle_row(
                ids, slots, cfg["pair_strength"], cfg["pair_connected"], cfg["boundary"], size).astype(np.float16)
        if "hashed_triangle" in output_arrays:
            try:
                slots = deterministic_placement(ids, cfg["candidates"], size)
            except LayoutOverflow as exc:
                raise LayoutOverflow(f"row {row}, syndrome weight {len(ids)}: {exc}") from exc
            output_arrays["hashed_triangle"][row] = triangle_row(
                ids, slots, cfg["pair_strength"], cfg["pair_connected"], cfg["boundary"], size).astype(np.float16)
    return start, stop, max_weight


def copy_base_arrays(source, output, shots, n, chunk_rows=65536):
    raw_syndromes = np.memmap(source / "syndromes.u8", mode="r", dtype="u1", shape=(shots, n))
    raw_labels = np.memmap(source / "labels.f64", mode="r", dtype="<f8", shape=(shots, 6))
    syndromes = np.lib.format.open_memmap(output / "syndrome.npy", mode="w+", dtype="uint8", shape=(shots, n))
    labels = np.lib.format.open_memmap(output / "labels.npy", mode="w+", dtype="float64", shape=(shots, 6))
    for start in range(0, shots, chunk_rows):
        stop = min(shots, start + chunk_rows)
        syndromes[start:stop] = raw_syndromes[start:stop]
        labels[start:stop] = raw_labels[start:stop]
        if not np.isfinite(labels[start:stop]).all() or (labels[start:stop, :2] < 0).any():
            raise ValueError(f"Invalid finite/nonnegative gap labels in rows {start}:{stop}")
        print(f"copied source arrays: {stop}/{shots}", flush=True)
    syndromes.flush(); labels.flush()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True, help="Raw export directory")
    parser.add_argument("--output", type=Path, required=True, help="New prepared dataset directory")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--features", nargs="+", choices=FEATURE_CHOICES,
                        default=list(FEATURE_CHOICES))
    parser.add_argument("--size", type=int, default=48)
    parser.add_argument("--hash-seed", type=int, default=17)
    parser.add_argument("--chunk-rows", type=int, default=4096)
    args = parser.parse_args()
    source = args.dataset.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {output}")
    if args.workers < 1 or args.chunk_rows < 1:
        raise ValueError("--workers and --chunk-rows must be positive")
    features = list(dict.fromkeys(args.features))
    if not features:
        raise ValueError("At least one feature must be selected")
    missing = [name for name in RAW_FILES if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Raw dataset is missing: {', '.join(missing)}")

    meta = json.loads((source / "metadata.json").read_text())
    shots, n = int(meta["shots"]), int(meta["num_detectors"])
    if shots < 1 or n < 1:
        raise ValueError("Source metadata must have positive shots and detector count")
    expected_sizes = {
        "syndromes.u8": shots * n,
        "labels.f64": shots * 6 * 8,
        "affinity.f64": (n + 1) * (n + 1) * 8,
        "hops.u8": (n + 1) * (n + 1),
    }
    for name, expected in expected_sizes.items():
        actual = (source / name).stat().st_size
        if actual != expected:
            raise ValueError(f"{name} has {actual} bytes; expected {expected}")

    output.mkdir(parents=True, exist_ok=False)
    shutil_metadata = json.loads(json.dumps(meta))
    (output / "metadata.json").write_text(json.dumps(shutil_metadata, indent=2) + "\n")

    boundary, relevant = graph_geometry(source, n)
    coordinates = load_coordinates(source, n)
    candidates = build_candidates(coordinates, relevant, size=args.size, hashseed=args.hash_seed)
    scale = float(meta["distance"]) * np.log(1. / float(meta["p"]))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Invalid distance/probability scale in source metadata")

    # Copy the fixed affinity matrix once, symmetrize before transforming, and
    # use this same transform for both triangle layouts.
    affinity = np.memmap(source / "affinity.f64", mode="r", dtype="<f8", shape=(n + 1, n + 1))
    raw = (np.asarray(affinity[:-1, :-1], dtype=np.float64) +
           np.asarray(affinity[:-1, :-1].T, dtype=np.float64)) * .5
    pair_connected = raw > 0
    pair_strength = np.zeros_like(raw)
    pair_strength[pair_connected] = np.clip(1. + np.log(np.maximum(raw[pair_connected], 1e-300)) / scale, 0., 1.)
    del raw

    side = args.size * (args.size + 1) // 2
    shapes = {"syndrome": (shots, n), "summary": (shots, len(SUMMARY_NAMES)),
              "packed_triangle": (shots, 2 * side), "hashed_triangle": (shots, 2 * side)}
    dtypes = {"syndrome": np.uint8, "summary": np.float32,
              "packed_triangle": np.float16, "hashed_triangle": np.float16}
    for feature in features:
        path = output / f"{feature}.npy"
        if feature != "syndrome":
            np.lib.format.open_memmap(path, mode="w+", dtype=dtypes[feature], shape=shapes[feature]).flush()
        else:
            # Syndrome is copied in bounded blocks and converted without a full
            # sample list or a duplicate 1M-row in-memory array.
            pass
    copy_base_arrays(source, output, shots, n)
    config = {
        "source": source, "output": output, "shots": shots, "n": n,
        "features": [f for f in features if f != "syndrome"], "size": args.size,
        "scale": scale, "boundary": boundary, "relevant": relevant,
        "candidates": candidates, "pair_strength": pair_strength,
        "pair_connected": pair_connected,
    }
    tasks = [(start, min(shots, start + args.chunk_rows))
             for start in range(0, shots, args.chunk_rows)]
    max_weight = 0
    if config["features"]:
        with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker,
                                 initargs=(config,)) as pool:
            futures = [pool.submit(process_chunk, bounds) for bounds in tasks]
            completed = 0
            for future in as_completed(futures):
                start, stop, weight = future.result()
                max_weight = max(max_weight, weight)
                completed += stop - start
                print(f"prepared feature rows: {completed}/{shots}; max syndrome weight={max_weight}", flush=True)
    else:
        max_weight = int(np.max(np.asarray(np.load(output / "syndrome.npy", mmap_mode="r")).sum(axis=1)))

    # The source has one fixed DEM, affinity/hop tables and coordinate layout;
    # include hashes of all raw inputs so training runs can verify provenance.
    source_hashes = {name: sha256_file(source / name) for name in RAW_FILES}
    layout = {
        "size": args.size,
        "hash_seed": args.hash_seed,
        "candidate_table_shape": list(candidates.shape),
        "candidate_slots_per_detector": int(candidates.shape[1]),
        "physical_slots": [0, args.size - 3],
        "boundary_slots": [args.size - 2, args.size - 1],
        "upper_triangle_order": "row-major upper triangle including diagonal; strength then connectivity mask",
        "triangle_channels": ["normalized_affinity_strength", "connected_or_active_mask"],
        "triangle_width": int(2 * side),
        "triangle_storage_dtype": "float16; load and cast to float32 for model input",
        "triangle_quantization": "IEEE float16 rounding; maximum absolute quantization error is bounded by half-precision rounding for each stored value",
        "hash_assignment": "coordinate-keyed deterministic 2-choice bucket layout plus two shared physical stash candidates; at most 10 candidates per detector; detector IDs are not inputs",
        "overflow": "hard error; active detectors are never dropped or merged",
    }
    protocol = {
        "schema": 1,
        "source_dataset": str(source),
        "source_metadata": meta,
        "source_sha256": source_hashes,
        "shots": shots,
        "num_detectors": n,
        "features": features,
        "feature_arrays": {name: {"file": f"{name}.npy", "shape": list(shapes[name]),
                                  "dtype": np.dtype(dtypes[name]).name}
                           for name in features if name != "syndrome"},
        "labels": {"file": "labels.npy", "shape": [shots, 6], "dtype": "float64",
                   "columns": meta.get("labels.f64", {}).get("columns", [])},
        "syndrome": {"file": "syndrome.npy", "shape": [shots, n], "dtype": "uint8"},
        "summary_names": SUMMARY_NAMES,
        "summary_shape": [shots, len(SUMMARY_NAMES)],
        "triangle_layout": layout,
        "log_distance_scale": scale,
        "boundary_distances": "separate graph distances to physical boundary classes 0 and 1; boundaries are isolated endpoints",
        "max_syndrome_weight": max_weight,
        "workers": args.workers,
        "chunk_rows": args.chunk_rows,
        "label_inputs_excluded_from_features": True,
    }
    text = json.dumps(protocol, indent=2, allow_nan=False) + "\n"
    (output / "datasetprotocol.json").write_text(text)
    (output / "protocol.json").write_text(text)
    print(json.dumps({"output": str(output), "shots": shots, "num_detectors": n,
                      "features": features, "max_syndrome_weight": max_weight,
                      "protocol": str(output / "datasetprotocol.json")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
