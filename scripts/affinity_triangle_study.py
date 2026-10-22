#!/usr/bin/env python3
"""Prepare and evaluate packed or hashed upper-triangle affinity features.

This script never trains models. ``prepare`` writes train/validation/test
feature arrays compatible with ``affinity_gap_study.py``'s MLP trainer;
``evaluate`` applies already-trained packed/hash triangle MLPs and the raw
fixed-detector syndrome MLP control to the prepared test split.
"""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import affinity_gap_study as study
from affinity_matrix_layout import (LayoutOverflow, build_candidates,
                                    deterministic_placement, packed_placement,
                                    triangle_features)


TOPOLOGY_FILES = ("circuit.stim", "raw.dem", "decomposed.dem", "graph_edges.csv", "coords.csv")
MATCH_FIELDS = ("distance", "rounds", "p", "num_detectors", "memory_basis",
                "include_opposite_basis_detectors")
SPLITS = ("train", "validation", "test")
FEATURES = ("packed_triangle", "hashed_triangle")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_protocol(reference):
    protocol_path = reference / "protocol.json"
    if not protocol_path.is_file():
        raise ValueError(f"Reference protocol missing: {protocol_path}")
    protocol = json.loads(protocol_path.read_text())
    if len(protocol.get("datasets", [])) != 3 or len(protocol.get("metadata", [])) != 3:
        raise ValueError("Reference protocol must describe train, validation, and test splits")
    if not protocol.get("topology_sha256"):
        raise ValueError("Reference protocol lacks topology hashes; use a current prepared dataset")
    return protocol


def verify_dataset(path, expected_metadata, topology_hashes, old_seeds, require_new_seed=False):
    path = Path(path).expanduser().resolve()
    if not (path / "metadata.json").is_file():
        raise ValueError(f"Dataset metadata missing: {path / 'metadata.json'}")
    meta = json.loads((path / "metadata.json").read_text())
    mismatch = [key for key in MATCH_FIELDS if meta.get(key) != expected_metadata.get(key)]
    if mismatch:
        raise ValueError(f"Dataset {path} differs from reference in: {', '.join(mismatch)}")
    if require_new_seed and meta.get("seed") in old_seeds:
        raise ValueError(f"Override test seed {meta.get('seed')} repeats a reference split seed")
    for name in TOPOLOGY_FILES:
        file_path = path / name
        if not file_path.is_file():
            raise ValueError(f"Dataset topology file missing: {file_path}")
        if digest(file_path) != topology_hashes.get(name):
            raise ValueError(f"Dataset topology hash mismatch: {name}")
    return path, meta


def load_coords(path, n):
    coords = np.full((n, 3), np.nan, dtype=np.float64)
    with (path / "coords.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            detector = int(row["detector_id"])
            if detector < 0 or detector >= n:
                raise ValueError(f"Invalid detector ID in coords.csv: {detector}")
            coords[detector] = [float(row[name]) for name in ("x", "y", "t")]
    return coords


def _hardlink(source, target):
    source = Path(source)
    if not source.is_file():
        raise ValueError(f"Missing reusable prepared artifact: {source}")
    os.link(source, target)


def verify_reuse(reuse, reference, candidate_slots, size, hash_seed):
    protocol_path = reuse / "protocol.json"
    if not protocol_path.is_file():
        raise ValueError(f"Reusable prepared protocol missing: {protocol_path}")
    old = json.loads(protocol_path.read_text())
    layout = old.get("triangle_layout", {})
    if layout.get("size") != size or layout.get("hash_seed") != hash_seed:
        raise ValueError("Reusable features use a different size or hash seed")
    existing = np.load(reuse / "candidate_slots.npy", allow_pickle=False)
    if not np.array_equal(existing, candidate_slots):
        raise ValueError("Reusable candidate-slot table differs from current coordinates/layout")
    current = load_protocol(reference)
    if old.get("log_distance_scale") != current.get("log_distance_scale"):
        raise ValueError("Reusable features use a different log-distance scale")
    for key in ("train", "validation"):
        if old["metadata"][SPLITS.index(key)] != current["metadata"][SPLITS.index(key)]:
            raise ValueError(f"Reusable {key} dataset metadata differs from the reference")
    for name in TOPOLOGY_FILES:
        if old.get("topology_sha256", {}).get(name) != current["topology_sha256"].get(name):
            raise ValueError(f"Reusable topology differs from current reference: {name}")


def prepare(args):
    reference = args.reference.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if output.exists():
        raise ValueError(f"Output path already exists: {output}")
    protocol = load_protocol(reference)
    original_paths = [Path(p).expanduser().resolve() for p in protocol["datasets"]]
    old_seeds = {meta["seed"] for meta in protocol["metadata"]}
    paths = list(original_paths)
    if args.test is not None:
        paths[2] = args.test.expanduser().resolve()
    verified = [verify_dataset(path, protocol["metadata"][0], protocol["topology_sha256"],
                               old_seeds, require_new_seed=(i == 2 and args.test is not None))
                for i, path in enumerate(paths)]
    metas = [item[1] for item in verified]
    if len({meta["seed"] for meta in metas}) != 3:
        raise ValueError("Train, validation, and test datasets must have distinct seeds")

    n = metas[0]["num_detectors"]
    boundary_distances, relevant = study.graph_geometry(paths[0], n)
    coords = load_coords(paths[0], n)
    candidates = build_candidates(coords, relevant, size=args.size, hashseed=args.hash_seed)
    reference_affinity = np.fromfile(paths[0] / "affinity.f64", dtype="<f8")
    reference_hops = np.fromfile(paths[0] / "hops.u8", dtype="u1")
    table_shape = (n + 1, n + 1)
    if reference_affinity.size != np.prod(table_shape) or reference_hops.size != np.prod(table_shape):
        raise ValueError("Reference affinity or hop lookup table has an invalid byte length")
    reference_affinity = reference_affinity.reshape(table_shape)
    reference_hops = reference_hops.reshape(table_shape)
    if args.reuse_training is not None:
        verify_reuse(args.reuse_training.expanduser().resolve(), reference, candidates,
                     args.size, args.hash_seed)

    output.mkdir(parents=True, exist_ok=False)
    np.save(output / "candidate_slots.npy", candidates)
    scale = float(protocol["log_distance_scale"])
    for split_index, split in enumerate(SPLITS):
        dataset_path = paths[split_index]
        if split_index < 2 and args.reuse_training is not None:
            reuse = args.reuse_training.expanduser().resolve()
            for suffix in ("labels", "summary", "syndrome", "placements", *FEATURES):
                _hardlink(reuse / f"{split}_{suffix}.npy", output / f"{split}_{suffix}.npy")
            print(f"{split}: reused prepared arrays ({metas[split_index]['shots']} shots)", flush=True)
            continue

        meta, syndromes, labels, affinity, hops = study.read_dataset(dataset_path)
        if not np.array_equal(affinity, reference_affinity) or not np.array_equal(hops, reference_hops):
            raise ValueError(f"Affinity lookup tables differ for {split} dataset")

        summaries = np.empty((len(syndromes), len(study.SUMMARY_NAMES)), dtype=np.float32)
        feature_arrays = {feature: np.lib.format.open_memmap(
            output / f"{split}_{feature}.npy", mode="w+", dtype="float32",
            shape=(len(syndromes), 2 * args.size * (args.size + 1) // 2)) for feature in FEATURES}
        placement_rows = []
        for row, syndrome in enumerate(syndromes):
            ids = np.flatnonzero(syndrome)
            try:
                packed_slots = packed_placement(ids, args.size)
                hashed_slots = deterministic_placement(ids, candidates, args.size)
            except LayoutOverflow as exc:
                raise LayoutOverflow(f"{split} shot {row}, weight {len(ids)}: {exc}; increase --size") from exc
            summaries[row] = study.summary_features(ids, affinity, boundary_distances,
                                                     relevant, scale)
            feature_arrays["packed_triangle"][row] = triangle_features(
                ids, affinity, boundary_distances, scale, args.size, packed_slots)
            feature_arrays["hashed_triangle"][row] = triangle_features(
                ids, affinity, boundary_distances, scale, args.size, hashed_slots)
            placement_rows.extend((row, int(detector), int(packed_slot), int(hashed_slot))
                                  for detector, packed_slot, hashed_slot
                                  in zip(ids, packed_slots, hashed_slots))
        for feature_array in feature_arrays.values():
            feature_array.flush()
        np.save(output / f"{split}_labels.npy", labels)
        np.save(output / f"{split}_syndrome.npy", syndromes.astype(np.float32))
        np.save(output / f"{split}_summary.npy", summaries)
        np.save(output / f"{split}_placements.npy",
                np.asarray(placement_rows, dtype=np.int32).reshape((-1, 4)))
        print(f"{split}: prepared {len(syndromes)} shots; max syndrome weight="
              f"{int(syndromes.sum(axis=1).max())}; wrote packed and hashed triangles",
              flush=True)

    # When training arrays were reused, load training syndromes for the fixed
    # novelty reference; source artifacts remain untouched.
    train_syndromes = np.load(output / "train_syndrome.npy", mmap_mode="r")
    train_hashes = {hashlib.sha256(np.packbits(row.astype(np.uint8)).tobytes()).digest()
                    for row in train_syndromes}
    for split_index, split in ((1, "validation"), (2, "test")):
        syndrome_file = output / f"{split}_syndrome.npy"
        syndromes = np.load(syndrome_file, mmap_mode="r")
        novel = np.asarray([hashlib.sha256(np.packbits(row.astype(np.uint8)).tobytes()).digest()
                            not in train_hashes for row in syndromes], dtype=bool)
        np.save(output / f"{split}_novel.npy", novel)
    topology = protocol["topology_sha256"]
    updated = json.loads(json.dumps(protocol))
    updated["datasets"] = [str(path) for path in paths]
    updated["metadata"] = metas
    updated["triangle_layout"] = {
        "size": args.size, "hash_seed": args.hash_seed,
        "candidate_slots_file": "candidate_slots.npy",
        "candidate_table_shape": list(candidates.shape),
        "boundaries": [args.size - 2, args.size - 1],
        "feature_names": list(FEATURES),
        "feature_shape": [2 * args.size * (args.size + 1) // 2],
        "triangle_order": "row-major upper triangle including diagonal; strength then mask",
        "placements_file_pattern": "{split}_placements.npy rows are [shot_index, detector_id, packed_slot, hashed_slot]",
        "overflow": "hard error; no active detector is dropped or merged",
        "reuse_training": str(args.reuse_training.resolve()) if args.reuse_training else None,
    }
    updated["topology_sha256"] = topology
    updated["feature_inputs"] = "syndrome and fixed DEM only; labels and decoder outputs excluded"
    updated["boundary_detector_slots"] = [args.size - 2, args.size - 1]
    updated["test_novel_shots"] = int(np.load(output / "test_novel.npy").sum())
    (output / "protocol.json").write_text(json.dumps(updated, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(output), "shots": {split: int(meta["shots"]) for split, meta in zip(SPLITS, metas)},
                      "size": args.size, "hash_seed": args.hash_seed,
                      "test_novel_shots": updated["test_novel_shots"]}, indent=2))


def evaluate(args):
    from evaluate_affinity_gap_holdout import parity_scores, predict_torch

    prepared = args.prepared.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if output.exists():
        raise ValueError(f"Output path already exists: {output}")
    protocol_path = prepared / "protocol.json"
    if not protocol_path.is_file():
        raise ValueError(f"Prepared protocol missing: {protocol_path}")
    protocol = json.loads(protocol_path.read_text())
    if not protocol.get("triangle_layout"):
        raise ValueError("Prepared data has no triangle layout metadata")
    if len(protocol.get("metadata", [])) != 3 or not protocol.get("topology_sha256"):
        raise ValueError("Prepared protocol lacks train/validation/test metadata or topology hashes")
    output.mkdir(parents=True, exist_ok=False)
    labels = np.load(prepared / "test_labels.npy")
    validation_labels = np.load(prepared / "validation_labels.npy")
    summaries = np.load(prepared / "test_summary.npy", mmap_mode="r")
    novel = np.load(prepared / "test_novel.npy")
    feature_arrays = {feature: np.load(prepared / f"test_{feature}.npy", mmap_mode="r")
                      for feature in FEATURES}
    feature_arrays["syndrome"] = np.load(prepared / "test_syndrome.npy", mmap_mode="r")
    feature_arrays["summary"] = summaries
    import torch
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    manifest = {"prepared": str(prepared), "models": []}
    for run in args.runs:
        run = run.expanduser().resolve()
        metrics_path = run / "metrics.json"
        if not metrics_path.is_file():
            raise ValueError(f"Missing model metrics: {metrics_path}")
        metrics = json.loads(metrics_path.read_text())
        model_name = metrics.get("model", run.name)
        accepted_models = ("mlp_packed_triangle", "mlp_hashed_triangle", "mlp_syndrome", "mlp_boundary")
        if model_name not in accepted_models:
            raise ValueError(f"Unsupported triangle-study control model {model_name!r} in {run}")
        if run.name in {entry["run_name"] for entry in manifest["models"]}:
            raise ValueError(f"Duplicate source run name: {run.name}")
        validation_path = run / "validation_predictions.npy"
        if not validation_path.is_file():
            raise ValueError(f"Missing validation predictions: {validation_path}")
        validation_predictions = np.load(validation_path)
        if validation_predictions.shape != (len(validation_labels), 2):
            raise ValueError(f"Validation prediction shape mismatch for {run}")
        checkpoint_path = run / "model.pt"
        if not checkpoint_path.is_file():
            raise ValueError(f"Missing model checkpoint: {checkpoint_path}")
        checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
        checkpoint_protocol = checkpoint.get("protocol")
        if not isinstance(checkpoint_protocol, dict):
            raise ValueError(f"Checkpoint has no prepared-data protocol: {checkpoint_path}")
        expected_layout = protocol["triangle_layout"]
        checkpoint_layout = checkpoint_protocol.get("triangle_layout")
        if model_name != "mlp_boundary":
            if not isinstance(checkpoint_layout, dict):
                raise ValueError(f"Checkpoint has no triangle layout metadata: {checkpoint_path}")
            for key in ("size", "hash_seed", "boundaries", "feature_shape"):
                if checkpoint_layout.get(key) != expected_layout.get(key):
                    raise ValueError(f"Checkpoint triangle layout mismatch for {key}: {checkpoint_path}")
        for key in ("topology_sha256", "log_distance_scale"):
            if checkpoint_protocol.get(key) != protocol.get(key):
                raise ValueError(f"Checkpoint prepared-data mismatch for {key}: {checkpoint_path}")
        for index, split in enumerate(("train", "validation")):
            if (len(checkpoint_protocol.get("metadata", [])) <= index or
                    checkpoint_protocol["metadata"][index] != protocol["metadata"][index]):
                raise ValueError(f"Checkpoint {split} dataset metadata mismatch: {checkpoint_path}")
        expected_feature = {"mlp_boundary": "summary"}.get(
            model_name, model_name.removeprefix("mlp_"))
        if checkpoint.get("feature") != expected_feature:
            raise ValueError(f"Checkpoint feature {checkpoint.get('feature')!r} does not match {model_name!r}")
        if expected_feature in FEATURES:
            expected_shape = tuple(expected_layout.get("feature_shape", ()))
        elif expected_feature == "syndrome":
            expected_shape = (int(protocol["metadata"][0]["num_detectors"]),)
        elif expected_feature == "summary":
            expected_shape = (len(study.SUMMARY_NAMES),)
        else:
            raise ValueError(f"Unsupported checkpoint feature: {expected_feature}")
        if expected_shape != tuple(feature_arrays[expected_feature].shape[1:]):
            raise ValueError(f"Prepared feature array shape disagrees with its protocol: {expected_feature}")

        predictions = predict_torch(checkpoint_path, feature_arrays, model_name)
        if predictions.shape != (len(labels), 2) or not np.isfinite(predictions).all():
            raise ValueError(f"Invalid test predictions for {run}")
        model_output = output / run.name
        model_output.mkdir()
        np.save(model_output / "test_predictions.npy", predictions)
        shutil.copy2(validation_path, model_output / "validation_predictions.npy")
        scores = study.score_predictions(labels[:, :2], predictions,
                                         validation_labels[:, :2], validation_predictions,
                                         novel, int(metrics.get("seed", 0)))
        parity_scores(scores, labels, predictions, summaries)
        scores.update(model=model_name, seed=metrics.get("seed"), source_run=str(run),
                      parameters=metrics.get("parameters"),
                      training_seconds=metrics.get("training_seconds"),
                      best_epoch=metrics.get("best_epoch"),
                      evaluation="frozen MLP; prepared test split; no retraining",
                      feature_layout=protocol["triangle_layout"])
        (model_output / "metrics.json").write_text(json.dumps(scores, indent=2, allow_nan=False) + "\n")
        manifest["models"].append({"model": model_name, "seed": metrics.get("seed"),
                                   "run_name": run.name, "source_run": str(run),
                                   "output_run": run.name})
    if not manifest["models"]:
        raise ValueError("No model runs were evaluated")
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print(f"Evaluated {len(manifest['models'])} runs; results are in {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--reference", required=True, type=Path)
    prepare_parser.add_argument("--output", required=True, type=Path)
    prepare_parser.add_argument("--size", type=int, default=48)
    prepare_parser.add_argument("--hash-seed", type=int, default=17)
    prepare_parser.add_argument("--test", type=Path, help="Optional fresh raw dataset override")
    prepare_parser.add_argument("--reuse-training", type=Path,
                                help="Optional compatible triangle-prepared directory; reuse train/validation artifacts by hardlink")
    evaluate_parser = commands.add_parser("evaluate")
    evaluate_parser.add_argument("--prepared", required=True, type=Path)
    evaluate_parser.add_argument("--runs", required=True, nargs="+", type=Path)
    evaluate_parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            if args.size < 12 or (args.size - 4) % 4:
                raise ValueError("--size must be at least 12 and have the form 4*m+4 (two boundary anchors, two stash slots)")
            prepare(args)
        else:
            evaluate(args)
    except (OSError, ValueError, LayoutOverflow) as exc:
        parser.exit(2, f"affinity_triangle_study: {exc}\n")


if __name__ == "__main__":
    main()
