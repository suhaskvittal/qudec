#!/usr/bin/env python3
"""Evaluate locked affinity-gap models on one fresh exported holdout dataset."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
import affinity_gap_study as study


TOPOLOGY_FILES = ("circuit.stim", "raw.dem", "decomposed.dem", "graph_edges.csv", "coords.csv")
MATCH_FIELDS = ("distance", "rounds", "p", "num_detectors", "memory_basis",
                "include_opposite_basis_detectors")


def load_protocol(reference):
    path = reference / "protocol.json"
    if not path.is_file():
        raise ValueError(f"Reference protocol not found: {path}")
    protocol = json.loads(path.read_text())
    if len(protocol.get("datasets", [])) != 3 or len(protocol.get("metadata", [])) != 3:
        raise ValueError("Reference protocol must describe train, validation, and original test splits")
    if not protocol.get("topology_sha256"):
        raise ValueError("Reference protocol lacks topology_sha256; prepare data with the current study script")
    return protocol


def verify_dataset(dataset, protocol):
    if not (dataset / "metadata.json").is_file():
        raise ValueError("Holdout metadata.json is missing; it marks completed dataset export")
    metadata = json.loads((dataset / "metadata.json").read_text())
    expected = protocol["metadata"][0]
    mismatches = [key for key in MATCH_FIELDS
                  if key not in metadata or key not in expected or metadata[key] != expected[key]]
    if mismatches:
        raise ValueError("Holdout differs from reference in: " + ", ".join(mismatches))
    old_seeds = {meta["seed"] for meta in protocol["metadata"]}
    if metadata.get("seed") in old_seeds:
        raise ValueError(f"Holdout seed {metadata.get('seed')} repeats an existing split seed")
    for name in TOPOLOGY_FILES:
        file_path = dataset / name
        if not file_path.is_file():
            raise ValueError(f"Holdout topology file missing: {file_path}")
        digest = hashlib.sha256(file_path.read_bytes()).hexdigest()
        if digest != protocol["topology_sha256"].get(name):
            raise ValueError(f"Holdout topology hash mismatch: {name}")
    return metadata


def verify_holdout_lock(reference, dataset_metadata, source_runs):
    lock_path = reference.parent / "holdout_lock.json"
    if not lock_path.is_file():
        return None
    lock = json.loads(lock_path.read_text())
    if dataset_metadata.get("seed") != lock.get("future_holdout_seed"):
        raise ValueError("Holdout seed does not match holdout_lock.json")
    resolved_runs = [path.expanduser().resolve() for path in source_runs]
    expected_names = set(lock.get("runs", []))
    if {path.name for path in resolved_runs} != expected_names:
        raise ValueError("Source runs do not match the model set fixed in holdout_lock.json")
    repo = SCRIPT_DIR.parent
    expected_hashes = lock.get("sha256", {})
    for source_run in resolved_runs:
        metrics_path = source_run / "metrics.json"
        if not metrics_path.is_file():
            raise ValueError(f"Missing source metrics: {metrics_path}")
        model_name = json.loads(metrics_path.read_text()).get("model", source_run.name)
        artifact_name = None if model_name == "constant" else (
            "model.pt" if model_name.startswith(("cnn_", "mlp_")) else "model.joblib")
        if artifact_name is None:
            continue
        artifact = source_run / artifact_name
        if not artifact.is_file():
            raise ValueError(f"Missing locked model artifact: {artifact}")
        key = f"{source_run.name}/{artifact_name}"
        matches = [(name, digest) for name, digest in expected_hashes.items()
                   if Path(name).name == artifact_name and Path(name).parent.name == source_run.name]
        if len(matches) != 1:
            raise ValueError(f"No unique locked hash for {key}")
        expected = matches[0][1]
        actual = hashlib.sha256(artifact.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Locked model hash mismatch: {key}")
    for name in ("scripts/affinity_gap_study.py", "main/export_affinity_dataset.cpp"):
        expected = expected_hashes.get(name)
        path = repo / name
        if expected and (not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected):
            lock["post_lock_source_hash_mismatches"] = lock.get("post_lock_source_hash_mismatches", []) + [name]
    return lock


def write_holdout_prepared(dataset, reference, protocol, output, reuse_prepared=None):
    meta, syndromes, labels, affinity, hops = study.read_dataset(dataset)
    if labels.shape[1] != 6:
        raise ValueError("Expected six-column exported labels")
    n = meta["num_detectors"]
    old_train = protocol["datasets"][0]
    boundary, relevant = study.graph_geometry(dataset, n)
    scale = float(protocol["log_distance_scale"])
    size = int(protocol["matrix_size"])

    prepared = output / "prepared"
    prepared.mkdir()
    if reuse_prepared is not None:
        reuse_prepared = reuse_prepared.expanduser().resolve()
        names = ("test_syndrome.npy", "test_labels.npy", "test_summary.npy", "test_hop.npy",
                 "test_weighted.npy", "test_boundary.npy", "test_context.npy", "validation_labels.npy")
        for name in names:
            source = reuse_prepared / name
            if not source.is_file():
                raise ValueError(f"Missing reusable prepared feature: {source}")
            os.link(source, prepared / name)
        reused_labels = np.load(prepared / "test_labels.npy", mmap_mode="r")
        reused_syndromes = np.load(prepared / "test_syndrome.npy", mmap_mode="r")
        if not np.array_equal(reused_labels, labels) or not np.array_equal(reused_syndromes, syndromes):
            raise ValueError("Reusable prepared arrays do not correspond to the holdout dataset")
        if not np.array_equal(np.load(prepared / "validation_labels.npy"),
                              np.load(reference / "validation_labels.npy")):
            raise ValueError("Reusable validation labels do not match the reference split")
        summaries = np.load(prepared / "test_summary.npy", mmap_mode="r")
    else:
        np.save(prepared / "test_syndrome.npy", syndromes.astype(np.float32))
        np.save(prepared / "test_labels.npy", labels)
        np.save(prepared / "validation_labels.npy", np.load(reference / "validation_labels.npy"))
        images = {kind: np.lib.format.open_memmap(
            prepared / f"test_{kind}.npy", mode="w+", dtype="float32",
            shape=(len(syndromes), 7 if kind == "context" else 3, size, size))
            for kind in ("hop", "weighted", "boundary", "context")}
        summaries = np.empty((len(syndromes), len(study.SUMMARY_NAMES)), dtype=np.float32)
        for index, syndrome in enumerate(syndromes):
            ids = np.flatnonzero(syndrome)
            summaries[index] = study.summary_features(ids, affinity, boundary, relevant, scale)
            images["hop"][index] = study.image_features(
                ids, hops, size, scale=meta["distance"], hop=True)
            images["weighted"][index] = study.image_features(ids, affinity, size, scale=scale)
            images["boundary"][index] = study.image_features(
                ids, affinity, size, boundary, scale=scale)
            images["context"][index] = study.context_features(ids, affinity, size, boundary, scale)
        for image in images.values():
            image.flush()
        np.save(prepared / "test_summary.npy", summaries)

    train_syndromes = np.load(reference / "train_syndrome.npy", mmap_mode="r")
    train_hashes = {hashlib.sha256(np.packbits(row.astype(np.uint8)).tobytes()).digest()
                    for row in train_syndromes}
    novel = np.asarray([hashlib.sha256(np.packbits(row).tobytes()).digest() not in train_hashes
                        for row in syndromes], dtype=bool)
    np.save(prepared / "test_novel.npy", novel)

    updated = json.loads(json.dumps(protocol))
    updated["datasets"][-1] = str(dataset.resolve())
    updated["metadata"][-1] = meta
    updated["test_novel_shots"] = int(novel.sum())
    updated["final_holdout"] = {
        "dataset": str(dataset.resolve()), "seed": meta["seed"],
        "note": "Fresh final holdout; source models and hyperparameters remain locked; no retraining.",
        "novel_syndromes_vs_original_train": int(novel.sum())}
    (prepared / "protocol.json").write_text(json.dumps(updated, indent=2, allow_nan=False) + "\n")
    return prepared, labels, summaries, novel, updated


def predict_torch(checkpoint_path, feature_arrays, model_name, batch_size=256):
    import torch

    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    saved_model = checkpoint.get("model")
    if saved_model is None and checkpoint.get("feature") in ("hop", "weighted", "boundary", "context"):
        # Earlier CNN checkpoints predate the explicit model-name field.
        saved_model = "cnn_" + checkpoint["feature"]
    if saved_model != model_name:
        raise ValueError(f"Checkpoint model {saved_model!r} does not match metrics model {model_name!r}")
    feature = checkpoint["feature"]
    array = feature_arrays[feature]
    target_mean = np.asarray(checkpoint["mean"], dtype=np.float64)
    target_std = np.asarray(checkpoint["std"], dtype=np.float64)
    if model_name.startswith("mlp_"):
        input_mean = np.asarray(checkpoint["input_mean"], dtype=np.float32)
        input_std = np.asarray(checkpoint["input_std"], dtype=np.float32)
        model = study.GapMLP(array.shape[1])
    else:
        model = study.GapCNN(channels=array.shape[1])
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    predictions = []
    with torch.no_grad():
        for start in range(0, len(array), batch_size):
            x = np.array(array[start:start + batch_size], dtype=np.float32, copy=True, order="C")
            if model_name.startswith("mlp_"):
                x = (x - input_mean) / input_std
            y = model(torch.from_numpy(x)).numpy()
            predictions.append(y * target_std + target_mean)
    if not predictions:
        raise ValueError("Cannot evaluate an empty holdout dataset")
    return np.maximum(np.concatenate(predictions), 0.)


def predict_sklearn(source_run, model_name, summary, syndrome, train_labels, batch_size=512):
    if model_name == "constant":
        return np.tile(train_labels[:, :2].mean(axis=0), (len(summary), 1))
    import joblib
    model = joblib.load(source_run / "model.joblib")
    if model_name == "weight":
        features = summary[:, :1]
    elif model_name in ("ridge_affinity", "trees_affinity"):
        features = summary[:, :15]
    elif model_name == "trees_basis":
        features = summary[:, :17]
    elif model_name == "trees_boundary":
        features = summary
    elif model_name == "trees_syndrome":
        features = syndrome
    else:
        raise ValueError(f"No inference path for model {model_name!r}")
    batches = [model.predict(features[i:i + batch_size])
               for i in range(0, len(features), batch_size)]
    if not batches:
        raise ValueError("Cannot evaluate an empty holdout dataset")
    return np.maximum(np.concatenate(batches), 0.)


def parity_scores(scores, labels, predictions, summary):
    weights = summary[:, study.SUMMARY_NAMES.index("weight")].astype(np.int64)
    for target_index, target in enumerate(("complementary", "cluster")):
        for parity, label in ((0, "even_hamming_weight"), (1, "odd_hamming_weight")):
            mask = weights % 2 == parity
            if mask.any():
                scores[target][label] = {
                    "shots": int(mask.sum()),
                    **study.regression_metrics(labels[mask, target_index], predictions[mask, target_index])}
            else:
                scores[target][label] = {"shots": 0, "mae_db": None, "rmse_db": None,
                                         "r2": None, "spearman": None}


def evaluate_runs(source_runs, output, prepared, labels, summaries, novel,
                  reference, protocol):
    import torch

    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    feature_arrays = {kind: np.load(prepared / f"test_{kind}.npy", mmap_mode="r")
                      for kind in ("hop", "weighted", "boundary", "context")}
    feature_arrays["summary"] = summaries
    syndrome = np.load(prepared / "test_syndrome.npy", mmap_mode="r")
    train_labels = np.load(reference / "train_labels.npy")
    validation_labels = np.load(reference / "validation_labels.npy")
    manifest = {"holdout_dataset": str(protocol["datasets"][-1]), "models": [], "excluded_runs": []}
    run_names = set()
    for source_run in source_runs:
        source_run = source_run.expanduser().resolve()
        metrics_path = source_run / "metrics.json"
        if not metrics_path.is_file():
            raise ValueError(f"Missing source metrics: {metrics_path}")
        source_metrics = json.loads(metrics_path.read_text())
        model_name = source_metrics.get("model", source_run.name)
        if "quantile" in model_name.lower():
            manifest["excluded_runs"].append({"run": str(source_run), "reason": "cancelled quantile variant"})
            continue
        if source_run.name in run_names or source_run.name == "prepared":
            raise ValueError(f"Source run directory names must be unique: {source_run.name}")
        run_names.add(source_run.name)
        validation_prediction_path = source_run / "validation_predictions.npy"
        if not validation_prediction_path.is_file():
            raise ValueError(f"Missing source validation predictions: {validation_prediction_path}")
        validation_predictions = np.load(validation_prediction_path)
        if validation_predictions.shape != (len(validation_labels), 2):
            raise ValueError(f"Validation prediction shape mismatch for {source_run}")
        if model_name.startswith("cnn_") or model_name.startswith("mlp_"):
            checkpoint_path = source_run / "model.pt"
            if not checkpoint_path.is_file():
                raise ValueError(f"Missing PyTorch checkpoint: {checkpoint_path}")
            predictions = predict_torch(checkpoint_path, feature_arrays, model_name)
        else:
            predictions = predict_sklearn(source_run, model_name, summaries, syndrome, train_labels)
        if predictions.shape != (len(labels), 2) or not np.isfinite(predictions).all():
            raise ValueError(f"Invalid holdout predictions for {source_run}")
        model_output = output / source_run.name
        model_output.mkdir()
        np.save(model_output / "test_predictions.npy", predictions)
        shutil.copy2(validation_prediction_path, model_output / "validation_predictions.npy")
        scores = study.score_predictions(labels[:, :2], predictions, validation_labels[:, :2],
                                         validation_predictions, novel,
                                         int(source_metrics.get("seed", 0)))
        parity_scores(scores, labels, predictions, summaries)
        scores.update(model=model_name, seed=source_metrics.get("seed"),
                      source_run=str(source_run), parameters=source_metrics.get("parameters"),
                      training_seconds=source_metrics.get("training_seconds"),
                      best_epoch=source_metrics.get("best_epoch"),
                      holdout_dataset=str(protocol["datasets"][-1]),
                      evaluation="locked model; fresh final holdout; no retraining")
        (model_output / "metrics.json").write_text(json.dumps(scores, indent=2, allow_nan=False) + "\n")
        manifest["models"].append({"model": model_name, "seed": source_metrics.get("seed"),
                                   "source_run": str(source_run), "output_run": model_output.name,
                                   "parameters": source_metrics.get("parameters")})
    if not manifest["models"]:
        raise ValueError("No eligible source models were evaluated")
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path,
                        help="Fresh exported dataset directory")
    parser.add_argument("--reference", required=True, type=Path,
                        help="Original prepared_d7 directory")
    parser.add_argument("--runs", required=True, nargs="+", type=Path,
                        help="Locked trained model run directories")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reuse-prepared", type=Path,
                        help="Reuse already-generated holdout prepared arrays by hardlink")
    args = parser.parse_args()

    dataset = args.dataset.expanduser().resolve()
    reference = args.reference.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if output.exists():
        parser.error(f"output path already exists: {output}")
    protocol = load_protocol(reference)
    metadata = verify_dataset(dataset, protocol)
    lock = verify_holdout_lock(reference, metadata, args.runs)
    output.mkdir(parents=True, exist_ok=False)
    prepared, labels, summaries, novel, updated_protocol = write_holdout_prepared(
        dataset, reference, protocol, output, args.reuse_prepared)
    evaluate_runs(args.runs, output, prepared, labels, summaries, novel, reference, updated_protocol)
    if lock is not None:
        (output / "run_manifest.json").write_text(json.dumps({
            **json.loads((output / "run_manifest.json").read_text()),
            "holdout_lock": str((reference.parent / "holdout_lock.json").resolve()),
            "holdout_lock_sha256": hashlib.sha256(
                (reference.parent / "holdout_lock.json").read_bytes()).hexdigest(),
            "post_lock_source_hash_mismatches": lock.get("post_lock_source_hash_mismatches", [])},
            indent=2, allow_nan=False) + "\n")
    print(f"Evaluated locked models on fresh holdout; results are in {output}")


if __name__ == "__main__":
    main()
