#!/usr/bin/env python3
"""Train/evaluate CPU MLPs for low-gap screening or gap regression.

Training reads only the designated train prefix and validation split. A test
split is never opened by the train command. Classification scores are
class-weighted BCE scores, not calibrated probabilities.
"""

import argparse
import copy
import json
import math
import os
from pathlib import Path
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from scipy.special import expit
from sklearn.metrics import average_precision_score, roc_auc_score
import torch
from torch import nn

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
import affinity_gap_study as study


FEATURES = ("syndrome", "summary", "packed_triangle", "hashed_triangle")
GAP_NAMES = ("complementary", "cluster")
RECALL_TARGETS = (0.95, 0.99, 0.999)
MATCH_FIELDS = ("distance", "rounds", "p", "num_detectors", "memory_basis",
                "include_opposite_basis_detectors")
STATIC_HASHES = ("circuit.stim", "raw.dem", "decomposed.dem", "graph_edges.csv",
                 "coords.csv", "affinity.f64", "hops.u8")


def find_protocol(path):
    path = Path(path).expanduser().resolve()
    for candidate in (path / "protocol.json", path.parent / "protocol.json"):
        if candidate.is_file():
            return json.loads(candidate.read_text()), candidate
    raise ValueError(f"No protocol.json found at or above {path}")


def split_file(path, split, name):
    """Resolve either a split directory or a flattened prepared directory."""
    path = Path(path).expanduser().resolve()
    candidates = (path / name, path / f"{split}_{name}", path / split / name)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise ValueError(f"Missing {split} {name} under {path}")


def load_memmap(path, expected_ndim=2):
    array = np.load(path, mmap_mode="r", allow_pickle=False)
    if not isinstance(array, np.ndarray) or array.ndim != expected_ndim:
        raise ValueError(f"Expected a {expected_ndim}D NumPy array at {path}")
    return array


def load_split(path, split, feature):
    labels = load_memmap(split_file(path, split, "labels.npy"))
    x = load_memmap(split_file(path, split, f"{feature}.npy"))
    if labels.shape[1] < 2 or x.shape[0] != labels.shape[0]:
        raise ValueError(f"Invalid {split} array shapes: labels={labels.shape}, features={x.shape}")
    gaps = np.asarray(labels[:, :2], dtype=np.float64)
    if not np.isfinite(gaps).all() or (gaps < 0).any():
        raise ValueError(f"{split} gap labels must be finite and nonnegative")
    return x, labels, gaps


def validate_protocols(train_protocol, validation_protocol):
    # prepare_gap_learning.py writes one independent exporter dataset per
    # protocol; dynamic labels/syndromes/seed hashes intentionally differ.
    train_meta = train_protocol.get("source_metadata")
    validation_meta = validation_protocol.get("source_metadata")
    if not isinstance(train_meta, dict) or not isinstance(validation_meta, dict):
        raise ValueError("Could not identify train and validation metadata in protocol files")
    if train_meta.get("seed") is None or validation_meta.get("seed") is None:
        raise ValueError("Train and validation source metadata must include seeds")
    if train_meta["seed"] == validation_meta["seed"]:
        raise ValueError("Training and validation seeds must differ")
    for key in MATCH_FIELDS:
        if train_meta.get(key) != validation_meta.get(key):
            raise ValueError(f"Train and validation dataset mismatch: {key}")
    train_hashes = train_protocol.get("source_sha256", {})
    validation_hashes = validation_protocol.get("source_sha256", {})
    for name in STATIC_HASHES:
        if train_hashes.get(name) is None or train_hashes.get(name) != validation_hashes.get(name):
            raise ValueError(f"Train and validation topology hash differs or is missing: {name}")
    for key in ("log_distance_scale", "triangle_layout"):
        if train_protocol.get(key) != validation_protocol.get(key):
            raise ValueError(f"Train and validation {key.replace('_', ' ')} differ")
    if train_protocol.get("features") != validation_protocol.get("features"):
        raise ValueError("Train and validation feature lists differ")
    return train_meta, validation_meta


def chunked_standardizer(x, count, floor, chunk_size=8192):
    """Two-pass float64 moments over a training prefix; never copies all X."""
    width = x.shape[1]
    total = np.zeros(width, dtype=np.float64)
    for start in range(0, count, chunk_size):
        block = np.asarray(x[start:min(start + chunk_size, count)], dtype=np.float64)
        if not np.isfinite(block).all():
            raise ValueError("Training features contain non-finite values")
        total += block.sum(axis=0, dtype=np.float64)
    mean = total / count
    squared = np.zeros(width, dtype=np.float64)
    for start in range(0, count, chunk_size):
        block = np.asarray(x[start:min(start + chunk_size, count)], dtype=np.float64)
        delta = block - mean
        squared += np.square(delta).sum(axis=0, dtype=np.float64)
    std = np.maximum(np.sqrt(squared / count), floor)
    return mean, std


def binary_metrics(labels, scores, threshold):
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=np.float64)
    predicted = scores >= threshold
    tp = int(np.count_nonzero(predicted & labels))
    fn = int(np.count_nonzero(~predicted & labels))
    fp = int(np.count_nonzero(predicted & ~labels))
    tn = int(np.count_nonzero(~predicted & ~labels))
    n = len(labels)
    nlow = tp + fn
    nhigh = fp + tn
    npredlow = tp + fp
    npredhigh = fn + tn
    recall = tp / nlow if nlow else None
    specificity = tn / nhigh if nhigh else None
    return {
        "threshold_on_low_score": float(threshold),
        "tp_low": tp, "fn_low_predicted_high": fn,
        "fp_high_predicted_low": fp, "tn_high": tn,
        "low_count": nlow, "high_count": nhigh,
        "predicted_low_expensive_fraction": npredlow / n if n else None,
        "predicted_high_accepted_fraction": npredhigh / n if n else None,
        "low_recall": recall,
        "low_precision": tp / npredlow if npredlow else None,
        "false_high_rate_among_low": fn / nlow if nlow else None,
        "false_high_fraction_all": fn / n if n else None,
        "false_high_fraction_among_predhigh": fn / npredhigh if npredhigh else None,
        "balanced_accuracy": ((recall + specificity) / 2 if recall is not None and specificity is not None
                              else None),
        "all_high_baseline_accuracy": nhigh / n if n else None,
        "all_high_baseline_low_recall": 0.0 if nlow else None,
    }


def ranking_metrics(labels, scores):
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=np.float64)
    positives = int(labels.sum())
    return {
        "low_count": positives,
        "high_count": int(len(labels) - positives),
        "low_prevalence": float(labels.mean()) if len(labels) else None,
        "roc_auc": (float(roc_auc_score(labels, scores))
                    if 0 < positives < len(labels) else None),
        "average_precision": (float(average_precision_score(labels, scores))
                              if positives else None),
    }


def low_scores(predictions, objective):
    # predict_model already applies sigmoid to classifier logits.
    return predictions if objective == "classification" else -predictions


def operating_thresholds(labels, predictions, objective, gap_threshold):
    output = []
    scores = low_scores(predictions, objective)
    low = labels <= gap_threshold
    for target in range(2):
        if objective == "classification":
            default = .5
        else:
            default = -gap_threshold
        item = {"target": GAP_NAMES[target], "default": {
            "kind": ("classification_score" if objective == "classification"
                     else "negative_predicted_gap_db"),
            "threshold": float(default),
            "metrics": binary_metrics(low[:, target], scores[:, target], default)}}
        for recall in RECALL_TARGETS:
            positive_scores = scores[low[:, target], target]
            if not len(positive_scores):
                threshold = None
                attained = None
                count = 0
            else:
                threshold = float(np.quantile(positive_scores, 1.0 - recall, method="lower"))
                attained = float(np.mean(positive_scores >= threshold))
                count = len(positive_scores)
            item[f"target_recall_{recall:g}"] = {
                "kind": "low_score" if objective == "classification" else "negative_predicted_gap_db",
                "threshold": threshold,
                "validation_low_count": count,
                "validation_recall_attained": attained,
                "metrics": (binary_metrics(low[:, target], scores[:, target], threshold)
                            if threshold is not None else None),
            }
        output.append(item)
    return output


def validation_metrics(labels, predictions, objective, threshold):
    result = {"objective": objective, "gap_threshold_db": threshold, "targets": {}}
    for j, name in enumerate(GAP_NAMES):
        low = labels[:, j] <= threshold
        scores = low_scores(predictions, objective)[:, j]
        target = {"ranking": ranking_metrics(low, scores)}
        if objective == "regression":
            target["regression"] = study.regression_metrics(labels[:, j], predictions[:, j])
        result["targets"][name] = target
    result["operating_points"] = operating_thresholds(labels, predictions, objective, threshold)
    return result


def standardize_batch(x, indices, mean32, std32):
    batch = np.array(x[indices], dtype=np.float32, copy=True, order="C")
    if not np.isfinite(batch).all():
        raise ValueError("Features contain non-finite values")
    batch -= mean32
    batch /= std32
    return batch


def predict_model(model, x, feature_mean, feature_std, objective, target_mean,
                  target_std, batch_size):
    predictions = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(x), batch_size):
            indices = np.arange(start, min(start + batch_size, len(x)))
            batch = torch.from_numpy(standardize_batch(x, indices, feature_mean, feature_std))
            raw = model(batch).cpu().numpy()
            if objective == "classification":
                batch_predictions = expit(raw).astype(np.float32)
            else:
                batch_predictions = np.maximum(raw * target_std + target_mean, 0.).astype(np.float32)
            if not np.isfinite(batch_predictions).all():
                raise ValueError("Model produced non-finite predictions")
            predictions.append(batch_predictions)
    if not predictions:
        return np.empty((0, 2), dtype=np.float32)
    return np.concatenate(predictions, axis=0)


def train(args):
    torch.set_num_threads(args.threads)
    torch.set_num_interop_threads(1)
    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)

    train_path = args.train.expanduser().resolve()
    validation_path = args.validation.expanduser().resolve()
    train_protocol, train_protocol_path = find_protocol(train_path)
    validation_protocol, validation_protocol_path = find_protocol(validation_path)
    train_meta, validation_meta = validate_protocols(train_protocol, validation_protocol)
    xtrain, train_labels, ytrain_all = load_split(train_path, "train", args.feature)
    xval, validation_labels, yval_all = load_split(validation_path, "validation", args.feature)
    if xtrain.shape[1] != xval.shape[1]:
        raise ValueError("Train and validation feature widths differ")
    ntrain = len(xtrain) if args.training_shots is None else args.training_shots
    if ntrain < 1 or ntrain > len(xtrain):
        raise ValueError(f"training-shots must be in [1, {len(xtrain)}]")
    if len(xval) == 0:
        raise ValueError("Validation split is empty")

    output = args.output.expanduser().resolve()
    if output.exists():
        raise ValueError(f"Output already exists: {output}")
    floor = .05 if args.feature.endswith("triangle") else 1e-6
    input_mean, input_std = chunked_standardizer(xtrain, ntrain, floor)
    if not np.isfinite(input_mean).all() or not np.isfinite(input_std).all():
        raise ValueError("Could not compute finite training feature statistics")
    input_mean32 = input_mean.astype(np.float32)
    input_std32 = input_std.astype(np.float32)

    if args.objective == "regression":
        ytrain = ytrain_all[:ntrain]
        target_mean = ytrain.mean(axis=0, dtype=np.float64)
        target_std = np.maximum(ytrain.std(axis=0, dtype=np.float64), 1.0)
        target_mean32, target_std32 = target_mean.astype(np.float32), target_std.astype(np.float32)
        train_targets = torch.from_numpy(((ytrain - target_mean) / target_std).astype(np.float32))
        val_targets = torch.from_numpy(((yval_all - target_mean) / target_std).astype(np.float32))
        pos_weight = None
        loss_fn = nn.MSELoss()
    else:
        low_train = ytrain_all[:ntrain] <= args.threshold_db
        positives = low_train.sum(axis=0, dtype=np.int64)
        negatives = ntrain - positives
        if (positives == 0).any() or (negatives == 0).any():
            raise ValueError("Classification needs both low- and high-gap examples for both targets")
        pos_weight_np = negatives.astype(np.float64) / positives.astype(np.float64)
        pos_weight = torch.from_numpy(pos_weight_np.astype(np.float32))
        train_targets = torch.from_numpy(low_train.astype(np.float32))
        val_targets = torch.from_numpy((yval_all <= args.threshold_db).astype(np.float32))
        target_mean = target_std = None
        target_mean32 = target_std32 = None
        loss_fn = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    model = study.GapMLP(xtrain.shape[1])
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    best_loss, best_epoch, best_state = math.inf, 0, None
    history = []
    started = time.perf_counter()
    for epoch in range(args.epochs):
        model.train()
        order = rng.permutation(ntrain)
        total_loss = 0.0
        for start in range(0, ntrain, args.batch_size):
            indices = order[start:start + args.batch_size]
            xb = torch.from_numpy(standardize_batch(xtrain, indices, input_mean32, input_std32))
            optimizer.zero_grad(set_to_none=True)
            loss = loss_fn(model(xb), train_targets[indices])
            loss.backward()
            optimizer.step()
            total_loss += float(loss.detach()) * len(indices)

        val_parts = []
        model.eval()
        with torch.no_grad():
            for start in range(0, len(xval), args.batch_size):
                indices = np.arange(start, min(start + args.batch_size, len(xval)))
                xb = torch.from_numpy(standardize_batch(xval, indices, input_mean32, input_std32))
                val_parts.append(model(xb).cpu())
        raw_val = torch.cat(val_parts, dim=0)
        val_loss = float(loss_fn(raw_val, val_targets).item())
        row = {"epoch": epoch + 1, "train_loss": total_loss / ntrain,
               "validation_loss": val_loss}
        history.append(row)
        print(json.dumps(row), flush=True)
        if val_loss < best_loss - 1e-8:
            best_loss, best_epoch = val_loss, epoch + 1
            best_state = copy.deepcopy(model.state_dict())
        elif epoch + 1 - best_epoch >= args.patience:
            break

    if best_state is None:
        raise RuntimeError("Training did not produce a best checkpoint")
    model.load_state_dict(best_state)
    validation_predictions = predict_model(
        model, xval, input_mean32, input_std32, args.objective,
        target_mean32, target_std32, args.batch_size)
    val_metrics = validation_metrics(yval_all, validation_predictions,
                                     args.objective, args.threshold_db)
    output.mkdir(parents=True, exist_ok=False)
    np.save(output / "validation_predictions.npy", validation_predictions)
    with (output / "history.json").open("x") as stream:
        json.dump(history, stream, indent=2, allow_nan=False)
    protocol = train_protocol
    checkpoint = {
        "state_dict": model.state_dict(), "model": "gap_mlp",
        "feature": args.feature, "objective": args.objective,
        "input_mean": input_mean.tolist(), "input_std": input_std.tolist(),
        "target_mean": None if target_mean is None else target_mean.tolist(),
        "target_std": None if target_std is None else target_std.tolist(),
        "pos_weight": None if pos_weight is None else pos_weight.tolist(),
        "threshold_db": args.threshold_db,
        "operating_points": val_metrics["operating_points"],
        "protocol": protocol,
        "feature_shape": [int(xtrain.shape[1])],
        "train_metadata": train_meta, "validation_metadata": validation_meta,
        "configuration": {k: str(v) if isinstance(v, Path) else v
                           for k, v in vars(args).items()},
        "training_prefix_rows": ntrain,
        "training_prefix_rule": "first N rows; only this prefix used for scaler and objective statistics",
        "validation_protocol_path": str(validation_protocol_path),
    }
    torch.save(checkpoint, output / "model.pt")
    metrics = {
        "model": f"gap_mlp_{args.objective}_{args.feature}",
        "objective": args.objective, "feature": args.feature,
        "seed": args.seed, "threads": args.threads,
        "training_prefix_rows": ntrain, "training_rows_available": len(xtrain),
        "validation_rows": len(xval), "training_seconds": time.perf_counter() - started,
        "best_epoch": best_epoch, "best_validation_loss": best_loss,
        "parameters": sum(parameter.numel() for parameter in model.parameters()),
        "threshold_db": args.threshold_db,
        "classification_pos_weight": None if pos_weight is None else pos_weight.tolist(),
        "validation": val_metrics,
        "feature_protocol_path": str(train_protocol_path),
        "train_protocol_path": str(train_protocol_path),
        "validation_protocol_path": str(validation_protocol_path),
        "feature_shape": [int(xtrain.shape[1])],
    }
    with (output / "metrics.json").open("x") as stream:
        json.dump(metrics, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(metrics, indent=2), flush=True)


def compare_protocol(checkpoint, data_protocol, feature, test_metadata):
    source_protocol = checkpoint.get("protocol")
    if not isinstance(source_protocol, dict):
        raise ValueError("Checkpoint lacks its source protocol")
    source_meta = checkpoint.get("train_metadata")
    validation_meta = checkpoint.get("validation_metadata")
    if not isinstance(source_meta, dict) or not isinstance(validation_meta, dict):
        raise ValueError("Checkpoint lacks train/validation source metadata")
    if not isinstance(test_metadata, dict):
        raise ValueError("Evaluation protocol lacks source metadata")
    for key in MATCH_FIELDS:
        if source_meta.get(key) != test_metadata.get(key):
            raise ValueError(f"Evaluation dataset mismatch for {key}")
    source_hashes = source_protocol.get("source_sha256", {})
    data_hashes = data_protocol.get("source_sha256", {})
    for name in STATIC_HASHES:
        if source_hashes.get(name) is None or source_hashes.get(name) != data_hashes.get(name):
            raise ValueError(f"Evaluation topology hash differs or is missing: {name}")
    if source_protocol.get("log_distance_scale") != data_protocol.get("log_distance_scale"):
        raise ValueError("Evaluation log-distance scale differs from training checkpoint")
    if test_metadata.get("seed") in {source_meta.get("seed"), validation_meta.get("seed")}:
        raise ValueError("Evaluation test seed overlaps checkpoint train/validation seed")
    if feature not in data_protocol.get("features", []):
        raise ValueError(f"Evaluation dataset does not contain feature {feature}")
    if feature.endswith("triangle"):
        source_layout = source_protocol.get("triangle_layout")
        data_layout = data_protocol.get("triangle_layout")
        if not isinstance(source_layout, dict) or not isinstance(data_layout, dict):
            raise ValueError("Triangle model/evaluation protocol lacks layout metadata")
        if source_layout != data_layout:
            raise ValueError("Evaluation triangle layout differs from training checkpoint")


def evaluate(args):
    data = args.data.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if output.exists():
        raise ValueError(f"Output already exists: {output}")
    data_protocol, _ = find_protocol(data)
    test_meta = data_protocol.get("source_metadata")
    if not isinstance(test_meta, dict) or test_meta.get("seed") is None:
        raise ValueError("Evaluation protocol must include source metadata and seed")
    labels = load_memmap(split_file(data, "test", "labels.npy"))
    ytest = np.asarray(labels[:, :2], dtype=np.float64)
    if not np.isfinite(ytest).all() or (ytest < 0).any():
        raise ValueError("Test gap labels must be finite and nonnegative")
    import torch
    torch.set_num_threads(args.threads)
    torch.set_num_interop_threads(1)
    output.mkdir(parents=True, exist_ok=False)
    manifest = {"data": str(data), "test_metadata": test_meta, "models": []}
    for run in args.runs:
        run = run.expanduser().resolve()
        checkpoint_path = run / "model.pt"
        if not checkpoint_path.is_file():
            raise ValueError(f"Missing checkpoint: {checkpoint_path}")
        checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
        feature = checkpoint.get("feature")
        objective = checkpoint.get("objective")
        if feature not in FEATURES or objective not in ("regression", "classification"):
            raise ValueError(f"Unsupported checkpoint feature/objective: {feature}/{objective}")
        compare_protocol(checkpoint, data_protocol, feature, test_meta)
        x = load_memmap(split_file(data, "test", f"{feature}.npy"))
        if x.shape[1] != int(checkpoint.get("feature_shape", [-1])[0]):
            raise ValueError(f"Test feature width mismatch for {feature}")
        model = study.GapMLP(x.shape[1])
        model.load_state_dict(checkpoint["state_dict"])
        predictions = predict_model(
            model, x, np.asarray(checkpoint["input_mean"], dtype=np.float32),
            np.asarray(checkpoint["input_std"], dtype=np.float32), objective,
            None if objective == "classification" else np.asarray(checkpoint["target_mean"], dtype=np.float32),
            None if objective == "classification" else np.asarray(checkpoint["target_std"], dtype=np.float32),
            args.batch_size)
        low_labels = ytest <= float(checkpoint["threshold_db"])
        scores = low_scores(predictions, objective)
        per_target = {}
        saved_operating = checkpoint.get("operating_points", [])
        for j, name in enumerate(GAP_NAMES):
            item = {"ranking": ranking_metrics(low_labels[:, j], scores[:, j])}
            if objective == "regression":
                item["regression"] = study.regression_metrics(ytest[:, j], predictions[:, j])
            item["operating_points"] = {}
            source_target = saved_operating[j] if len(saved_operating) > j else {}
            for key, point in source_target.items():
                if not isinstance(point, dict):
                    continue
                cutoff = point.get("threshold")
                item["operating_points"][key] = {
                    "validation_threshold": cutoff,
                    "validation_low_count": point.get("validation_low_count"),
                    "validation_recall_attained": point.get("validation_recall_attained"),
                    "test": (binary_metrics(low_labels[:, j], scores[:, j], cutoff)
                             if cutoff is not None else None),
                }
            per_target[name] = item
        scores_json = {
            "model": run.name, "feature": feature, "objective": objective,
            "threshold_db": checkpoint["threshold_db"],
            "test_rows": len(ytest), "targets": per_target,
            "test_prediction_semantics": ("class-weighted BCE low-gap score; not a calibrated probability"
                                          if objective == "classification" else "nonnegative predicted gap in dB"),
            "operating_thresholds_source": "validation only; frozen in checkpoint",
        }
        model_output = output / run.name
        model_output.mkdir()
        np.save(model_output / "test_predictions.npy", predictions)
        with (model_output / "metrics.json").open("x") as stream:
            json.dump(scores_json, stream, indent=2, allow_nan=False)
            stream.write("\n")
        manifest["models"].append({"run": str(run), "output": run.name,
                                   "feature": feature, "objective": objective})
        print(f"evaluated {run.name}: {feature}/{objective}, {len(ytest)} shots", flush=True)
    if not manifest["models"]:
        raise ValueError("No model runs were evaluated")
    with (output / "run_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, allow_nan=False)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    train_parser = subparsers.add_parser("train")
    train_parser.add_argument("--train", required=True, type=Path,
                              help="train split directory or flattened prepared directory")
    train_parser.add_argument("--validation", required=True, type=Path,
                              help="validation split directory or flattened prepared directory")
    train_parser.add_argument("--output", required=True, type=Path)
    train_parser.add_argument("--feature", required=True, choices=FEATURES)
    train_parser.add_argument("--objective", required=True, choices=("regression", "classification"))
    train_parser.add_argument("--training-shots", type=int,
                              help="use the first N training rows; default uses all")
    train_parser.add_argument("--seed", type=int, default=101)
    train_parser.add_argument("--threads", type=int, default=4)
    train_parser.add_argument("--epochs", type=int, default=20)
    train_parser.add_argument("--patience", type=int, default=5)
    train_parser.add_argument("--threshold-db", type=float, default=40.)
    train_parser.add_argument("--batch-size", type=int, default=1024)
    train_parser.add_argument("--lr", type=float, default=.001)

    eval_parser = subparsers.add_parser("evaluate")
    eval_parser.add_argument("--data", required=True, type=Path)
    eval_parser.add_argument("--runs", required=True, nargs="+", type=Path)
    eval_parser.add_argument("--output", required=True, type=Path)
    eval_parser.add_argument("--threads", type=int, default=4)
    eval_parser.add_argument("--batch-size", type=int, default=4096)
    args = parser.parse_args()
    try:
        if args.command == "train":
            if args.training_shots is not None and args.training_shots < 1:
                raise ValueError("training-shots must be positive")
            if min(args.threads, args.epochs, args.patience, args.batch_size) < 1:
                raise ValueError("threads, epochs, patience, and batch size must be positive")
            if not math.isfinite(args.threshold_db) or args.threshold_db < 0:
                raise ValueError("threshold-db must be finite and nonnegative")
            if not math.isfinite(args.lr) or args.lr <= 0:
                raise ValueError("lr must be finite and positive")
            train(args)
        else:
            if min(args.threads, args.batch_size) < 1:
                raise ValueError("threads and batch size must be positive")
            evaluate(args)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"train_gap_confidence: {exc}\n")


if __name__ == "__main__":
    main()
