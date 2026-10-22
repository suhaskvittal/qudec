#!/usr/bin/env python3
"""CPU study of predecode affinity features and unsigned surface-code gap labels.

Contributor: OpenAI Codex. See docs/affinity-gap-study.md for the protocol.
All generated files belong under .agentic_artifacts/ by default.
"""

import os
_ARTIFACT_ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             ".agentic_artifacts", "affinity_gap")
os.makedirs(os.path.join(_ARTIFACT_ROOT, "tmp"), exist_ok=True)
os.environ.setdefault("TMPDIR", os.path.join(_ARTIFACT_ROOT, "tmp"))
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits
import torch
from torch import nn


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def read_dataset(path):
    meta = json.loads((path / "metadata.json").read_text())
    n, shots = meta["num_detectors"], meta["shots"]
    arrays = []
    for name, dtype, shape in [("syndromes.u8", "u1", (shots, n)),
                               ("labels.f64", "<f8", (shots, 6)),
                               ("affinity.f64", "<f8", (n + 1, n + 1)),
                               ("hops.u8", "u1", (n + 1, n + 1))]:
        if (path / name).stat().st_size != np.prod(shape) * np.dtype(dtype).itemsize:
            raise ValueError(f"Truncated or incompatible {path / name}")
        arrays.append(np.fromfile(path / name, dtype=dtype).reshape(shape))
    if not np.isfinite(arrays[1]).all() or (arrays[1][:, :2] < 0).any():
        raise ValueError("Invalid unsigned gap labels")
    return meta, *arrays


def graph_geometry(path, n):
    """Precompute distances to distinct boundaries without allowing boundary transit.

    Exported graph components are not yet merged. Independent parallel component
    probabilities combine by XOR, matching PyMatching's independent graph model.
    No shot, correction, growth radius, or decoded endpoint enters this function.
    """
    probabilities = {}
    with (path / "graph_edges.csv").open() as stream:
        for row in csv.DictReader(stream):
            u, v = int(row["detector_u"]), int(row["detector_v"])
            if v == -1:
                v = n + int(row["boundary_class"])
            key = tuple(sorted((u, v)))
            p = 1 / (1 + np.exp(float(row["weight_natural_log"])))
            old = probabilities.get(key, 0.)
            probabilities[key] = old + p - 2 * old * p
    rows, cols, weights = [], [], []
    # Boundary vertices have outgoing edges only. A path can start there but
    # cannot leave a physical component via one boundary and reenter elsewhere.
    for (u, v), p in probabilities.items():
        w = np.log((1 - p) / p)
        rows.append(v); cols.append(u); weights.append(w)
        if v < n:
            rows.append(u); cols.append(v); weights.append(w)
    graph = csr_matrix((weights, (rows, cols)), shape=(n + 2, n + 2))
    boundary = dijkstra(graph, directed=True, indices=[n, n + 1])[:, :n]
    relevant = np.isfinite(boundary[1])
    return boundary, relevant


SUMMARY_NAMES = ["weight", "weight_parity", "mean_nearest_log_distance",
                 "max_nearest_log_distance", "mean_nearest_second_margin",
                 "min_nearest_second_margin", "mean_effective_partners",
                 "max_effective_partners", "isolated_rows", "mean_pair_log_distance",
                 "min_pair_log_distance", "max_pair_log_distance",
                 "mean_shared_boundary_log_distance", "min_shared_boundary_log_distance",
                 "max_shared_boundary_log_distance", "memory_weight", "opposite_weight",
                 "min_boundary0_distance", "min_boundary1_distance",
                 "mean_boundary0_distance", "mean_boundary1_distance",
                 "min_boundary_asymmetry", "mean_boundary_asymmetry",
                 "max_boundary_asymmetry"]


def summary_features(ids, affinity, boundary, relevant, scale):
    k = len(ids)
    if not k:
        return np.zeros(len(SUMMARY_NAMES), dtype=np.float32)
    raw = affinity[np.ix_(ids, ids)].copy()
    distances = np.where(raw > 0, -np.log(np.maximum(raw, 1e-300)) / scale, 2.)
    np.fill_diagonal(distances, 2.)
    nearest = np.sort(distances, axis=1)[:, :2]
    # For a singleton, use two equal sentinels rather than a missing column.
    if k == 1:
        nearest = np.full((1, 2), 2.)
    margin = nearest[:, 1] - nearest[:, 0]
    sums = raw.sum(axis=1)
    eff = np.divide(sums ** 2, (raw ** 2).sum(axis=1),
                    out=np.zeros(k), where=(raw ** 2).sum(axis=1) > 0)
    pairs = distances[np.triu_indices(k, 1)]
    pairs = pairs[pairs < 2.]
    shared = affinity[ids, -1]
    shared = -np.log(np.maximum(shared, 1e-300)) / scale
    base = [k, k % 2, nearest[:, 0].mean(), nearest[:, 0].max(), margin.mean(), margin.min(),
            eff.mean(), eff.max(), np.count_nonzero(sums == 0),
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


def image_features(ids, table, size, boundary=None, scale=1., hop=False):
    """Pack ascending active IDs after persistent boundary anchor(s), with masks.

    Channel 0 is normalized path strength (one minus clipped log distance).
    Channel 1 distinguishes connected, active pairs from padding/disconnection.
    Channel 2 marks active node diagonals and the identities of boundary anchors.
    """
    anchors = 2 if boundary is not None else 1
    k = len(ids)
    if k + anchors > size:
        raise ValueError(f"{k} active detectors exceed capacity {size - anchors}; increase --size")
    values = table[np.ix_(ids, ids)]
    if hop:
        connected = values < 31
        strength = np.clip(1 - values / scale, 0., 1.)
    else:
        connected = values > 0
        strength = np.clip(1 + np.log(np.maximum(values, 1e-300)) / scale, 0., 1.)
    out = np.zeros((3, size, size), dtype=np.float32)
    sl = slice(anchors, anchors + k)
    out[0, sl, sl] = strength * connected
    out[1, sl, sl] = connected
    # Explicit self presence includes singletons even though affinity diagonal=0.
    active = np.arange(anchors, anchors + k)
    out[0, active, active] = 1.
    out[1, active, active] = 1.
    out[2, active, active] = 1.
    for b in range(anchors):
        out[0, b, b] = out[1, b, b] = 1.
        out[2, b, b] = -1. if b == 0 else -.5
        if boundary is None:
            vals = table[ids, -1]
            good = vals < 31 if hop else vals > 0
            strength_b = (1 - vals / scale) if hop else (1 + np.log(np.maximum(vals, 1e-300)) / scale)
        else:
            vals = boundary[b, ids]
            good = np.isfinite(vals)
            strength_b = 1 - vals / scale
        strength_b = np.clip(strength_b, 0., 1.) * good
        out[0, b, sl] = out[0, sl, b] = strength_b
        out[1, b, sl] = out[1, sl, b] = good
    return out


def context_features(ids, table, size, boundary, scale):
    """Expose both endpoints' boundary costs at every active pair pixel.

    Convolutions on event-rank axes have no fixed geometric neighborhoods.
    Broadcasting costs avoids forcing a local kernel to propagate information
    from the two anchor rows across the whole matrix.
    """
    base = image_features(ids, table, size, boundary, scale)
    out = np.zeros((7, size, size), dtype=np.float32)
    out[:3] = base
    k = len(ids) + 2
    fields = np.zeros((2, size), dtype=np.float32)
    fields[:, :2] = np.eye(2)
    fields[:, 2:k] = np.clip(1 - boundary[:, ids] / scale, 0., 1.)
    out[3:5, :k, :k] = fields[:, :k, None]
    out[5:7, :k, :k] = fields[:, None, :k]
    return out


def prepare(args):
    args.output.mkdir(parents=True, exist_ok=False)
    datasets = [read_dataset(p) for p in [args.train, args.validation, args.test]]
    metas = [d[0] for d in datasets]
    for key in ["distance", "rounds", "p", "num_detectors", "memory_basis", "include_opposite_basis_detectors"]:
        if len({str(m[key]) for m in metas}) != 1:
            raise ValueError(f"Incompatible dataset field: {key}")
    if len({m["seed"] for m in metas}) != 3:
        raise ValueError("Train, validation and test must use distinct seeds")
    for d in datasets[1:]:
        if not np.array_equal(d[3], datasets[0][3]) or not np.array_equal(d[4], datasets[0][4]):
            raise ValueError("Affinity lookup tables differ between splits")
    topology_hashes = {}
    for name in ["circuit.stim", "raw.dem", "decomposed.dem", "graph_edges.csv", "coords.csv"]:
        hashes = [hashlib.sha256((p / name).read_bytes()).hexdigest()
                  for p in [args.train, args.validation, args.test]]
        if len(set(hashes)) != 1:
            raise ValueError(f"Fixed circuit/graph differs between splits: {name}")
        topology_hashes[name] = hashes[0]
    n = metas[0]["num_detectors"]
    boundary, relevant = graph_geometry(args.train, n)
    scale = metas[0]["distance"] * np.log(1 / metas[0]["p"])
    duplicate_sets = []
    for split, path, dataset in zip(["train", "validation", "test"],
                                     [args.train, args.validation, args.test], datasets):
        meta, syndromes, labels, affinity, hops = dataset
        np.save(args.output / f"{split}_labels.npy", labels)
        np.save(args.output / f"{split}_syndrome.npy", syndromes.astype(np.float32))
        images = {kind: np.lib.format.open_memmap(args.output / f"{split}_{kind}.npy", mode="w+",
                  dtype="float32", shape=(len(syndromes), 7 if kind == "context" else 3, args.size, args.size))
                  for kind in ["hop", "weighted", "boundary", "context"]}
        summaries = np.empty((len(syndromes), len(SUMMARY_NAMES)), dtype=np.float32)
        t0 = time.perf_counter()
        for i, s in enumerate(syndromes):
            ids = np.flatnonzero(s)
            summaries[i] = summary_features(ids, affinity, boundary, relevant, scale)
            images["hop"][i] = image_features(ids, hops, args.size, scale=metas[0]["distance"], hop=True)
            images["weighted"][i] = image_features(ids, affinity, args.size, scale=scale)
            images["boundary"][i] = image_features(ids, affinity, args.size, boundary, scale=scale)
            images["context"][i] = context_features(ids, affinity, args.size, boundary, scale)
        for array in images.values():
            array.flush()
        np.save(args.output / f"{split}_summary.npy", summaries)
        hashes = {hashlib.sha256(np.packbits(s).tobytes()).digest() for s in syndromes}
        duplicate_sets.append(hashes)
        print(f"{split}: {len(syndromes)} shots, max weight={syndromes.sum(1).max()}, "
              f"{len(hashes)} unique syndromes; features {time.perf_counter()-t0:.1f}s", flush=True)
    # Repeated syndromes are physically legitimate, especially empty syndromes;
    # report overlap and additionally score novel test syndromes.
    train_hashes = duplicate_sets[0]
    novel = np.asarray([hashlib.sha256(np.packbits(s).tobytes()).digest() not in train_hashes
                        for s in datasets[2][1]])
    np.save(args.output / "test_novel.npy", novel)
    dump(args.output / "protocol.json", {"datasets": [str(p.resolve()) for p in [args.train, args.validation, args.test]],
        "metadata": metas, "topology_sha256": topology_hashes,
        "matrix_size": args.size, "log_distance_scale": scale,
        "summary_names": SUMMARY_NAMES, "test_novel_shots": int(novel.sum()),
        "boundary_relevant_detectors": int(relevant.sum()),
        "feature_inputs": "syndrome and fixed DEM only; labels, truth, predictions, decoder timing excluded",
        "shared_boundary": "always present; symmetrized from detector-to-boundary affinity",
        "normalization": "fixed d*ln(1/p) for log-affinity/log-odds, d for hops; no per-shot normalization"})


class GapCNN(nn.Module):
    def __init__(self, channels=3):
        super().__init__()
        self.features = nn.Sequential(nn.Conv2d(channels, 8, 3, stride=2, padding=1), nn.ReLU(),
            nn.AvgPool2d(2), nn.Conv2d(8, 16, 3, padding=1), nn.ReLU(),
            nn.AvgPool2d(2), nn.AdaptiveAvgPool2d((6, 6)), nn.Flatten())
        self.head = nn.Sequential(nn.Linear(16 * 6 * 6, 64), nn.ReLU(), nn.Linear(64, 2))

    def forward(self, x):
        return self.head(self.features(x))


class GapMLP(nn.Module):
    """Small neural control for whether summary extraction is already sufficient."""
    def __init__(self, inputs):
        super().__init__()
        self.layers = nn.Sequential(nn.Linear(inputs, 64), nn.ReLU(),
                                    nn.Linear(64, 64), nn.ReLU(), nn.Linear(64, 2))

    def forward(self, x):
        return self.layers(x)


def regression_metrics(y, pred):
    if not len(y):
        return {key: None for key in ["mae_db", "rmse_db", "r2", "spearman",
                "abs_error_p99_db", "max_abs_error_db", "max_overestimate_db"]}
    err = pred - y
    denom = np.sum((y - y.mean()) ** 2)
    rho = float(spearmanr(y, pred).statistic) if np.std(pred) > 1e-9 and denom > 0 else float("nan")
    return {"mae_db": float(np.abs(err).mean()), "rmse_db": float(np.sqrt((err ** 2).mean())),
            "r2": float(1 - np.sum(err ** 2) / denom) if denom > 0 else None,
            "spearman": rho if np.isfinite(rho) else None,
            "abs_error_p99_db": float(np.quantile(np.abs(err), .99)),
            "max_abs_error_db": float(np.abs(err).max()),
            "max_overestimate_db": float(np.maximum(err, 0).max())}


def score_predictions(y, pred, validation_y, validation_pred, novel, seed):
    result = {}
    for j, name in enumerate(["complementary", "cluster"]):
        yy, pp = y[:, j], pred[:, j]
        value = regression_metrics(yy, pp)
        value["novel_syndromes"] = regression_metrics(yy[novel], pp[novel]) if novel.any() else None
        # Threshold is defined by validation labels only, then frozen. The
        # prediction threshold targets >=95% validation recall of low-gap shots.
        threshold = float(np.quantile(validation_y[:, j], .10))
        val_hard = validation_y[:, j] <= threshold
        pred_threshold = float(np.quantile(validation_pred[val_hard, j], .95))
        hard, route = yy <= threshold, pp <= pred_threshold
        tp = np.count_nonzero(hard & route)
        value["routing"] = {"gap_threshold_db": threshold, "prediction_threshold_db": pred_threshold,
            "hard_count": int(hard.sum()), "hard_prevalence": float(hard.mean()),
            "recall": float(tp / hard.sum()) if hard.any() else None,
            "precision": float(tp / route.sum()) if route.any() else None,
            "routed_fraction": float(route.mean()),
            "roc_auc": float(roc_auc_score(hard, -pp)) if 0 < hard.sum() < len(hard) else None,
            "average_precision": float(average_precision_score(hard, -pp)) if hard.any() else None}
        # Paired pointwise errors are saved for later model comparisons; these
        # intervals describe a fixed trained model, not training-run variation.
        rng = np.random.default_rng(seed)
        boot = [np.abs(pp[idx] - yy[idx]).mean()
                for idx in rng.integers(0, len(yy), size=(300, len(yy)))]
        value["mae_95pct_bootstrap_db"] = np.quantile(boot, [.025, .975]).tolist()
        result[name] = value
    return result


def train(args):
    torch.set_num_threads(args.threads)
    torch.set_num_interop_threads(1)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    args.output.mkdir(parents=True, exist_ok=False)
    protocol = json.loads((args.prepared / "protocol.json").read_text())
    ys = [np.load(args.prepared / f"{s}_labels.npy") for s in ["train", "validation", "test"]]
    ytrain, yval, ytest = [y[:, :2] for y in ys]
    novel = np.load(args.prepared / "test_novel.npy")
    # Validation controls stopping; test labels are accessed for scoring only.
    t0 = time.perf_counter()
    history = []
    if args.model.startswith("cnn_") or args.model.startswith("mlp_"):
        feature = ("summary" if args.model == "mlp_boundary" else
                   args.model.removeprefix("mlp_") if args.model.startswith("mlp_") else
                   args.model.removeprefix("cnn_"))
        xs = [np.load(args.prepared / f"{s}_{feature}.npy", mmap_mode="r")
              for s in ["train", "validation", "test"]]
        input_mean = input_std = None
        if args.model.startswith("mlp_"):
            # Triangle entries have a known [0,1] scale. Avoid amplifying rare
            # slots (including hash stash slots) using tiny sample deviations.
            floor = .05 if feature.endswith("triangle") else 1e-6
            input_mean, input_std = xs[0].mean(0), np.maximum(xs[0].std(0), floor)
            xs = [((x - input_mean) / input_std).astype(np.float32) for x in xs]
            model = GapMLP(xs[0].shape[1])
        else:
            model = GapCNN(channels=xs[0].shape[1])
        mean, std = ytrain.mean(axis=0), ytrain.std(axis=0)
        std = np.maximum(std, 1.)
        target = torch.from_numpy(((ytrain - mean) / std).astype(np.float32))
        val_target = torch.from_numpy(((yval - mean) / std).astype(np.float32))
        optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
        best, best_epoch, best_state = float("inf"), 0, None
        rng = np.random.default_rng(args.seed)

        def objective(prediction, truth):
            return nn.functional.mse_loss(prediction, truth)

        def predict(x):
            model.eval()
            with torch.no_grad():
                return np.concatenate([model(torch.from_numpy(np.array(x[i:i+args.batch_size]))).numpy()
                    for i in range(0, len(x), args.batch_size)]) * std + mean

        for epoch in range(args.epochs):
            model.train()
            order = rng.permutation(len(xs[0]))
            train_loss = 0.
            for start in range(0, len(order), args.batch_size):
                ix = order[start:start+args.batch_size]
                x = torch.from_numpy(np.array(xs[0][ix]))
                optimizer.zero_grad()
                loss = objective(model(x), target[ix])
                loss.backward()
                optimizer.step()
                train_loss += float(loss.detach()) * len(ix)
            vp = predict(xs[1])
            val_loss = float(objective(torch.from_numpy(((vp - mean) / std).astype(np.float32)), val_target))
            history.append({"epoch": epoch + 1, "train_loss": train_loss / len(order), "validation_loss": val_loss})
            print(json.dumps(history[-1]), flush=True)
            if val_loss < best - 1e-5:
                best, best_epoch, best_state = val_loss, epoch + 1, copy.deepcopy(model.state_dict())
            elif epoch + 1 - best_epoch >= args.patience:
                break
        model.load_state_dict(best_state)
        predval, predtest = predict(xs[1]), predict(xs[2])
        torch.save({"state_dict": model.state_dict(), "mean": mean.tolist(), "std": std.tolist(),
                    "feature": feature, "model": args.model,
                    "input_mean": None if input_mean is None else input_mean.tolist(),
                    "input_std": None if input_std is None else input_std.tolist(),
                    "protocol": protocol}, args.output / "model.pt")
        parameters = sum(p.numel() for p in model.parameters())
        dump(args.output / "history.json", history)
    else:
        xs = [np.load(args.prepared / f"{s}_summary.npy") for s in ["train", "validation", "test"]]
        if args.model == "constant":
            predval = np.tile(ytrain.mean(0), (len(yval), 1))
            predtest = np.tile(ytrain.mean(0), (len(ytest), 1))
            model = None
        else:
            if args.model == "weight":
                xs = [x[:, :1] for x in xs]
            elif args.model.endswith("affinity"):
                xs = [x[:, :15] for x in xs]
            elif args.model == "trees_basis":
                xs = [x[:, :17] for x in xs]
            elif args.model == "trees_syndrome":
                xs = [np.load(args.prepared / f"{s}_syndrome.npy") for s in ["train", "validation", "test"]]
            if args.model == "ridge_affinity":
                model = make_pipeline(StandardScaler(), Ridge(alpha=10.))
            else:
                model = ExtraTreesRegressor(n_estimators=128, min_samples_leaf=8,
                                             max_features="sqrt" if args.model == "trees_syndrome" else 1.,
                                             n_jobs=args.threads, random_state=args.seed)
            with threadpool_limits(limits=args.threads):
                model.fit(xs[0], ytrain)
                predval, predtest = model.predict(xs[1]), model.predict(xs[2])
            import joblib
            joblib.dump(model, args.output / "model.joblib")
            if hasattr(model, "feature_importances_") and xs[0].shape[1] == len(SUMMARY_NAMES):
                dump(args.output / "feature_importance.json", dict(zip(SUMMARY_NAMES, model.feature_importances_.tolist())))
        parameters, best_epoch = None, None
    elapsed = time.perf_counter() - t0
    # Physical gaps are nonnegative. Apply the same projection to every model.
    predval, predtest = np.maximum(predval, 0.), np.maximum(predtest, 0.)
    np.save(args.output / "validation_predictions.npy", predval)
    np.save(args.output / "test_predictions.npy", predtest)
    scores = score_predictions(ytest, predtest, yval, predval, novel, args.seed)
    weights = np.load(args.prepared / "test_summary.npy", mmap_mode="r")[:, 0].astype(int)
    for j, name in enumerate(["complementary", "cluster"]):
        for parity, label in [(0, "even_hamming_weight"), (1, "odd_hamming_weight")]:
            mask = weights % 2 == parity
            scores[name][label] = {"shots": int(mask.sum()),
                **regression_metrics(ytest[mask, j], predtest[mask, j])}
    scores.update(model=args.model, seed=args.seed, training_seconds=elapsed,
                  parameters=parameters, best_epoch=best_epoch, torch_version=torch.__version__,
                  configuration={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                  threads=args.threads, protocol=str((args.prepared / "protocol.json").resolve()))
    dump(args.output / "metrics.json", scores)
    print(json.dumps(scores, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    for name in ["train", "validation", "test", "output"]:
        p.add_argument("--" + name, required=True, type=Path)
    p.add_argument("--size", type=int, default=48)
    p = sub.add_parser("train")
    p.add_argument("--prepared", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--model", required=True, choices=["constant", "weight", "ridge_affinity", "trees_affinity",
        "trees_basis", "trees_boundary", "trees_syndrome", "cnn_hop", "cnn_weighted", "cnn_boundary", "cnn_context", "mlp_boundary", "mlp_syndrome", "mlp_packed_triangle", "mlp_hashed_triangle"])
    p.add_argument("--seed", type=int, default=101)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--epochs", type=int, default=40)
    p.add_argument("--patience", type=int, default=8)
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--lr", type=float, default=.001)
    args = parser.parse_args()
    if args.command == "prepare":
        if args.size < 8:
            parser.error("--size must be at least 8 for the CNN pooling layers")
        prepare(args)
    else:
        if min(args.threads, args.epochs, args.patience, args.batch_size) < 1:
            parser.error("threads, epochs, patience, batch size must be positive")
        train(args)


if __name__ == "__main__":
    main()
