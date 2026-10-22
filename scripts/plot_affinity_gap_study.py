#!/usr/bin/env python3
"""Summarize and plot prepared affinity-gap prediction runs."""

import argparse
import csv
import json
import os
from pathlib import Path

_ARTIFACT_ROOT = Path(__file__).resolve().parents[1] / ".agentic_artifacts" / "affinity_gap"
(_ARTIFACT_ROOT / "mpl").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_ARTIFACT_ROOT / "mpl"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr


TARGETS = ("complementary", "cluster")
BASELINES = ("constant", "weight", "trees_affinity", "trees_boundary")
CSV_METRICS = ("mae_db", "rmse_db", "r2", "spearman", "lowgap_recall",
               "lowgap_precision", "routed_fraction", "roc_auc", "average_precision")


def load_run(path, prepared):
    path = path.expanduser().resolve()
    metrics_path = path / "metrics.json"
    predictions_path = path / "test_predictions.npy"
    validation_path = path / "validation_predictions.npy"
    for required in (metrics_path, predictions_path, validation_path):
        if not required.is_file():
            raise ValueError(f"Missing run artifact: {required}")
    metrics = json.loads(metrics_path.read_text())
    if metrics.get("protocol") and Path(metrics["protocol"]).resolve() != (prepared / "protocol.json").resolve():
        raise ValueError(f"Run uses a different prepared dataset: {path}")
    protocol = json.loads((prepared / "protocol.json").read_text())
    if metrics.get("holdout_dataset") and Path(metrics["holdout_dataset"]).resolve() != Path(protocol["datasets"][-1]).resolve():
        raise ValueError(f"Run uses a different holdout dataset: {path}")
    test_pred = np.load(predictions_path)
    val_pred = np.load(validation_path)
    if test_pred.ndim != 2 or test_pred.shape[1] != 2 or val_pred.ndim != 2 or val_pred.shape[1] != 2:
        raise ValueError(f"Predictions must be N x 2 arrays: {path}")
    return {"path": path, "metrics": metrics, "test_pred": test_pred, "val_pred": val_pred,
            "model": str(metrics.get("model", path.name)), "seed": metrics.get("seed", "")}


def write_metrics_csv(runs, path):
    columns = ["run", "model", "training_seed", "target", *CSV_METRICS,
               "novel_mae_db", "novel_rmse_db", "novel_r2", "novel_spearman"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for run in runs:
            for target in TARGETS:
                scores = run["metrics"].get(target, {})
                routing = scores.get("routing", {})
                novel = scores.get("novel_syndromes") or {}
                row = {"run": str(run["path"]), "model": run["model"],
                       "training_seed": run["seed"], "target": target,
                       "mae_db": scores.get("mae_db"), "rmse_db": scores.get("rmse_db"),
                       "r2": scores.get("r2"), "spearman": scores.get("spearman"),
                       "lowgap_recall": routing.get("recall"),
                       "lowgap_precision": routing.get("precision"),
                       "routed_fraction": routing.get("routed_fraction"),
                       "roc_auc": routing.get("roc_auc"),
                       "average_precision": routing.get("average_precision"),
                       "novel_mae_db": novel.get("mae_db"),
                       "novel_rmse_db": novel.get("rmse_db"),
                       "novel_r2": novel.get("r2"),
                       "novel_spearman": novel.get("spearman")}
                writer.writerow(row)


def select_best_cnn(runs, validation_y):
    cnn = [run for run in runs if run["model"].startswith("cnn_")]
    if not cnn:
        return None
    for run in cnn:
        if run["val_pred"].shape != validation_y.shape:
            raise ValueError(f"Validation prediction shape mismatch for {run['path']}")
        run["validation_mae"] = float(np.mean(np.abs(run["val_pred"] - validation_y)))
    return min(cnn, key=lambda run: (run["validation_mae"], run["model"], str(run["path"])))


def plot_metrics(runs, output):
    models = sorted({run["model"] for run in runs})
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(1, len(models))))
    fig, axes = plt.subplots(1, 2, figsize=(max(11, len(models) * .95), 5), constrained_layout=True)
    for j, target in enumerate(TARGETS):
        ax = axes[j]
        x = np.arange(len(models))
        width = .36
        for offset, metric, label, color in [(-width / 2, "mae_db", "MAE", "#2878b5"),
                                             (width / 2, "rmse_db", "RMSE", "#e07a27")]:
            means, stds = [], []
            for model in models:
                values = [run["metrics"][target].get(metric) for run in runs
                          if run["model"] == model and run["metrics"].get(target, {}).get(metric) is not None]
                means.append(float(np.mean(values)) if values else np.nan)
                stds.append(float(np.std(values, ddof=1)) if len(values) > 1 else 0.)
            ax.bar(x + offset, means, width, yerr=stds, capsize=2, label=label, color=color,
                   alpha=.88, error_kw={"linewidth": .8})
        ax.set_title(f"{target.capitalize()} gap")
        ax.set_ylabel("Error (dB)")
        ax.set_xticks(x, models, rotation=45, ha="right")
        ax.grid(axis="y", alpha=.25)
        if j == 0:
            ax.legend(frameon=False)
    fig.suptitle("Test regression metrics; bars show run mean ± seed standard deviation")
    fig.savefig(output, dpi=170)
    plt.close(fig)


def equal_population_calibration(y, pred, bins=12):
    order = np.argsort(pred, kind="stable")
    groups = np.array_split(order, min(bins, len(order)))
    groups = [g for g in groups if len(g)]
    return (np.asarray([pred[g].mean() for g in groups]),
            np.asarray([y[g].mean() for g in groups]),
            np.asarray([y[g].std() / np.sqrt(len(g)) for g in groups]))


def plot_calibration(best, test_y, output, seed=42, max_points=5000):
    rng = np.random.default_rng(seed)
    n = len(test_y)
    selected = np.arange(n) if n <= max_points else np.sort(rng.choice(n, max_points, replace=False))
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    for j, target in enumerate(TARGETS):
        y, pred = test_y[:, j], best["test_pred"][:, j]
        ax = axes[j]
        ax.scatter(y[selected], pred[selected], s=8, alpha=.18, rasterized=True)
        px, ty, uncertainty = equal_population_calibration(y, pred)
        ax.errorbar(ty, px, xerr=uncertainty, fmt="o-", color="#c23b22", lw=1.5,
                    ms=4, label="Equal-population bins (x error: SE)")
        lo = min(float(np.min(y)), float(np.min(pred)))
        hi = max(float(np.max(y)), float(np.max(pred)))
        ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="Ideal")
        ax.set_xlabel("True gap (dB)")
        ax.set_ylabel("Predicted gap (dB)")
        ax.set_title(target.capitalize())
        ax.grid(alpha=.2)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"Validation-selected CNN test predictions: {best['model']} seed={best['seed']}")
    fig.savefig(output, dpi=170)
    plt.close(fig)


def routing_curve(y, pred, validation_y, target_index, points=101):
    lowgap_cutoff = float(np.quantile(validation_y[:, target_index], .10))
    hard = y[:, target_index] <= lowgap_cutoff
    if not hard.any():
        return np.linspace(0, 1, points), np.full(points, np.nan)
    thresholds = np.quantile(pred[:, target_index], np.linspace(0, 1, points))
    routed_fraction, recall = [], []
    for threshold in thresholds:
        route = pred[:, target_index] <= threshold
        routed_fraction.append(float(route.mean()))
        recall.append(float(np.count_nonzero(route & hard) / hard.sum()))
    x = np.asarray(routed_fraction)
    r = np.asarray(recall)
    unique_x, inverse = np.unique(x, return_inverse=True)
    best_recall = np.zeros(len(unique_x))
    np.maximum.at(best_recall, inverse, r)
    if unique_x[0] > 0:
        unique_x = np.insert(unique_x, 0, 0.)
        best_recall = np.insert(best_recall, 0, 0.)
    return unique_x, best_recall


def plot_routing(runs, test_y, validation_y, output):
    families = sorted({run["model"] for run in runs})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for j, target in enumerate(TARGETS):
        ax = axes[j]
        for model in families:
            curves = [routing_curve(test_y, run["test_pred"], validation_y, j)
                      for run in runs if run["model"] == model]
            if not curves:
                continue
            grid = np.linspace(0, 1, 101)
            interpolated = [np.interp(grid, x, y, left=y[0], right=y[-1]) for x, y in curves]
            ax.plot(grid, np.mean(interpolated, axis=0), label=model, lw=1.5)
        ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Fraction routed to second decoder",
               ylabel="Recall of low-gap shots", title=target.capitalize())
        ax.grid(alpha=.25)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Descriptive test routing curves; low-gap cutoff fixed from validation labels")
    fig.savefig(output, dpi=170)
    plt.close(fig)


def plot_weight_and_histograms(test_y, syndrome, output_hex, output_hist):
    weights = syndrome.sum(axis=1)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for j, target in enumerate(TARGETS):
        hb = axes[j].hexbin(weights, test_y[:, j], gridsize=40, mincnt=1, bins="log", cmap="viridis")
        axes[j].set(title=target.capitalize(), xlabel="Syndrome weight", ylabel="True gap (dB)")
        fig.colorbar(hb, ax=axes[j], label="Shots (log color scale)")
    fig.suptitle("Syndrome weight versus decoder gap")
    fig.savefig(output_hex, dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for j, target in enumerate(TARGETS):
        axes[j].hist(test_y[:, j], bins=50, color=("#2878b5" if j == 0 else "#e07a27"), alpha=.85)
        axes[j].set(title=target.capitalize(), xlabel="Gap (dB)", ylabel="Shots")
        axes[j].grid(axis="y", alpha=.2)
    fig.suptitle("Test-set gap distributions")
    fig.savefig(output_hist, dpi=170)
    plt.close(fig)


def plot_affinity_intuition(test_y, summary, summary_names, output):
    required = ("weight", "mean_effective_partners", "mean_nearest_second_margin")
    missing = [name for name in required if name not in summary_names]
    if missing:
        raise ValueError("protocol.summary_names missing required fields: " + ", ".join(missing))
    if summary.ndim != 2 or summary.shape[1] != len(summary_names) or len(summary) != len(test_y):
        raise ValueError("test_summary shape does not match protocol names/test labels")
    columns = {name: summary[:, summary_names.index(name)] for name in required}
    nonzero = columns["weight"] > 0
    if not nonzero.any():
        raise ValueError("test set has no nonzero syndrome weight for fixed-weight plot")
    values, counts = np.unique(columns["weight"][nonzero], return_counts=True)
    mode_weight = values[np.argmax(counts)]
    fixed = columns["weight"] == mode_weight
    fixed_n = int(fixed.sum())

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    features = ("mean_effective_partners", "mean_nearest_second_margin")
    for col, feature in enumerate(features):
        x = columns[feature]
        all_hex = axes[0, col].hexbin(x, test_y[:, 0], gridsize=45, mincnt=1,
                                      bins="log", cmap="viridis")
        axes[0, col].set_title(f"All test shots: {feature}")
        axes[0, col].set_ylabel("Complementary gap (dB)")
        axes[0, col].set_xlabel(feature.replace("_", " "))
        fig.colorbar(all_hex, ax=axes[0, col], label="Shots (log color scale)")

        fixed_hex = axes[1, col].hexbin(x[fixed], test_y[fixed, 0], gridsize=35, mincnt=1,
                                        bins="log", cmap="viridis")
        rho = float(spearmanr(x[fixed], test_y[fixed, 0]).statistic) if fixed_n > 1 else float("nan")
        rho_text = f"{rho:.3f}" if np.isfinite(rho) else "n/a"
        axes[1, col].set_title(f"Syndrome weight {mode_weight:g} (n={fixed_n}); Spearman ρ={rho_text}")
        axes[1, col].set_ylabel("Complementary gap (dB)")
        axes[1, col].set_xlabel(feature.replace("_", " "))
        fig.colorbar(fixed_hex, ax=axes[1, col], label="Shots (log color scale)")
    fig.suptitle("Affinity summaries versus complementary gap")
    fig.savefig(output, dpi=170)
    plt.close(fig)


def metric_table(runs):
    lines = ["| Model | Runs | Complementary MAE / RMSE (dB) | Cluster MAE / RMSE (dB) |",
             "|---|---:|---:|---:|"]
    for model in sorted({run["model"] for run in runs}):
        selected = [run for run in runs if run["model"] == model]
        values = []
        for target in TARGETS:
            mae = [run["metrics"][target]["mae_db"] for run in selected
                   if run["metrics"].get(target, {}).get("mae_db") is not None]
            rmse = [run["metrics"][target]["rmse_db"] for run in selected
                    if run["metrics"].get(target, {}).get("rmse_db") is not None]
            if mae:
                values.append(f"{np.mean(mae):.3g} / {np.mean(rmse):.3g}")
            else:
                values.append("n/a")
        lines.append(f"| {model} | {len(selected)} | {values[0]} | {values[1]} |")
    return "\n".join(lines)


def write_report(output, prepared, protocol, runs, best, test_novel):
    data_meta = protocol.get("metadata", [{}])[-1]
    dataset_paths = protocol.get("datasets", [])
    rows = ["# Affinity gap prediction study", "",
            f"Prepared data: `{prepared.resolve()}`",
            f"Protocol: `{(prepared / 'protocol.json').resolve()}`",
            f"Distance/rounds/noise: `{data_meta.get('distance', 'unknown')}` / `"
            f"{data_meta.get('rounds', 'unknown')}` / `p={data_meta.get('p', 'unknown')}`.",
            f"Memory basis: `{data_meta.get('memory_basis', 'unknown')}`; opposite-basis detectors: "
            f"`{data_meta.get('include_opposite_basis_detectors', 'unknown')}`.",
            "Labels are unsigned complementary MWPM gap and cluster-gap dB from uncorrelated PyMatching. "
            "Cluster gap retains the two physical boundary classes; shared boundary features remain in the fixed DEM input.",
            "This report evaluates one prepared code distance and makes no cross-distance generalization claim.", "",
            "## Input datasets", ""]
    rows.extend(f"- `{p}`" for p in dataset_paths)
    rows.extend(["", "## Test regression summary", "", metric_table(runs), "",
                 "Requested baselines present: " + ", ".join(
                     f"`{name}` ({'yes' if name in {run['model'] for run in runs} else 'no'})"
                     for name in BASELINES) + ".",
                 f"Test shots: {len(test_novel)}; novel syndromes relative to training: "
                 f"{int(test_novel.sum())} ({100 * test_novel.mean():.2f}%).",
                 f"Run directories: {len(runs)}."])
    if best:
        rows.extend(["", "## Calibration model selection", "",
                     f"Selected `{best['model']}` seed `{best['seed']}` by lowest mean validation MAE across both gaps "
                     f"({best['validation_mae']:.4g} dB). Test labels were not used to select this CNN."])
    rows.extend(["", "## Figures", "",
                 "- `metrics_comparison.png`: test MAE/RMSE by model, with run-to-run standard deviation.",
                 "- `best_cnn_calibration.png`: true versus predicted test gaps and equal-population calibration bins.",
                 "- `routing_curves.png`: descriptive recall-versus-routed-fraction curves. The low-gap label cutoff is "
                 "fixed from validation labels; the curves use test predictions and labels for retrospective description "
                 "only and do not select an operating threshold.",
                 "- `syndrome_weight_gap.png` and `gap_histograms.png`: test label diagnostics.",
                 "- `affinity_intuition.png`: effective-partner count and nearest-second margin versus complementary gap, "
                 "before and after restricting to the test set's most common nonzero syndrome weight.", "",
                 "## Run artifacts", ""])
    rows.extend(f"- `{run['model']}` seed `{run['seed']}`: `{run['path']}`" for run in runs)
    (output / "report.md").write_text("\n".join(rows) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--runs", required=True, nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    prepared = args.prepared.expanduser().resolve()
    output = args.output.expanduser().resolve()
    protocol_path = prepared / "protocol.json"
    required = [protocol_path, prepared / "test_labels.npy", prepared / "validation_labels.npy",
                prepared / "test_novel.npy", prepared / "test_syndrome.npy", prepared / "test_summary.npy"]
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        parser.error("missing prepared artifacts: " + ", ".join(missing))
    if output.exists():
        parser.error(f"output path already exists: {output}")
    protocol = json.loads(protocol_path.read_text())
    test_labels = np.load(prepared / "test_labels.npy")
    validation_labels = np.load(prepared / "validation_labels.npy")
    test_novel = np.load(prepared / "test_novel.npy").astype(bool)
    syndrome = np.load(prepared / "test_syndrome.npy", mmap_mode="r")
    test_summary = np.load(prepared / "test_summary.npy", mmap_mode="r")
    if test_labels.ndim != 2 or test_labels.shape[1] < 2 or validation_labels.ndim != 2 or validation_labels.shape[1] < 2:
        parser.error("prepared labels must have at least two columns")
    test_y, validation_y = test_labels[:, :2], validation_labels[:, :2]
    if len(test_y) != len(test_novel) or len(test_y) != len(syndrome):
        parser.error("test labels, novelty mask, and syndromes have inconsistent row counts")
    runs = [load_run(p, prepared) for p in args.runs]
    if not runs:
        parser.error("provide at least one run directory")
    for run in runs:
        if run["test_pred"].shape != test_y.shape:
            parser.error(f"test prediction shape mismatch: {run['path']}")
    output.mkdir(parents=True, exist_ok=False)
    write_metrics_csv(runs, output / "metrics.csv")
    plot_metrics(runs, output / "metrics_comparison.png")
    best = select_best_cnn(runs, validation_y)
    if best:
        plot_calibration(best, test_y, output / "best_cnn_calibration.png")
    plot_routing(runs, test_y, validation_y, output / "routing_curves.png")
    plot_weight_and_histograms(test_y, syndrome, output / "syndrome_weight_gap.png",
                               output / "gap_histograms.png")
    summary_names = protocol.get("summary_names", [])
    plot_affinity_intuition(test_y, test_summary, summary_names, output / "affinity_intuition.png")
    write_report(output, prepared, protocol, runs, best, test_novel)
    print(f"Wrote affinity-gap study report and figures to {output}")


if __name__ == "__main__":
    main()
