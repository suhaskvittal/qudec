#!/usr/bin/env python3
"""Contributor: OpenAI Codex. Analyze per-rank PDZ + Libra gzip records."""

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import struct

import numpy as np


HEADER = struct.Struct("<8sIIQQII")
FEATURES = np.dtype([
    ("input_hw", "<u4"), ("post_hw", "<u4"), ("gap_db", "<f8"),
    ("ensemble", "u1"), ("changed", "u1"), ("initial", "u1"),
    ("final", "u1"), ("truth", "u1"), ("bypassed", "u1"),
])


def load_records(directory):
    """Validate gzip integrity, packed syndromes, flags, and per-rank headers."""
    files = sorted(directory.glob("rank_*.bin.gz"))
    if not files:
        raise ValueError(f"No rank files in {directory}")
    parts, ranks, geometry = [], {}, None
    popcount = np.array([i.bit_count() for i in range(256)], dtype=np.uint8)
    for path in files:
        with gzip.open(path, "rb") as source:
            header = source.read(HEADER.size)
            if len(header) != HEADER.size:
                raise ValueError(f"Truncated header: {path}")
            magic, version, width, detectors, observables, rank, header_size = HEADER.unpack(header)
            nbytes = (detectors + 7) // 8
            if (magic != b"PDZLIBR1" or version not in (1, 2) or header_size != 40
                    or observables != 1 or width != 24 + 2 * nbytes):
                raise ValueError(f"Unsupported format: {path}")
            if rank in ranks or (geometry is not None and geometry != (detectors, observables)):
                raise ValueError(f"Duplicate rank or inconsistent detector geometry: {path}")
            geometry = detectors, observables
            payload = source.read()  # Reading to EOF also checks the gzip CRC.
        if len(payload) % width:
            raise ValueError(f"Truncated record: {path}")
        dtype = np.dtype({
            "names": FEATURES.names + ("input_bits", "post_bits"),
            "formats": [FEATURES.fields[n][0] for n in FEATURES.names]
                       + [("u1", (nbytes,)), ("u1", (nbytes,))],
            "offsets": [0, 4, 8, 16, 17, 18, 19, 20, 21, 24, 24 + nbytes],
            "itemsize": width,
        })
        records = np.frombuffer(payload, dtype=dtype)
        for hw, bits in (("input_hw", "input_bits"), ("post_hw", "post_bits")):
            if np.any(popcount[records[bits]].sum(axis=1) != records[hw]):
                raise ValueError(f"Recorded Hamming weight does not match packed bits: {path}")
            if detectors % 8 and np.any(records[bits][:, -1] >> (detectors % 8)):
                raise ValueError(f"Nonzero syndrome padding: {path}")
        if np.any(records["post_bits"] & ~records["input_bits"]):
            raise ValueError(f"PDZ introduced a detector event: {path}")
        for name in ("ensemble", "changed", "initial", "final", "truth", "bypassed"):
            if np.any(records[name] > 1):
                raise ValueError(f"Invalid boolean field {name}: {path}")
        expected = (records["ensemble"] != 0) & (records["initial"] != records["final"])
        if np.any(expected != records["changed"].astype(bool)):
            raise ValueError(f"Inconsistent frame-change flag: {path}")
        bypassed = records["bypassed"] != 0
        if version == 1 and np.any(bypassed):
            raise ValueError(f"Unexpected bypass flag in version 1 file: {path}")
        if np.any(bypassed & ((records["ensemble"] != 0) | (records["initial"] != records["final"]))):
            raise ValueError(f"Bypassed shot did not return original PyMatching: {path}")
        measured_gap = records["gap_db"][~bypassed]
        if (np.any(~np.isfinite(measured_gap)) or np.any(measured_gap < 0)
                or np.any(~np.isnan(records["gap_db"][bypassed]))):
            raise ValueError(f"Invalid complementary gap: {path}")
        part = np.empty(len(records), dtype=FEATURES)
        for name in FEATURES.names:
            part[name] = records[name]
        parts.append(part)
        ranks[rank] = len(records)
    return np.concatenate(parts), ranks, geometry


def group_table(data, keys, names):
    key_array = np.column_stack(keys).astype(np.int64)
    values, inverse = np.unique(key_array, axis=0, return_inverse=True)
    total = np.bincount(inverse)
    ensemble = np.bincount(inverse, weights=data["ensemble"]).astype(np.int64)
    changed = np.bincount(inverse, weights=data["changed"]).astype(np.int64)
    corrected = (data["changed"] != 0) & (data["initial"] != data["truth"])
    harmed = (data["changed"] != 0) & (data["final"] != data["truth"])
    helped = np.bincount(inverse, weights=corrected).astype(np.int64)
    hurt = np.bincount(inverse, weights=harmed).astype(np.int64)
    bypassed = np.bincount(inverse, weights=data["bypassed"]).astype(np.int64)
    nonempty_residual = (data["post_hw"] > 0) & (data["bypassed"] == 0)
    gap_count = np.bincount(inverse, weights=nonempty_residual)
    gap_sum = np.bincount(inverse, weights=np.where(nonempty_residual, data["gap_db"], 0))
    rows = []
    for i, key in enumerate(values):
        row = dict(zip(names, map(int, key)))
        row.update(shots=int(total[i]), ensemble=int(ensemble[i]), changes=int(changed[i]),
                   ensemble_probability=float(ensemble[i] / total[i]),
                   change_probability_given_ensemble=float(changed[i] / ensemble[i]) if ensemble[i] else None,
                   helped=int(helped[i]), harmed=int(hurt[i]), bypassed=int(bypassed[i]),
                   mean_gap_db_nonempty_residual=float(gap_sum[i] / gap_count[i]) if gap_count[i] else None)
        rows.append(row)
    return rows


def write_csv(path, rows):
    with path.open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def correlation(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if not len(x):
        return None
    x, y = x - x.mean(), y - y.mean()
    denominator = math.sqrt(float(x @ x) * float(y @ y))
    return float(x @ y / denominator) if denominator else None


def controlled_logistic(data):
    """Fit ensemble probability to initial HW and removal fraction, on grouped counts.

    No regularization or threshold selection. Newton steps use backtracking;
    report nonconvergence rather than interpreting a divergent fit.
    """
    removed = data["input_hw"].astype(np.int64) - data["post_hw"]
    rows = group_table(data, [data["input_hw"], removed], ["hw", "removed"])
    hw = np.array([r["hw"] for r in rows], dtype=float)
    removed = np.array([r["removed"] for r in rows], dtype=float)
    n = np.array([r["shots"] for r in rows], dtype=float)
    y = np.array([r["ensemble"] for r in rows], dtype=float)
    if y.sum() == 0 or y.sum() == n.sum():
        return {"converged": False, "reason": "No variation in ensemble use"}
    # Slopes correspond to ten input events and ten percentage points removed.
    x = np.column_stack([np.ones(len(hw)), hw / 10, (removed / hw) * 10])
    beta = np.array([math.log(y.sum() / (n.sum() - y.sum())), 0., 0.])
    converged = False
    for _ in range(100):
        z = x @ beta
        probability = 1 / (1 + np.exp(-np.clip(z, -700, 700)))
        fisher = x.T @ ((n * probability * (1 - probability))[:, None] * x)
        try:
            step = np.linalg.solve(fisher, x.T @ (y - n * probability))
        except np.linalg.LinAlgError:
            return {"converged": False, "reason": "Singular information matrix"}
        if np.max(np.abs(step)) < 1e-8:
            converged = True
            break
        old_loss = float(np.sum(n * np.logaddexp(0, z) - y * z))
        scale = 1.
        while scale > 1e-8:
            candidate = beta + scale * step
            new_z = x @ candidate
            if np.sum(n * np.logaddexp(0, new_z) - y * new_z) <= old_loss:
                beta = candidate
                break
            scale *= .5
        else:
            break
    result = {"converged": converged, "coefficients": beta.tolist()}
    if converged:
        se = np.sqrt(np.diag(np.linalg.inv(fisher)))
        result["standard_errors"] = se.tolist()
        result["odds_ratio_per_10_initial_bits"] = math.exp(float(beta[1]))
        result["odds_ratio_per_10_percentage_points_removed"] = math.exp(float(beta[2]))
        result["removal_odds_ratio_95_percent_interval"] = [
            math.exp(float(beta[2] - 1.96 * se[2])), math.exp(float(beta[2] + 1.96 * se[2]))]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    args = parser.parse_args()
    directory = args.run_directory
    data, ranks, geometry = load_records(directory)
    metadata_path = directory / "metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    if not len(data):
        raise ValueError("Run contains no decoded shots")
    if metadata:
        expected = ((data["post_hw"] > 0) & (data["gap_db"] < metadata["gap_threshold_db"])
                    & (metadata["ensemble_size"] > 0) & (data["bypassed"] == 0))
        if np.any(expected != data["ensemble"].astype(bool)):
            raise ValueError("Ensemble flags disagree with recorded gap and configuration")
        if "pdz_removal_threshold" in metadata:
            threshold = metadata["pdz_removal_threshold"]
            expected_bypass = (data["post_hw"] == 0) | ((data["input_hw"] > 0) & (
                data["input_hw"].astype(np.int64) - data["post_hw"] >= threshold * data["input_hw"]))
            if np.any(expected_bypass != data["bypassed"].astype(bool)):
                raise ValueError("Bypass flags disagree with configured PDZ removal threshold")
    if metadata.get("sampler", "monte_carlo") != "monte_carlo":
        raise ValueError("Physical-frequency analysis requires Monte Carlo records, not rare-event samples")

    removed = data["input_hw"].astype(np.int64) - data["post_hw"]
    nonempty = data["input_hw"] > 0
    fraction = removed[nonempty] / data["input_hw"][nonempty]
    # Analyze ensemble associations only where Libra was actually called.
    eligible = nonempty & (data["bypassed"] == 0)
    active = data[eligible]
    fit_fraction = removed[eligible] / active["input_hw"]
    # Bucket -1 is an empty input; 0 is no removal; 1..10 are (0,.1]..(.9,1].
    fraction_bucket = np.full(len(data), -1, dtype=int)
    fraction_bucket[nonempty] = np.ceil(fraction * 10 - 1e-12).astype(int)
    tables = {
        "by_input_hw": group_table(data, [data["input_hw"]], ["input_hw"]),
        "by_post_pdz_hw": group_table(data, [data["post_hw"]], ["post_hw"]),
        "by_removed_bits": group_table(data, [removed], ["removed_bits"]),
        "by_removal_fraction": group_table(data, [fraction_bucket], ["fraction_bucket"]),
        "by_input_hw_and_removed_bits": group_table(data, [data["input_hw"], removed], ["input_hw", "removed_bits"]),
        "by_initial_hw_band_and_fraction": group_table(data, [data["input_hw"] // 5, fraction_bucket], ["input_hw_band", "fraction_bucket"]),
        "by_initial_hw_band_and_removal_band": group_table(
            data, [data["input_hw"] // 10, np.where(fraction_bucket < 0, -1,
                                                  np.maximum(fraction_bucket - 1, 0) // 2)],
            ["input_hw_band", "removal_band"]),
    }
    for name, rows in tables.items():
        write_csv(directory / f"analysis_{name}.csv", rows)

    initial_hw = active["input_hw"].astype(int)
    n_by_hw = np.bincount(initial_hw)
    fraction_mean = np.bincount(initial_hw, weights=fit_fraction) / np.maximum(n_by_hw, 1)
    ensemble_mean = np.bincount(initial_hw, weights=active["ensemble"]) / np.maximum(n_by_hw, 1)
    summary = {
        "shots": len(data), "rank_shots": ranks, "detectors": geometry[0],
        "bypassed_shots": int(data["bypassed"].sum()),
        "ensemble_shots": int(data["ensemble"].sum()),
        "ensemble_frame_changes": int(data["changed"].sum()),
        "pymatching_failures": int(np.sum(data["initial"] != data["truth"])),
        "pipeline_failures": int(np.sum(data["final"] != data["truth"])),
        "ensemble_changes_helped": int(np.sum((data["changed"] != 0) & (data["initial"] != data["truth"]))),
        "ensemble_changes_harmed": int(np.sum((data["changed"] != 0) & (data["final"] != data["truth"]))),
        "mean_input_hw": float(data["input_hw"].mean()), "mean_post_hw": float(data["post_hw"].mean()),
        "mean_fraction_removed_nonempty": float(fraction.mean()) if len(fraction) else None,
        "correlation_removed_bits_and_ensemble_nonempty": correlation(removed[eligible], active["ensemble"]),
        "correlation_removal_fraction_and_ensemble_nonempty": correlation(fit_fraction, active["ensemble"]),
        "correlation_within_exact_initial_hw": correlation(fit_fraction - fraction_mean[initial_hw],
                                                           active["ensemble"] - ensemble_mean[initial_hw]),
        "logistic_controlled_for_initial_hw": controlled_logistic(active),
    }
    summary_path = directory / "analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    lines = ["PDZ + Libra recording analysis", f"Shots: {len(data):,}; ranks: {len(ranks)}",
             f"Libra bypassed: {summary['bypassed_shots']:,} shots (no gap measured)",
             f"Ensemble calls: {summary['ensemble_shots']:,}; ensemble frame changes: {summary['ensemble_frame_changes']:,}",
             f"Initial PyMatching failures: {summary['pymatching_failures']}; pipeline failures: {summary['pipeline_failures']}",
             f"Ensemble changes helped/harmed relative to PyMatching: {summary['ensemble_changes_helped']}/{summary['ensemble_changes_harmed']}",
             "", "Removal fraction (empty input separately; fractions of removed syndrome bits)",
             "Range          Shots     Ensemble    P(ensemble)    Changed    Helped/Harmed"]
    for row in tables["by_removal_fraction"]:
        bucket = row["fraction_bucket"]
        label = "empty input" if bucket == -1 else "0%" if bucket == 0 else f"({10*(bucket-1)},{10*bucket}]%"
        lines.append(f"{label:12} {row['shots']:9,d} {row['ensemble']:10,d} {row['ensemble_probability']:13.6g}"
                     f" {row['changes']:10,d} {row['helped']:5}/{row['harmed']}")
    lines.extend(["", "Post-PDZ Hamming weight", "HW     Shots      Ensemble   P(ensemble)    Changed"])
    for row in tables["by_post_pdz_hw"]:
        if row["ensemble"]:
            lines.append(f"{row['post_hw']:3} {row['shots']:10,d} {row['ensemble']:10,d}"
                         f" {row['ensemble_probability']:13.6g} {row['changes']:10,d}")
    lines.extend(["", "Initial-HW strata and removal fraction",
                  "Initial HW    Removed fraction     Shots    Ensemble     P(ensemble)"])
    for row in tables["by_initial_hw_band_and_removal_band"]:
        if row["removal_band"] < 0:
            continue
        hw, band = row["input_hw_band"], row["removal_band"]
        label = "[0,20]%" if band == 0 else f"({20*band},{20*(band+1)}]%"
        lines.append(f"{10*hw:3}-{10*hw+9:<3}       {label:12} {row['shots']:10,d} {row['ensemble']:10,d}"
                     f" {row['ensemble_probability']:14.6g}")
    lines.extend(["", "Associations (nonempty inputs where Libra was called; bypassed shots excluded):",
                  f"Within exact initial HW correlation: {summary['correlation_within_exact_initial_hw']}",
                  "Logistic model controls for initial HW:",
                  json.dumps(summary["logistic_controlled_for_initial_hw"], indent=2),
                  "", "Caveats: associations do not establish causation or a safe bypass rule.",
                  "Changed means final frame differs from initial PyMatching on an ensemble shot.",
                  "Empty residual syndromes have placeholder gap=0 and do not ensemble.",
                  "Detailed exact-HW and removal-count tables are in analysis_*.csv."])
    report = "\n".join(lines) + "\n"
    (directory / "analysis.txt").write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
