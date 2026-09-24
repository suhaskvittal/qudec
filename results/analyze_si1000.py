#!/usr/bin/env python3
"""OpenAI GPT-6-Sol, 24 September 2026: validate and summarize captured SI1000 shots."""

import json
import statistics
import sys
from collections import Counter
from pathlib import Path


def main(path: Path) -> None:
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    assert rows[0]["type"] == "metadata"
    assert rows[-1]["type"] == "summary"
    meta, summary = rows[0], rows[-1]
    shots = rows[1:-1]
    assert all(row["type"] == "shot" for row in shots)
    assert summary["completed_shots"] == meta["requested_shots"]
    assert len(shots) == summary["only_pymatching_errors"]
    assert summary["pymatching_errors"] == (summary["both_errors"] + len(shots) +
                                             summary.get("uncaptured_raw_only_pymatching", 0))
    assert summary["tesseract_errors"] == summary["both_errors"] + summary["only_tesseract_errors"]
    if "paper_scored_tesseract_errors" in summary:
        assert summary["paper_scored_tesseract_errors"] >= summary["tesseract_errors"]
        assert summary["paper_scored_only_pymatching_errors"] == len(shots)
        assert all(not row["tesseract_low_confidence"] for row in shots)
    assert len({row["shot_index"] for row in shots}) == len(shots)
    assert all(0 <= row["shot_index"] < summary["completed_shots"] for row in shots)
    assert all(row["pymatching_prediction"] != row["truth"] and
               row["tesseract_prediction"] == row["truth"] for row in shots)
    assert all(row["detector_ids"] == sorted(set(row["detector_ids"])) and
               all(0 <= det < meta["detector_count"] for det in row["detector_ids"])
               for row in shots)
    assert sum(row["tesseract_low_confidence"] for row in shots) == summary["captured_low_confidence"]
    assert sum(len(row["detector_ids"]) for row in shots) == summary["captured_detector_total"]
    assert sum(row["opposite_basis_detection_count"] for row in shots) == summary["captured_opposite_detector_total"]

    # OpenAI GPT-6-Sol: Check the X-check ID windows implied by the saved
    # head/body/tail circuit and report empirical event patterns.
    first_body = (meta["distance"] ** 2 - 1) // 2
    per_body = meta["distance"] ** 2 - 1
    opposite_per_body = per_body // 2
    basis_counts = Counter()
    round_counts = Counter()
    for shot in shots:
        opposite = 0
        for det in shot["detector_ids"]:
            if det < first_body:
                basis_counts["first_Z"] += 1
            elif det < first_body + (meta["rounds"] - 1) * per_body:
                body_offset = det - first_body
                round_index = body_offset // per_body + 2
                basis = "X" if body_offset % per_body >= opposite_per_body else "Z"
                basis_counts[f"body_{basis}"] += 1
                round_counts[round_index] += 1
                opposite += basis == "X"
            else:
                basis_counts["final_Z"] += 1
        assert opposite == shot["opposite_basis_detection_count"]

    sizes = [len(row["detector_ids"]) for row in shots]
    print(json.dumps({
        "completed_shots": summary["completed_shots"],
        "pymatching_errors": summary["pymatching_errors"],
        "tesseract_errors": summary["tesseract_errors"],
        "only_pymatching_errors": len(shots),
        "only_tesseract_errors": summary["only_tesseract_errors"],
        "both_errors": summary["both_errors"],
        "tesseract_low_confidence_all": summary["tesseract_low_confidence"],
        "tesseract_low_confidence_captured": summary["captured_low_confidence"],
        "truth_counts_captured": dict(Counter(row["truth"] for row in shots)),
        "detectors_per_captured_shot": {
            "min": min(sizes) if sizes else None,
            "median": statistics.median(sizes) if sizes else None,
            "max": max(sizes) if sizes else None,
            "mean": statistics.mean(sizes) if sizes else None,
        },
        "fired_detector_counts_by_basis": dict(basis_counts),
        "fired_detector_counts_by_body_round": dict(sorted(round_counts.items())),
    }, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
