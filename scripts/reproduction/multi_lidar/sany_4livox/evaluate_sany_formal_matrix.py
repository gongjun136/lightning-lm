#!/usr/bin/env python3
"""Evaluate every validated SANY matrix run against the 114 Voxel-SLAM proxy trajectory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from statistics import median
from typing import Any

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_elapsed(value: str) -> float:
    parts = [float(item) for item in value.split(":")]
    if len(parts) == 2:
        return parts[0] * 60.0 + parts[1]
    if len(parts) == 3:
        return parts[0] * 3600.0 + parts[1] * 60.0 + parts[2]
    raise ValueError(f"invalid GNU time elapsed value: {value}")


def parse_gnu_time(path: Path) -> dict[str, float | int | None]:
    text = path.read_text(encoding="utf-8", errors="replace")
    def match(pattern: str) -> str | None:
        result = re.search(pattern, text, flags=re.MULTILINE)
        return result.group(1).strip() if result else None
    elapsed = match(r"^\s*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)$")
    user = match(r"User time \(seconds\):\s*(\S+)$")
    system = match(r"System time \(seconds\):\s*(\S+)$")
    rss = match(r"Maximum resident set size \(kbytes\):\s*(\d+)$")
    return {
        "wall_time_s": parse_elapsed(elapsed) if elapsed else None,
        "user_time_s": float(user) if user else None,
        "system_time_s": float(system) if system else None,
        "cpu_time_s": float(user) + float(system) if user and system else None,
        "max_rss_kib": int(rss) if rss else None,
    }


def load_proxy_first_eight_columns(path: Path) -> np.ndarray:
    rows = []
    with path.open(encoding="utf-8-sig") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_number}: expected at least 8 columns")
            values = [float(value) for value in fields[:8]]
            if not all(np.isfinite(values)):
                raise ValueError(f"{path}:{line_number}: non-finite proxy pose")
            rows.append(values)
    result = np.asarray(rows, dtype=np.float64)
    if len(result) < 3 or np.any(np.diff(result[:, 0]) <= 0.0):
        raise ValueError("proxy timestamps are not strictly increasing")
    norms = np.linalg.norm(result[:, 4:8], axis=1)
    if np.any(norms < 1e-12) or np.any(np.abs(norms - 1.0) > 1e-3):
        raise ValueError("proxy contains invalid quaternions")
    result[:, 4:8] /= norms[:, None]
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-root", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--m3dgr-evaluator-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.m3dgr_evaluator_dir.resolve()))
    import evaluate_m3dgr_runs as evaluator

    reference = load_proxy_first_eight_columns(args.reference)
    rows: list[dict[str, Any]] = []
    for cell_dir in sorted((args.experiment_root / "runs").glob("*/repeat_[0-9][0-9]")):
        selection_path = cell_dir / "selected_attempt.json"
        if not selection_path.is_file():
            continue
        selection = json.loads(selection_path.read_text(encoding="utf-8"))
        attempt = Path(selection["selected_attempt"])
        manifest = json.loads((attempt / "completion_manifest.json").read_text(encoding="utf-8"))
        if not manifest["output_validation"]["passed"]:
            continue
        trajectory_path = attempt / "results" / "trajectory_lidar114.tum"
        estimate = evaluator.load_estimated_tum(trajectory_path)
        metrics = evaluator.evaluate_full_trajectory(estimate, reference, 0.20, True)
        resource = parse_gnu_time(attempt / "gnu_time.txt")
        row = {
            "variant": manifest["variant"],
            "repeat": int(manifest["repeat"]),
            "attempt": int(manifest["attempt"]),
            "trajectory": str(trajectory_path),
            "trajectory_sha256": sha256(trajectory_path),
            "trajectory_pose_count": len(estimate),
            "associated_pose_count": metrics["association"]["associated_pose_count"],
            "associated_pose_ratio": metrics["association"]["associated_pose_ratio"],
            "ape_translation_rmse_m": metrics["ape"]["translation_m"]["rmse"],
            "ape_translation_median_m": metrics["ape"]["translation_m"]["median"],
            "ape_translation_p95_m": metrics["ape"]["translation_m"]["p95"],
            "ape_rotation_rmse_deg": metrics["ape"]["rotation_deg"]["rmse"],
            "ape_rotation_median_deg": metrics["ape"]["rotation_deg"]["median"],
            "ape_rotation_p95_deg": metrics["ape"]["rotation_deg"]["p95"],
            "map_point_count": manifest["output_validation"]["map_point_count"],
            "frame_contract_ratio": manifest["output_validation"]["frame_contract_ratio"],
            **resource,
        }
        rows.append(row)
    if len(rows) != 21:
        raise SystemExit(f"expected 21 validated cells, found {len(rows)}")

    variants: dict[str, Any] = {}
    for variant in sorted({row["variant"] for row in rows}):
        selected = [row for row in rows if row["variant"] == variant]
        if sorted(row["repeat"] for row in selected) != [1, 2, 3]:
            raise SystemExit(f"{variant}: repeats are incomplete")
        variants[variant] = {
            "runs": 3,
            "unique_trajectory_sha256_count": len({row["trajectory_sha256"] for row in selected}),
            "ape_translation_rmse_m_median": median(row["ape_translation_rmse_m"] for row in selected),
            "ape_translation_rmse_m_range": [
                min(row["ape_translation_rmse_m"] for row in selected),
                max(row["ape_translation_rmse_m"] for row in selected),
            ],
            "ape_rotation_rmse_deg_median": median(row["ape_rotation_rmse_deg"] for row in selected),
            "map_point_count_median": median(row["map_point_count"] for row in selected),
            "wall_time_s_median": median(row["wall_time_s"] for row in selected),
            "cpu_time_s_median": median(row["cpu_time_s"] for row in selected),
            "max_rss_kib_median": median(row["max_rss_kib"] for row in selected),
        }
    result = {
        "reference": str(args.reference),
        "reference_sha256": sha256(args.reference),
        "reference_role": "114同源 Voxel-SLAM 代理轨迹，仅用于相对一致性比较，不是真值或绝对精度依据",
        "alignment": "fixed-scale rigid SE(3), no Sim(3), estimate interpolation at proxy timestamps",
        "evaluator": str((args.m3dgr_evaluator_dir / "evaluate_m3dgr_runs.py").resolve()),
        "evaluator_sha256": sha256(args.m3dgr_evaluator_dir / "evaluate_m3dgr_runs.py"),
        "runs": rows,
        "variants": variants,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "sany_proxy_trajectory_evaluation.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    with (args.output_dir / "sany_proxy_trajectory_runs.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"runs": len(rows), "variants": variants}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
