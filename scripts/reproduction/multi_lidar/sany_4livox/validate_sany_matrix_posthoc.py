#!/usr/bin/env python3
"""Independent strict post-hoc audit for all selected SANY matrix attempts."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_file_record(record: dict[str, Any], errors: list[str]) -> None:
    path = Path(record["path"])
    if not path.is_file():
        errors.append(f"frozen file is missing: {path}")
        return
    if "size" in record and path.stat().st_size != int(record["size"]):
        errors.append(f"frozen file size changed: {path}")
    if "size_bytes" in record and path.stat().st_size != int(record["size_bytes"]):
        errors.append(f"frozen file size changed: {path}")
    if sha256(path).lower() != str(record["sha256"]).lower():
        errors.append(f"frozen file SHA changed: {path}")


def trajectory_time_contract(times: list[float], label: str, errors: list[str]) -> dict[str, Any]:
    array = np.asarray(times, dtype=np.float64)
    intervals = np.diff(array)
    duration = float(array[-1] - array[0]) if len(array) >= 2 else 0.0
    rate = float((len(array) - 1) / duration) if duration > 0.0 else 0.0
    maximum_gap = float(np.max(intervals)) if len(intervals) else math.inf
    if duration < 105.0 or not 9.5 <= rate <= 10.5 or maximum_gap > 0.20:
        errors.append(f"{label}: trajectory duration/rate/gap contract failed")
    return {"count": len(array), "duration_s": duration, "average_rate_hz": rate, "maximum_gap_s": maximum_gap}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-root", type=Path, required=True)
    parser.add_argument("--tool-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.tool_dir.resolve()))
    import run_sany_formal_matrix as matrix
    import analyze_sany_map_cross_source as map_qc

    fingerprint_path = args.experiment_root / "formal_matrix" / "fingerprint.json"
    fingerprint = json.loads(fingerprint_path.read_text(encoding="utf-8"))
    expected_fingerprint = fingerprint["experiment_fingerprint"]
    errors: list[str] = []
    canonical = dict(fingerprint)
    canonical.pop("experiment_fingerprint")
    recomputed = hashlib.sha256(
        json.dumps(canonical, ensure_ascii=False, sort_keys=True).encode("utf-8")
    ).hexdigest()
    if recomputed != expected_fingerprint:
        errors.append("fingerprint JSON canonical digest is inconsistent")
    check_file_record(fingerprint["base_config"], errors)
    for record in fingerprint["configs"]:
        check_file_record(record, errors)
    for key in ("binary", "template", "controller", "install_setup"):
        check_file_record(fingerprint[key], errors)
    for record in fingerprint["bag"]["files"]:
        check_file_record(record, errors)

    run_audits = []
    for cell_dir in sorted((args.experiment_root / "runs").glob("*/repeat_[0-9][0-9]")):
        selection_path = cell_dir / "selected_attempt.json"
        if not selection_path.is_file():
            continue
        selection = json.loads(selection_path.read_text(encoding="utf-8"))
        attempt = Path(selection["selected_attempt"])
        completion = json.loads((attempt / "completion_manifest.json").read_text(encoding="utf-8"))
        variant = completion["variant"]
        repeat = int(completion["repeat"])
        run_errors: list[str] = []
        if completion["experiment_fingerprint"] != expected_fingerprint:
            run_errors.append("completion fingerprint mismatch")

        lidar_path = attempt / "results" / "trajectory_lidar114.tum"
        rear_path = attempt / "results" / "trajectory_rear_axle.tum"
        lidar_times, lidar_errors = matrix.read_tum(lidar_path)
        rear_times, rear_errors = matrix.read_tum(rear_path)
        run_errors.extend(lidar_errors)
        run_errors.extend(rear_errors)
        lidar_time = trajectory_time_contract(lidar_times, "lidar114", run_errors)
        rear_time = trajectory_time_contract(rear_times, "rear_axle", run_errors)
        if len(lidar_times) != len(rear_times) or np.max(np.abs(np.asarray(lidar_times) - np.asarray(rear_times))) > 1e-6:
            run_errors.append("lidar/rear trajectory timestamps differ")

        frame_path = attempt / "results" / "frame_stats.csv"
        import csv
        with frame_path.open(encoding="utf-8", newline="") as stream:
            frame_rows = list(csv.DictReader(stream))
        frame_audit: dict[str, Any] = {"count": len(frame_rows)}
        if variant.startswith("C0_"):
            if frame_rows:
                run_errors.append("single-lidar cell has multi-lidar frame rows")
        else:
            begins = np.asarray([float(row["begin_time"]) for row in frame_rows])
            ends = np.asarray([float(row["end_time"]) for row in frame_rows])
            end_intervals = np.diff(ends)
            duration = float(ends[-1] - ends[0]) if len(ends) >= 2 else 0.0
            rate = float((len(ends) - 1) / duration) if duration > 0.0 else 0.0
            conservation_failures = 0
            for row in frame_rows:
                source_sum = sum(int(row[f"points_lidar_{lidar_id}"]) for lidar_id in range(4))
                conservation_failures += int(source_sum != int(row["merged_points"]))
            frame_audit.update({
                "duration_s": duration, "average_rate_hz": rate,
                "maximum_end_gap_s": float(np.max(end_intervals)) if len(end_intervals) else None,
                "point_conservation_failure_count": conservation_failures,
            })
            if (
                len(frame_rows) < 1100 or np.any(np.diff(begins) <= 0.0) or np.any(end_intervals <= 0.0)
                or duration < 115.0 or not 9.5 <= rate <= 10.5
                or (len(end_intervals) and np.max(end_intervals) > 0.20) or conservation_failures
            ):
                run_errors.append("frame_stats duration/rate/monotonicity/point-conservation contract failed")

        map_path = attempt / "results" / "map_lio.pcd"
        arrays = map_qc.load_binary_compressed_pcd(map_path)
        required_fields = {"x", "y", "z", "intensity", "time", "lidar_id"}
        if set(arrays) != required_fields:
            run_errors.append(f"PCD fields differ from public six-field schema: {sorted(arrays)}")
        xyz = np.column_stack((arrays["x"], arrays["y"], arrays["z"]))
        intensity = arrays["intensity"]
        point_time = arrays["time"]
        ids = arrays["lidar_id"].astype(np.int64)
        nonfinite = int(np.count_nonzero(~np.all(np.isfinite(xyz), axis=1)))
        nonfinite += int(np.count_nonzero(~np.isfinite(intensity))) + int(np.count_nonzero(~np.isfinite(point_time)))
        invalid_ids = int(np.count_nonzero((ids < 0) | (ids > 3)))
        dropped = next((item[2] for item in matrix.VARIANTS if item[0] == variant), None)
        expected_ids = {0} if variant.startswith("C0_") else ({0, 1, 2, 3} - ({dropped} if dropped is not None else set()))
        actual_ids = {int(value) for value in np.unique(ids)}
        counts_by_id = {str(lidar_id): int(np.count_nonzero(ids == lidar_id)) for lidar_id in range(4)}
        if nonfinite or invalid_ids or actual_ids != expected_ids or any(counts_by_id[str(lidar_id)] < 10_000 for lidar_id in expected_ids):
            run_errors.append("PCD finite/source-ID/content contract failed")
        map_audit = {
            "point_count": len(ids), "nonfinite_field_count": nonfinite,
            "invalid_lidar_id_count": invalid_ids, "actual_lidar_ids": sorted(actual_ids),
            "points_by_lidar": counts_by_id, "sha256": sha256(map_path),
        }

        evidence_files = [
            attempt / "gnu_time.txt", attempt / "controller.stdout.log", attempt / "controller.stderr.log",
            attempt / "logs" / "run_frontend_offline.stderr.log",
        ]
        evidence = []
        for path in evidence_files:
            if not path.is_file():
                run_errors.append(f"resource/log evidence is missing: {path.name}")
            else:
                evidence.append({"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)})
        audit = {
            "variant": variant, "repeat": repeat, "attempt": str(attempt),
            "passed": not run_errors, "errors": run_errors,
            "trajectory_lidar114": lidar_time, "trajectory_rear_axle": rear_time,
            "frame_stats": frame_audit, "map": map_audit, "evidence_files": evidence,
        }
        (attempt / "posthoc_validation.json").write_text(
            json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        run_audits.append(audit)
        errors.extend(f"{variant}/repeat_{repeat:02d}: {error}" for error in run_errors)

    if len(run_audits) != 21:
        errors.append(f"expected 21 selected attempts, found {len(run_audits)}")
    result = {
        "created_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "validator": str(Path(__file__).resolve()), "validator_sha256": sha256(Path(__file__).resolve()),
        "experiment_fingerprint": expected_fingerprint,
        "passed": not errors, "errors": errors, "run_count": len(run_audits), "runs": run_audits,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": not errors, "runs": len(run_audits), "errors": errors[:10]}, ensure_ascii=False))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
