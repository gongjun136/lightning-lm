#!/usr/bin/env python3
"""Summarize one offline localization run and optionally compare against a TUM reference."""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from pathlib import Path
from statistics import mean, median


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--localization-csv", required=True, type=Path)
    parser.add_argument("--trajectory-tum", required=True, type=Path)
    parser.add_argument("--lidar-loc-tum", type=Path)
    parser.add_argument("--reference-tum", type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--errors-csv", required=True, type=Path)
    parser.add_argument("--max-association-dt", type=float, default=0.05)
    parser.add_argument("--max-speed", type=float, default=0.0)
    parser.add_argument("--max-z-range", type=float, default=0.0)
    parser.add_argument("--max-axis-range", type=float, default=0.0)
    return parser.parse_args()


def finite_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def percentile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    pos = (len(ordered) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - pos) + ordered[hi] * (pos - lo)


def stats(values: list[float]) -> dict[str, float | int | None]:
    values = [v for v in values if math.isfinite(v)]
    if not values:
        return {"count": 0, "min": None, "mean": None, "median": None, "p95": None, "max": None}
    return {
        "count": len(values),
        "min": min(values),
        "mean": mean(values),
        "median": median(values),
        "p95": percentile(values, 0.95),
        "max": max(values),
    }


def read_tum(path: Path) -> list[tuple[float, tuple[float, float, float]]]:
    poses: list[tuple[float, tuple[float, float, float]]] = []
    if not path or not path.exists():
        return poses
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 8:
                continue
            values = [finite_float(part) for part in parts[:4]]
            if any(value is None for value in values):
                continue
            poses.append((values[0], (values[1], values[2], values[3])))  # type: ignore[arg-type]
    poses.sort(key=lambda item: item[0])
    return poses


def associate_errors(
    estimate: list[tuple[float, tuple[float, float, float]]],
    reference: list[tuple[float, tuple[float, float, float]]],
    max_dt: float,
) -> list[dict[str, float]]:
    if not estimate or not reference:
        return []
    ref_times = [item[0] for item in reference]
    rows: list[dict[str, float]] = []
    for timestamp, pos in estimate:
        idx = bisect.bisect_left(ref_times, timestamp)
        candidates = []
        if idx < len(reference):
            candidates.append(reference[idx])
        if idx > 0:
            candidates.append(reference[idx - 1])
        if not candidates:
            continue
        ref_time, ref_pos = min(candidates, key=lambda item: abs(item[0] - timestamp))
        dt = timestamp - ref_time
        if abs(dt) > max_dt:
            continue
        dx = pos[0] - ref_pos[0]
        dy = pos[1] - ref_pos[1]
        dz = pos[2] - ref_pos[2]
        rows.append(
            {
                "timestamp": timestamp,
                "reference_timestamp": ref_time,
                "dt": dt,
                "dx": dx,
                "dy": dy,
                "dz": dz,
                "error_3d": math.sqrt(dx * dx + dy * dy + dz * dz),
                "error_xy": math.sqrt(dx * dx + dy * dy),
            }
        )
    return rows


def summarize_motion(
    poses: list[tuple[float, tuple[float, float, float]]],
    max_speed: float,
    max_z_range: float,
    max_axis_range: float,
) -> dict[str, object]:
    if not poses:
        return {"poses": 0, "passed": False, "failure_reasons": ["empty_trajectory"]}

    xs = [pose[1][0] for pose in poses]
    ys = [pose[1][1] for pose in poses]
    zs = [pose[1][2] for pose in poses]
    steps: list[float] = []
    speeds: list[float] = []
    for previous, current in zip(poses, poses[1:]):
        dt = current[0] - previous[0]
        delta = math.sqrt(sum((a - b) ** 2 for a, b in zip(current[1], previous[1])))
        steps.append(delta)
        if dt > 0.0:
            speeds.append(delta / dt)

    start_end_distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(poses[-1][1], poses[0][1])))
    x_range = max(xs) - min(xs)
    y_range = max(ys) - min(ys)
    z_range = max(zs) - min(zs)
    speed_violations = sum(1 for speed in speeds if max_speed > 0.0 and speed > max_speed)
    reasons: list[str] = []
    if max_speed > 0.0 and speed_violations:
        reasons.append("max_speed")
    if max_z_range > 0.0 and z_range > max_z_range:
        reasons.append("z_range")
    if max_axis_range > 0.0 and (x_range > max_axis_range or y_range > max_axis_range):
        reasons.append("xy_axis_range")

    return {
        "poses": len(poses),
        "duration_s": poses[-1][0] - poses[0][0],
        "start_end_distance_m": start_end_distance,
        "x_range_m": x_range,
        "y_range_m": y_range,
        "z_range_m": z_range,
        "step_m": stats(steps),
        "speed_mps": stats(speeds),
        "speed_violation_count": speed_violations,
        "passed": not reasons,
        "failure_reasons": reasons,
    }


def write_errors(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["timestamp", "reference_timestamp", "dt", "dx", "dy", "dz", "error_3d", "error_xy"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summarize_localization_csv(path: Path) -> dict[str, object]:
    rows = []
    with path.open(newline="", encoding="utf-8", errors="replace") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            rows.append(row)
    total = len(rows)
    valid = sum(1 for row in rows if row.get("lidar_loc_valid") == "1")
    status_counts: dict[str, int] = {}
    for row in rows:
        status = row.get("status", "")
        status_counts[status] = status_counts.get(status, 0) + 1
    confidences = [v for row in rows if (v := finite_float(row.get("confidence"))) is not None]
    processing_ms = [v for row in rows if (v := finite_float(row.get("processing_ms"))) is not None]
    active_chunks = [v for row in rows if (v := finite_float(row.get("active_map_chunks"))) is not None]
    iterations = [v for row in rows if (v := finite_float(row.get("match_iterations"))) is not None]
    odom_delta = [v for row in rows if (v := finite_float(row.get("loc_odom_delta"))) is not None]
    relocalization_attempts = sum(1 for row in rows if row.get("relocalization_attempted") == "1")
    relocalization_candidates = sum(1 for row in rows if row.get("relocalization_candidate_found") == "1")
    relocalization_accepts = sum(1 for row in rows if row.get("relocalization_accepted") == "1")
    return {
        "frames": total,
        "valid_frames": valid,
        "valid_ratio": valid / total if total else None,
        "status_counts": status_counts,
        "confidence": stats(confidences),
        "processing_ms": stats(processing_ms),
        "active_map_chunks": stats(active_chunks),
        "match_iterations": stats(iterations),
        "loc_odom_delta": stats(odom_delta),
        "relocalization": {
            "attempts": relocalization_attempts,
            "candidates": relocalization_candidates,
            "accepts": relocalization_accepts,
        },
    }


def main() -> int:
    args = parse_args()
    localization_summary = summarize_localization_csv(args.localization_csv)
    estimate = read_tum(args.trajectory_tum)
    lidar_loc_estimate = read_tum(args.lidar_loc_tum) if args.lidar_loc_tum else []
    reference = read_tum(args.reference_tum) if args.reference_tum else []
    errors = associate_errors(estimate, reference, args.max_association_dt)
    write_errors(args.errors_csv, errors)
    error_3d = [row["error_3d"] for row in errors]
    error_xy = [row["error_xy"] for row in errors]
    payload = {
        "schema_version": 2,
        "localization_csv": str(args.localization_csv.resolve()),
        "trajectory_tum": str(args.trajectory_tum.resolve()),
        "lidar_loc_tum": str(args.lidar_loc_tum.resolve()) if args.lidar_loc_tum else None,
        "reference_tum": str(args.reference_tum.resolve()) if args.reference_tum else None,
        "max_association_dt": args.max_association_dt,
        "localization": localization_summary,
        "trajectory": {
            "estimate_poses": len(estimate),
            "reference_poses": len(reference),
            "associated_pairs": len(errors),
            "error_3d": stats(error_3d),
            "error_xy": stats(error_xy),
            "motion": summarize_motion(
                estimate, args.max_speed, args.max_z_range, args.max_axis_range
            ),
        },
        "lidar_loc_trajectory": {
            "estimate_poses": len(lidar_loc_estimate),
            "motion": summarize_motion(
                lidar_loc_estimate, args.max_speed, args.max_z_range, args.max_axis_range
            ),
        },
        "physical_thresholds": {
            "max_speed_mps": args.max_speed,
            "max_z_range_m": args.max_z_range,
            "max_x_or_y_range_m": args.max_axis_range,
            "start_end_distance_is_diagnostic_only": True,
        },
    }
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.summary_json.with_name(args.summary_json.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(args.summary_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
