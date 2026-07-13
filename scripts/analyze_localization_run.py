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
    parser.add_argument("--reference-tum", type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--errors-csv", required=True, type=Path)
    parser.add_argument("--max-association-dt", type=float, default=0.05)
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
    }


def main() -> int:
    args = parse_args()
    localization_summary = summarize_localization_csv(args.localization_csv)
    estimate = read_tum(args.trajectory_tum)
    reference = read_tum(args.reference_tum) if args.reference_tum else []
    errors = associate_errors(estimate, reference, args.max_association_dt)
    write_errors(args.errors_csv, errors)
    error_3d = [row["error_3d"] for row in errors]
    error_xy = [row["error_xy"] for row in errors]
    payload = {
        "schema_version": 1,
        "localization_csv": str(args.localization_csv.resolve()),
        "trajectory_tum": str(args.trajectory_tum.resolve()),
        "reference_tum": str(args.reference_tum.resolve()) if args.reference_tum else None,
        "max_association_dt": args.max_association_dt,
        "localization": localization_summary,
        "trajectory": {
            "estimate_poses": len(estimate),
            "reference_poses": len(reference),
            "associated_pairs": len(errors),
            "error_3d": stats(error_3d),
            "error_xy": stats(error_xy),
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
