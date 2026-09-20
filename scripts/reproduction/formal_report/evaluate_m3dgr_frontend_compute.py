#!/usr/bin/env python3
"""Evaluate standardized per-frame compute timing for the M3DGR LIO matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Iterable


SEQUENCES = ("Grass02", "Outdoor04", "Z-Rough-Road01", "Dark01")
METHODS = ("lightning_lm", "fastlio", "fastlivo2_lio", "voxel_slam_frontend")
LIGHTNING_STAGES = (
    "preprocess_ms",
    "imu_undistort_ms",
    "downsample_ms",
    "match_setup_ms",
    "scan_match_ms",
    "map_update_ms",
)
DEADLINE_MS = 100.0


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "m3dgr_frontend_compute",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "analysis" / "m3dgr_frontend_compute",
    )
    parser.add_argument("--sequences", default=",".join(SEQUENCES))
    parser.add_argument("--methods", default=",".join(METHODS))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmup-frames", type=int, default=10)
    parser.add_argument("--deadline-ms", type=float, default=DEADLINE_MS)
    return parser.parse_args()


def parse_metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    return values


def parse_benchmark_log(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        marker = "LIO_BENCH_FRAME "
        if marker not in line:
            continue
        fields: dict[str, object] = {"source_line": line_number}
        for token in line.split(marker, 1)[1].strip().split():
            key, separator, value = token.partition("=")
            if not separator:
                continue
            if key in {"method", "phase"}:
                fields[key] = value
            elif key in {"input_points", "output_points"}:
                fields[key] = int(value)
            else:
                fields[key] = float(value)
        rows.append(fields)
    return rows


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def describe(values: Iterable[float], prefix: str) -> dict[str, float | int]:
    data = [float(value) for value in values]
    if not data:
        raise ValueError(f"cannot describe empty values for {prefix}")
    median = statistics.median(data)
    return {
        f"{prefix}_count": len(data),
        f"{prefix}_mean": statistics.fmean(data),
        f"{prefix}_std": statistics.stdev(data) if len(data) > 1 else 0.0,
        f"{prefix}_median": median,
        f"{prefix}_p95": percentile(data, 0.95),
        f"{prefix}_p99": percentile(data, 0.99),
        f"{prefix}_max": max(data),
        f"{prefix}_mad": statistics.median(abs(value - median) for value in data),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def log_path(run_dir: Path, method: str) -> Path:
    name = "algorithm.stderr.log" if method == "lightning_lm" else "algorithm.log"
    return run_dir / "logs" / name


def run_contract(metadata: dict[str, str], method: str) -> bool:
    if method == "lightning_lm":
        return (
            metadata.get("completion") == "reached_final_lidar"
            and metadata.get("algorithm_rc") == "0"
            and metadata.get("watchdog_status") == "completed"
            and metadata.get("invalid_count") == "0"
            and metadata.get("nonmonotonic_count") == "0"
        )
    return (
        metadata.get("rosbag_play_rc") == "0"
        and metadata.get("algorithm_contract") == "passed"
        and metadata.get("subscriber_contract") == "passed"
        and metadata.get("algorithm_launch_rc") in {"0", "137"}
    )


def main() -> int:
    args = parse_args()
    sequences = tuple(item.strip() for item in args.sequences.split(",") if item.strip())
    methods = tuple(item.strip() for item in args.methods.split(",") if item.strip())
    if set(sequences) - set(SEQUENCES) or set(methods) - set(METHODS):
        raise SystemExit("unknown sequence or method")
    if args.repeats < 1 or args.warmup_frames < 0 or args.deadline_ms <= 0:
        raise SystemExit("repeats/deadline must be positive and warmup non-negative")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []
    warnings: list[str] = []
    frame_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []
    included_by_run: dict[tuple[str, str, int], list[dict[str, object]]] = {}

    for sequence in sequences:
        for method in methods:
            for repeat in range(1, args.repeats + 1):
                run_dir = args.runs_root / sequence / method / f"repeat_{repeat:02d}"
                metadata_path = run_dir / "run_metadata.txt"
                timing_log = log_path(run_dir, method)
                run_key = f"{sequence}/{method}/repeat_{repeat:02d}"
                if not metadata_path.is_file() or not timing_log.is_file():
                    failures.append(f"missing run artifacts: {run_key}")
                    continue
                metadata = parse_metadata(metadata_path)
                if not run_contract(metadata, method):
                    failures.append(f"failed run contract: {run_key}")
                if method != "lightning_lm" and (
                    metadata.get("completion") != "reached_final_lidar"
                    or metadata.get("algorithm_launch_rc") != "0"
                    or any(
                        int(metadata.get(field, "0")) > 0
                        for field in (
                            "invalid_count", "nonmonotonic_count", "excessive_output_gap_count",
                            "shutdown_abnormal_exit_count",
                        )
                    )
                ):
                    warnings.append(
                        f"ROS trajectory-output audit issues (outside compute boundary): {run_key}; "
                        f"completion={metadata.get('completion', 'unknown')}, "
                        f"launch_rc={metadata.get('algorithm_launch_rc', 'unknown')}, "
                        f"invalid={metadata.get('invalid_count', '0')}, "
                        f"nonmonotonic={metadata.get('nonmonotonic_count', '0')}, "
                        f"excessive_gaps={metadata.get('excessive_output_gap_count', '0')}, "
                        f"shutdown_abnormal={metadata.get('shutdown_abnormal_exit_count', '0')}"
                    )
                try:
                    records = parse_benchmark_log(timing_log)
                except (ValueError, OverflowError) as exc:
                    failures.append(f"invalid timing log {run_key}: {exc}")
                    continue
                tracking = [row for row in records if row.get("phase") == "tracking"]
                if not tracking:
                    failures.append(f"no tracking timing rows: {run_key}")
                    continue

                timestamps = [float(row["timestamp_s"]) for row in tracking if "timestamp_s" in row]
                if len(timestamps) != len(tracking) or any(b <= a for a, b in zip(timestamps, timestamps[1:])):
                    failures.append(f"non-monotonic or missing timing timestamps: {run_key}")

                valid_rows: list[dict[str, object]] = []
                tracking_seen = 0
                for frame_index, record in enumerate(records):
                    record_method = record.get("method")
                    if record_method != method:
                        failures.append(f"method mismatch at {run_key}:{record.get('source_line')}")
                        continue
                    required = ("timestamp_s", "preprocess_ms", "core_update_ms", "total_ms", "input_points", "output_points")
                    if any(field not in record for field in required):
                        failures.append(f"missing timing field at {run_key}:{record.get('source_line')}")
                        continue
                    numeric_values = [
                        float(value) for key, value in record.items()
                        if key.endswith("_ms") or key == "timestamp_s"
                    ]
                    if not numeric_values or not all(math.isfinite(value) for value in numeric_values):
                        failures.append(f"non-finite timing value at {run_key}:{record.get('source_line')}")
                        continue
                    if float(record["total_ms"]) <= 0 or float(record["preprocess_ms"]) < 0 or float(record["core_update_ms"]) < 0:
                        failures.append(f"out-of-range timing value at {run_key}:{record.get('source_line')}")
                        continue
                    identity_error = abs(
                        float(record["total_ms"])
                        - float(record["preprocess_ms"])
                        - float(record["core_update_ms"])
                    )
                    if identity_error > 0.05:
                        failures.append(f"total identity mismatch ({identity_error:.6f} ms): {run_key}:{record.get('source_line')}")
                        continue
                    if method == "lightning_lm" and record.get("phase") == "tracking":
                        if any(stage not in record for stage in LIGHTNING_STAGES):
                            failures.append(f"missing Lightning stage at {run_key}:{record.get('source_line')}")
                            continue
                        stage_error = abs(
                            float(record["core_update_ms"])
                            - sum(float(record[stage]) for stage in LIGHTNING_STAGES[1:])
                        )
                        if stage_error > 0.05:
                            failures.append(f"Lightning stage identity mismatch ({stage_error:.6f} ms): {run_key}")
                            continue

                    is_tracking = record.get("phase") == "tracking"
                    tracking_index = tracking_seen
                    if is_tracking:
                        tracking_seen += 1
                    included = is_tracking and tracking_index >= args.warmup_frames
                    row = {
                        "sequence": sequence,
                        "method": method,
                        "repeat": repeat,
                        "frame_index": frame_index,
                        "tracking_index": tracking_index if is_tracking else None,
                        "is_warmup": is_tracking and not included,
                        "included": included,
                        **record,
                    }
                    frame_rows.append(row)
                    if included:
                        valid_rows.append(row)

                if not valid_rows:
                    failures.append(f"no included rows after warmup: {run_key}")
                    continue
                included_by_run[(sequence, method, repeat)] = valid_rows
                trajectory_lines = int(metadata.get("trajectory_lines", "0"))
                coverage = len(tracking) / trajectory_lines if trajectory_lines else 0.0
                if coverage < 0.8:
                    failures.append(f"timing coverage below 80% ({coverage:.3f}): {run_key}")
                elif coverage < 0.95:
                    warnings.append(f"timing coverage below 95% ({coverage:.3f}): {run_key}")

                totals = [float(row["total_ms"]) for row in valid_rows]
                preprocess = [float(row["preprocess_ms"]) for row in valid_rows]
                core = [float(row["core_update_ms"]) for row in valid_rows]
                total_stats = describe(totals, "total_ms")
                mean_total = float(total_stats["total_ms_mean"])
                run_rows.append({
                    "sequence": sequence,
                    "method": method,
                    "repeat": repeat,
                    "raw_record_count": len(records),
                    "tracking_record_count": len(tracking),
                    "included_record_count": len(valid_rows),
                    "trajectory_lines": trajectory_lines,
                    "timing_coverage": coverage,
                    "warmup_frames": args.warmup_frames,
                    **total_stats,
                    **describe(preprocess, "preprocess_ms"),
                    **describe(core, "core_update_ms"),
                    "deadline_ms": args.deadline_ms,
                    "deadline_miss_count": sum(value > args.deadline_ms for value in totals),
                    "deadline_miss_rate": sum(value > args.deadline_ms for value in totals) / len(totals),
                    "compute_utilization": mean_total / args.deadline_ms,
                    "effective_throughput_hz": 1000.0 / mean_total,
                })

    expected_runs = len(sequences) * len(methods) * args.repeats
    if len(run_rows) != expected_runs:
        failures.append(f"incomplete run summary: expected {expected_runs}, got {len(run_rows)}")

    for row in run_rows:
        peers = [
            item for item in run_rows
            if item["sequence"] == row["sequence"] and item["method"] == row["method"]
        ]
        expected_tracking_frames = max(int(item["tracking_record_count"]) for item in peers)
        repeat_frame_coverage = int(row["tracking_record_count"]) / expected_tracking_frames
        row["expected_tracking_frames_from_repeats"] = expected_tracking_frames
        row["repeat_frame_coverage"] = repeat_frame_coverage
        run_key = f"{row['sequence']}/{row['method']}/repeat_{int(row['repeat']):02d}"
        if repeat_frame_coverage < 0.98:
            failures.append(
                f"timing frames below 98% of same-cell maximum ({repeat_frame_coverage:.3f}): {run_key}"
            )
        elif repeat_frame_coverage < 0.995:
            warnings.append(
                f"timing frames below 99.5% of same-cell maximum ({repeat_frame_coverage:.3f}): {run_key}"
            )

    sequence_rows: list[dict[str, object]] = []
    for sequence in sequences:
        for method in methods:
            group = [row for row in run_rows if row["sequence"] == sequence and row["method"] == method]
            if not group:
                continue
            result: dict[str, object] = {
                "sequence": sequence,
                "method": method,
                "successful_repeats": len(group),
                "included_frames": sum(int(row["included_record_count"]) for row in group),
            }
            for field in (
                "total_ms_mean", "total_ms_median", "total_ms_p95", "total_ms_p99", "total_ms_max",
                "preprocess_ms_mean", "core_update_ms_mean", "deadline_miss_rate",
                "compute_utilization", "effective_throughput_hz", "timing_coverage",
            ):
                values = [float(row[field]) for row in group]
                result[f"{field}_repeat_mean"] = statistics.fmean(values)
                result[f"{field}_repeat_std"] = statistics.stdev(values) if len(values) > 1 else 0.0
                result[f"{field}_repeat_worst"] = max(values)
            sequence_rows.append(result)

    overall_rows: list[dict[str, object]] = []
    for method in methods:
        pooled = [
            row for (sequence, item_method, repeat), rows in included_by_run.items()
            if item_method == method for row in rows
        ]
        if not pooled:
            continue
        totals = [float(row["total_ms"]) for row in pooled]
        run_group = [row for row in run_rows if row["method"] == method]
        result = {
            "method": method,
            "sequence_count": len({row["sequence"] for row in run_group}),
            "run_count": len(run_group),
            "included_frames": len(pooled),
            **describe(totals, "total_ms"),
            **describe((float(row["preprocess_ms"]) for row in pooled), "preprocess_ms"),
            **describe((float(row["core_update_ms"]) for row in pooled), "core_update_ms"),
        }
        mean_total = float(result["total_ms_mean"])
        result.update({
            "deadline_ms": args.deadline_ms,
            "deadline_miss_count": sum(value > args.deadline_ms for value in totals),
            "deadline_miss_rate": sum(value > args.deadline_ms for value in totals) / len(totals),
            "compute_utilization": mean_total / args.deadline_ms,
            "effective_throughput_hz": 1000.0 / mean_total,
            "run_mean_ms_mean": statistics.fmean(float(row["total_ms_mean"]) for row in run_group),
            "run_mean_ms_std": statistics.stdev(float(row["total_ms_mean"]) for row in run_group) if len(run_group) > 1 else 0.0,
        })
        overall_rows.append(result)

    lightning_stage_rows: list[dict[str, object]] = []
    for sequence_label in (*sequences, "ALL"):
        lightning_frames = [
            row for (sequence, method, repeat), rows in included_by_run.items()
            if method == "lightning_lm" and (sequence_label == "ALL" or sequence == sequence_label)
            for row in rows
        ]
        if not lightning_frames:
            continue
        total_mean = statistics.fmean(float(row["total_ms"]) for row in lightning_frames)
        for stage in LIGHTNING_STAGES:
            values = [float(row[stage]) for row in lightning_frames]
            stage_stats = describe(values, "stage_ms")
            lightning_stage_rows.append({
                "sequence": sequence_label,
                "stage": stage,
                **stage_stats,
                "total_ms_mean": total_mean,
                "share_of_total_mean": float(stage_stats["stage_ms_mean"]) / total_mean,
            })

    write_csv(args.output_dir / "frame_timing.csv", frame_rows)
    write_csv(args.output_dir / "run_timing_summary.csv", run_rows)
    write_csv(args.output_dir / "sequence_method_summary.csv", sequence_rows)
    write_csv(args.output_dir / "overall_method_summary.csv", overall_rows)
    write_csv(args.output_dir / "lightning_stage_summary.csv", lightning_stage_rows)

    repo = Path(__file__).resolve().parents[3]
    code = repo.parents[1]
    instrumentation = {
        "lightning_lm": repo / "src" / "core" / "lio" / "laser_mapping.cc",
        "fastlio": code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_fastlio" / "src" / "FAST_LIO" / "src" / "laserMapping.cpp",
        "fastlivo2_lio": code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_fastlivo2" / "src" / "FAST-LIVO2" / "src" / "LIVMapper.cpp",
        "voxel_slam_frontend": code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_voxel_slam" / "src" / "Voxel-SLAM" / "VoxelSLAM" / "src" / "voxelslam.cpp",
    }
    binaries = {
        "lightning_lm": repo / "bin" / "run_frontend_offline",
        "fastlio": code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_fastlio" / "devel" / "lib" / "fast_lio" / "fastlio_mapping",
        "fastlivo2_lio": code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_fastlivo2" / "devel" / "lib" / "fast_livo" / "fastlivo_mapping",
        "voxel_slam_frontend": code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_voxel_slam" / "devel" / "lib" / "voxel_slam" / "voxelslam",
    }
    missing_provenance = [str(path) for path in (*instrumentation.values(), *binaries.values()) if not path.is_file()]
    if missing_provenance:
        failures.extend(f"missing provenance artifact: {path}" for path in missing_provenance)
    run_manifest = args.runs_root / "_state" / "experiment_manifest.json"
    manifest = {
        "schema_version": 1,
        "evaluator": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        "run_manifest": {
            "path": str(run_manifest.resolve()),
            "sha256": sha256(run_manifest) if run_manifest.is_file() else None,
        },
        "instrumentation": {
            method: {"path": str(path.resolve()), "sha256": sha256(path) if path.is_file() else None}
            for method, path in instrumentation.items()
        },
        "binaries": {
            method: {"path": str(path.resolve()), "sha256": sha256(path) if path.is_file() else None}
            for method, path in binaries.items()
        },
        "parameters": {
            "sequences": sequences,
            "methods": methods,
            "repeats": args.repeats,
            "warmup_frames_per_run": args.warmup_frames,
            "deadline_ms": args.deadline_ms,
            "timing_boundary": (
                "LiDAR message preprocessing compute plus synchronized LIO core update; excludes bag wait, "
                "ROS publication/trajectory serialization, final optimization, map export and disk output"
            ),
        },
    }
    write_json(args.output_dir / "analysis_manifest.json", manifest)

    validation = {
        "status": "failed" if failures else ("passed_with_warnings" if warnings else "passed"),
        "expected_run_count": expected_runs,
        "evaluated_run_count": len(run_rows),
        "frame_row_count": len(frame_rows),
        "included_frame_count": sum(int(row["included_record_count"]) for row in run_rows),
        "failures": failures,
        "warnings": warnings,
        "checks": [
            "compute-run contract: playback/algorithm/subscribers passed; ROS launch exit 0 or runner shutdown 137",
            "required fields and finite non-negative timings",
            "strictly increasing measurement timestamps",
            "total_ms = preprocess_ms + core_update_ms within 0.05 ms",
            "Lightning core_update_ms equals detailed core stages within 0.05 ms",
            "timing coverage at least 80% of recorded trajectory poses",
            "timing frames at least 98% of the same sequence/method repeat maximum",
            "ten tracking frames removed per run as warm-up by default",
        ],
    }
    write_json(args.output_dir / "validation.json", validation)
    write_json(
        args.output_dir / "metrics.json",
        {
            "run_summary": run_rows,
            "sequence_method_summary": sequence_rows,
            "overall_method_summary": overall_rows,
            "lightning_stage_summary": lightning_stage_rows,
        },
    )
    print(json.dumps(validation, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
