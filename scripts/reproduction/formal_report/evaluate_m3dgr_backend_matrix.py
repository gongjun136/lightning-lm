#!/usr/bin/env python3
"""Evaluate the formal M3DGR backend matrix on common official-GT support."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from evaluate_m3dgr_frontend_matrix import (
    GT_PATHS,
    SEQUENCES,
    align_se3,
    gt_orientation_available,
    interpolate,
    load_tum,
    mean_std,
    normalize_gt,
    parse_metadata,
    read_resources,
    rpe,
    sha256,
    stats,
    support_mask,
    unavailable_stats,
    write_csv,
)


METHODS = (
    "lightning_frontend",
    "lightning_legacy",
    "lightning_new_backend",
    "voxel_slam_full",
)

BACKEND_ACTIVITY_FIELDS = (
    "keyframes",
    "local_ba_attempts",
    "local_ba_accepted",
    "btc_descriptors",
    "btc_candidates",
    "loops_accepted",
    "loops_applied",
    "loops_graph_inliers",
    "hba_runs",
    "hba_accepted",
    "local_ba_time_ms",
    "btc_time_ms",
    "pose_graph_time_ms",
    "hba_time_ms",
)


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "m3dgr_backend",
    )
    parser.add_argument(
        "--frontend-runs-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "m3dgr_frontend_v2",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "analysis" / "m3dgr_backend",
    )
    parser.add_argument(
        "--inventory", type=Path,
        default=repo.parents[1] / "_m3dgr_work" / "bench" / "inventory" / "bag_inventory.json",
    )
    parser.add_argument("--sequences", default=",".join(SEQUENCES))
    parser.add_argument("--methods", default=",".join(METHODS))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--max-interpolation-gap", type=float, default=0.2)
    return parser.parse_args()


def trajectory_path(run_dir: Path, method: str) -> Path:
    if method == "lightning_frontend":
        return run_dir / "trajectory_mid360.tum"
    if method == "voxel_slam_full":
        return run_dir / "results" / "trajectory_slam_opt.tum"
    return run_dir / "results" / "trajectory_slam_opt.tum"


def method_run_dir(args: argparse.Namespace, sequence: str, method: str, repeat: int) -> Path:
    if method == "lightning_frontend":
        return args.frontend_runs_root / sequence / "lightning_lm" / f"repeat_{repeat:02d}"
    return args.runs_root / sequence / method / f"repeat_{repeat:02d}"


def run_contract(metadata: dict[str, str], method: str) -> tuple[bool, list[str]]:
    warnings: list[str] = []
    valid_trajectory = (
        metadata.get("invalid_count") == "0"
        and metadata.get("nonmonotonic_count") == "0"
        and metadata.get("excessive_output_gap_count") == "0"
    )
    if method == "lightning_frontend":
        return (
            valid_trajectory
            and metadata.get("completion") == "reached_final_lidar"
            and metadata.get("method") == "lightning_lm"
            and metadata.get("algorithm_rc") == "0"
            and metadata.get("watchdog_status") == "completed"
        ), warnings
    if method == "voxel_slam_full":
        try:
            output_frames = int(metadata.get("trajectory_lines", "0"))
            expected_frames = int(metadata.get("expected_lidar_frames", "0"))
            tail_gap_s = float(metadata["expected_last_lidar_end_s"]) - float(metadata["last_stamp"])
        except (KeyError, TypeError, ValueError):
            return False, warnings
        coverage = output_frames / expected_frames if expected_frames > 0 else 0.0
        passed = (
            valid_trajectory
            and metadata.get("method") == "voxel_slam_full_backend"
            and metadata.get("pcd_count") == metadata.get("trajectory_lines")
            and coverage >= 0.98
            and 0.0 <= tail_gap_s <= 2.0
        )
        if passed and metadata.get("completion") != "reached_final_lidar":
            warnings.append(
                f"Voxel-SLAM initialization tail: coverage={coverage:.6f}, tail_gap_s={tail_gap_s:.6f}"
            )
        return passed, warnings
    return (
        valid_trajectory
        and metadata.get("completion") == "reached_final_lidar"
        and metadata.get("method") == "lightning_lm_offline_slam_map_export"
        and metadata.get("algorithm_rc") == "0"
        and metadata.get("watchdog_status") == "completed"
        and int(metadata.get("map_chunk_count", "0")) > 0
    ), warnings


def parse_backend_summary(path: Path) -> dict[str, int | float | str]:
    """Parse the flat scalar YAML emitted by BackendPipeline without a YAML dependency."""
    parsed: dict[str, int | float | str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#") or ":" not in raw_line:
            continue
        key, raw_value = (item.strip() for item in raw_line.split(":", 1))
        try:
            value: int | float | str = int(raw_value)
        except ValueError:
            try:
                value = float(raw_value)
            except ValueError:
                value = raw_value
        parsed[key] = value
    return parsed


def parse_legacy_activity(path: Path) -> dict[str, int]:
    """Extract auditable legacy loop/optimization counters from its own log vocabulary."""
    text = path.read_text(encoding="utf-8", errors="replace")
    loop_counts = [
        int(value) for value in re.findall(r"optimize finished, loops:\s*(\d+)", text)
    ]
    return {
        "legacy_loop_success_events": len(re.findall(r"loop_closing\.cc:\d+\]\s+success:", text)),
        "legacy_pose_graph_optimization_events": len(loop_counts),
        "legacy_final_reported_loop_edges": loop_counts[-1] if loop_counts else 0,
    }


def main() -> int:
    args = parse_args()
    sequences = tuple(item for item in args.sequences.split(",") if item)
    methods = tuple(item for item in args.methods.split(",") if item)
    if set(sequences) - set(SEQUENCES) or set(methods) - set(METHODS):
        raise SystemExit("unknown sequence or method")
    inventory_rows = json.loads(args.inventory.read_text(encoding="utf-8"))
    inventory = {row["sequence"]: row for row in inventory_rows}

    runs: dict[str, dict[tuple[str, int], tuple[Path, np.ndarray, dict[str, str]]]] = {}
    failures: list[str] = []
    warnings: list[str] = []
    for sequence in sequences:
        sequence_runs = {}
        for method in methods:
            for repeat in range(1, args.repeats + 1):
                run_dir = method_run_dir(args, sequence, method, repeat)
                metadata_path = run_dir / "run_metadata.txt"
                trajectory = trajectory_path(run_dir, method)
                resources = run_dir / "resource_summary.json"
                if not all(path.is_file() for path in (metadata_path, trajectory, resources)):
                    failures.append(f"missing artifacts: {sequence}/{method}/repeat_{repeat:02d}")
                    continue
                metadata = parse_metadata(metadata_path)
                passed, run_warnings = run_contract(metadata, method)
                if not passed:
                    failures.append(f"failed run contract: {sequence}/{method}/repeat_{repeat:02d}")
                    continue
                warnings.extend(
                    f"{sequence}/{method}/repeat_{repeat:02d}: {warning}" for warning in run_warnings
                )
                sequence_runs[(method, repeat)] = (run_dir, load_tum(trajectory), metadata)
        runs[sequence] = sequence_runs

    run_rows: list[dict[str, object]] = []
    aligned_sample_rows: list[dict[str, object]] = []
    backend_activity_rows: list[dict[str, object]] = []
    diagnostic_files: list[Path] = []
    result_sequences: dict[str, object] = {}
    gt_preprocessing: dict[str, dict[str, int | str]] = {}
    for sequence in sequences:
        sequence_runs = runs[sequence]
        required = [(method, repeat) for method in methods for repeat in range(1, args.repeats + 1)]
        if not all(key in sequence_runs for key in required):
            failures.append(f"incomplete common-support matrix: {sequence}")
            continue
        gt, gt_diagnostics = normalize_gt(load_tum(GT_PATHS[sequence], require_monotonic=False))
        gt_preprocessing[sequence] = gt_diagnostics
        if (
            int(gt_diagnostics["raw_negative_transition_count"]) > 0
            or int(gt_diagnostics["duplicate_timestamp_rows_removed"]) > 0
        ):
            warnings.append(
                f"{sequence}: official GT normalized by stable timestamp sort and exact-time deduplication "
                f"({gt_diagnostics['raw_negative_transition_count']} negative transitions, "
                f"{gt_diagnostics['duplicate_timestamp_rows_removed']} duplicate rows removed)"
            )
        common = np.ones(len(gt), dtype=bool)
        for _, trajectory, _ in sequence_runs.values():
            common &= support_mask(trajectory, gt[:, 0], args.max_interpolation_gap)
        common_gt = gt[common]
        if len(common_gt) < 100:
            failures.append(f"insufficient common GT support: {sequence} ({len(common_gt)})")
            continue
        truth_position = common_gt[:, 1:4]
        truth_rotation = Rotation.from_quat(common_gt[:, 4:8])
        orientation_available = gt_orientation_available(common_gt)
        gt_diagnostics["orientation_available"] = orientation_available
        if not orientation_available:
            warnings.append(
                f"{sequence}: official RTK quaternions are all identity placeholders; "
                "rotation metrics are N/A and translation RPE uses aligned global displacements"
            )
        sequence_payload: dict[str, object] = {
            "common_gt_pose_count": len(common_gt),
            "common_gt_coverage": len(common_gt) / len(gt),
            "gt_orientation_available": orientation_available,
            "methods": {},
        }
        for (method, repeat), (run_dir, trajectory, metadata) in sequence_runs.items():
            estimate_position, estimate_rotation = interpolate(trajectory, common_gt[:, 0])
            alignment_rotation, alignment_translation = align_se3(estimate_position, truth_position)
            aligned_position = alignment_rotation.apply(estimate_position) + alignment_translation
            aligned_rotation = alignment_rotation * estimate_rotation
            translation_error = np.linalg.norm(aligned_position - truth_position, axis=1)
            rotation_error = (
                np.degrees((truth_rotation.inv() * aligned_rotation).magnitude())
                if orientation_available else None
            )
            ate_translation = stats(translation_error, "m")
            ate_rotation = stats(rotation_error, "deg") if rotation_error is not None else unavailable_stats("deg")
            rpe_1m = rpe(
                aligned_position, aligned_rotation, truth_position, truth_rotation, 1.0, orientation_available
            )
            rpe_10m = rpe(
                aligned_position, aligned_rotation, truth_position, truth_rotation, 10.0, orientation_available
            )
            resources = read_resources(run_dir / "resource_summary.json")
            wall_time = float(metadata["wall_time_s"])
            sensor_duration = float(metadata["sensor_duration_s"])
            output_frames = int(metadata["trajectory_lines"])
            expected_frames = int(inventory[sequence]["lidar_count"])
            row = {
                "sequence": sequence,
                "method": method,
                "repeat": repeat,
                "common_gt_pose_count": len(common_gt),
                "common_gt_coverage": len(common_gt) / len(gt),
                "gt_orientation_available": orientation_available,
                "ate_rmse_m": ate_translation["rmse_m"],
                "ate_median_m": ate_translation["median_m"],
                "ate_p95_m": ate_translation["p95_m"],
                "ate_max_m": ate_translation["max_m"],
                "ate_rotation_rmse_deg": ate_rotation["rmse_deg"],
                "rpe_1m_translation_rmse_m": rpe_1m["translation"]["rmse_m"],
                "rpe_1m_rotation_rmse_deg": rpe_1m["rotation"]["rmse_deg"],
                "rpe_10m_translation_rmse_m": rpe_10m["translation"]["rmse_m"],
                "rpe_10m_rotation_rmse_deg": rpe_10m["rotation"]["rmse_deg"],
                "wall_time_s": wall_time,
                "realtime_factor": wall_time / sensor_duration,
                "output_frames": output_frames,
                "expected_lidar_frames": expected_frames,
                "output_deficit_ratio": max(0.0, 1.0 - output_frames / expected_frames),
                "wall_ms_per_output_frame": wall_time * 1000.0 / output_frames,
                **resources,
            }
            run_rows.append(row)
            payload = {
                **row,
                "ate_translation": ate_translation,
                "ate_rotation": ate_rotation,
                "rpe_1m": rpe_1m,
                "rpe_10m": rpe_10m,
            }
            sequence_payload["methods"].setdefault(method, []).append(payload)
            aligned_sample_rows.extend({
                "sequence": sequence,
                "method": method,
                "repeat": repeat,
                "timestamp_s": float(timestamp),
                "truth_x_m": float(truth[0]),
                "truth_y_m": float(truth[1]),
                "truth_z_m": float(truth[2]),
                "estimate_x_m": float(estimate[0]),
                "estimate_y_m": float(estimate[1]),
                "estimate_z_m": float(estimate[2]),
                "translation_error_m": float(error),
                "rotation_error_deg": (
                    float(rotation_error_value) if rotation_error_value is not None else None
                ),
            } for timestamp, truth, estimate, error, rotation_error_value in zip(
                common_gt[:, 0], truth_position, aligned_position, translation_error,
                rotation_error if rotation_error is not None else [None] * len(common_gt), strict=True
            ))
        result_sequences[sequence] = sequence_payload

        for repeat in range(1, args.repeats + 1):
            legacy_dir = sequence_runs[("lightning_legacy", repeat)][0]
            legacy_log = legacy_dir / "logs" / "algorithm.stderr.log"
            if not legacy_log.is_file():
                failures.append(f"missing legacy diagnostics: {sequence}/repeat_{repeat:02d}")
            else:
                diagnostic_files.append(legacy_log)
                backend_activity_rows.append({
                    "sequence": sequence,
                    "method": "lightning_legacy",
                    "repeat": repeat,
                    **parse_legacy_activity(legacy_log),
                })

            new_dir = sequence_runs[("lightning_new_backend", repeat)][0]
            new_summary = new_dir / "data" / "new_map" / "backend_diagnostics" / "backend_summary.yaml"
            if not new_summary.is_file():
                failures.append(f"missing new-backend diagnostics: {sequence}/repeat_{repeat:02d}")
            else:
                diagnostic_files.append(new_summary)
                parsed = parse_backend_summary(new_summary)
                missing_fields = [field for field in BACKEND_ACTIVITY_FIELDS if field not in parsed]
                if missing_fields:
                    failures.append(
                        f"incomplete new-backend diagnostics: {sequence}/repeat_{repeat:02d}: "
                        + ",".join(missing_fields)
                    )
                backend_activity_rows.append({
                    "sequence": sequence,
                    "method": "lightning_new_backend",
                    "repeat": repeat,
                    **{field: parsed.get(field) for field in BACKEND_ACTIVITY_FIELDS},
                })

    numeric_fields = [
        "ate_rmse_m", "ate_median_m", "ate_p95_m", "ate_max_m", "ate_rotation_rmse_deg",
        "rpe_1m_translation_rmse_m", "rpe_1m_rotation_rmse_deg",
        "rpe_10m_translation_rmse_m", "rpe_10m_rotation_rmse_deg",
        "wall_time_s", "realtime_factor", "output_deficit_ratio", "wall_ms_per_output_frame",
        "mean_cpu_cores", "peak_cpu_cores", "p95_cpu_cores", "mean_rss_mb", "peak_rss_mb",
    ]
    summary_rows: list[dict[str, object]] = []
    for sequence in sequences:
        for method in methods:
            group = [row for row in run_rows if row["sequence"] == sequence and row["method"] == method]
            summary: dict[str, object] = {
                "sequence": sequence, "method": method, "successful_repeats": len(group)
            }
            for field in numeric_fields:
                values = [float(item[field]) for item in group if item.get(field) is not None]
                mean, std = mean_std(values)
                summary[f"{field}_mean"] = mean
                summary[f"{field}_std"] = std
                summary[f"{field}_worst"] = max(values, default=None)
            summary_rows.append(summary)

    activity_numeric_fields = tuple(
        sorted({
            key
            for row in backend_activity_rows
            for key, value in row.items()
            if key not in {"sequence", "method", "repeat"} and isinstance(value, (int, float))
        })
    )
    backend_activity_summary_rows: list[dict[str, object]] = []
    for sequence in sequences:
        for method in ("lightning_legacy", "lightning_new_backend"):
            group = [
                row for row in backend_activity_rows
                if row["sequence"] == sequence and row["method"] == method
            ]
            summary: dict[str, object] = {
                "sequence": sequence, "method": method, "successful_repeats": len(group)
            }
            for field in activity_numeric_fields:
                values = [float(item[field]) for item in group if item.get(field) is not None]
                mean, std = mean_std(values)
                summary[f"{field}_mean"] = mean
                summary[f"{field}_std"] = std
                summary[f"{field}_worst"] = max(values, default=None)
            backend_activity_summary_rows.append(summary)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "run_metrics.csv", run_rows)
    write_csv(args.output_dir / "summary_metrics.csv", summary_rows)
    write_csv(args.output_dir / "aligned_samples.csv", aligned_sample_rows)
    write_csv(args.output_dir / "backend_activity.csv", backend_activity_rows)
    write_csv(args.output_dir / "backend_activity_summary.csv", backend_activity_summary_rows)
    (args.output_dir / "metrics.json").write_text(
        json.dumps(
            {"sequences": result_sequences, "failures": failures, "warnings": warnings},
            ensure_ascii=False,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    run_manifest = args.runs_root / "_state" / "experiment_manifest.json"
    frontend_run_manifest = args.frontend_runs_root / "_state" / "experiment_manifest.json"
    analysis_manifest = {
        "schema_version": 1,
        "evaluator": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
        "run_manifest": {
            "path": str(run_manifest.resolve()),
            "sha256": sha256(run_manifest) if run_manifest.is_file() else None,
        },
        "frontend_run_manifest": {
            "path": str(frontend_run_manifest.resolve()),
            "sha256": sha256(frontend_run_manifest) if frontend_run_manifest.is_file() else None,
        },
        "backend_diagnostics": {
            str(path.resolve()): sha256(path) for path in sorted(set(diagnostic_files))
        },
        "inventory": {"path": str(args.inventory.resolve()), "sha256": sha256(args.inventory)},
        "ground_truth": {
            sequence: {
                "path": str(GT_PATHS[sequence].resolve()),
                "sha256": sha256(GT_PATHS[sequence]),
                "preprocessing": gt_preprocessing.get(sequence),
            }
            for sequence in sequences
        },
        "parameters": {
            "sequences": sequences, "methods": methods, "repeats": args.repeats,
            "max_interpolation_gap_s": args.max_interpolation_gap,
            "alignment": "fixed-scale SE(3)",
            "ground_truth_orientation_contract": (
                "all-identity RTK quaternions are placeholders: rotation metrics N/A; "
                "translation RPE uses aligned global displacement vectors"
            ),
        },
    }
    (args.output_dir / "analysis_manifest.json").write_text(
        json.dumps(analysis_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    validation = {
        "status": "failed" if failures else ("passed_with_warnings" if warnings else "passed"),
        "required_run_count": len(sequences) * len(methods) * args.repeats,
        "evaluated_run_count": len(run_rows),
        "failures": failures,
        "warnings": warnings,
        "alignment": "fixed-scale SE(3)",
        "ground_truth_orientation_contract": (
            "all-identity RTK quaternions are placeholders: rotation metrics N/A; "
            "translation RPE uses aligned global displacement vectors"
        ),
        "common_support": (
            "intersection of interpolable official-GT timestamps across the Lightning frontend, "
            "legacy backend, new backend and Voxel-SLAM full runs per sequence"
        ),
        "backend_activity": {
            "legacy": "loop success/pose-graph optimization/final edge counters parsed from algorithm.stderr.log",
            "new": "BA/BTC/pose-graph/HBA counters parsed from backend_summary.yaml",
        },
    }
    (args.output_dir / "validation.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(validation, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
