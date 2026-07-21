#!/usr/bin/env python3
"""Evaluate the controlled Lightning-LM backend matrix on common RTK support."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from evaluate_m3dgr_frontend_matrix import (
    GT_PATHS, SEQUENCES, align_se3, gt_orientation_available, interpolate,
    load_tum, mean_std, normalize_gt, parse_metadata, read_resources, rpe,
    sha256, stats, support_mask, unavailable_stats, write_csv,
)


METHODS = ("frontend_only", "legacy_controlled", "local_ba_only", "loop_pgo", "full_backend")
BACKEND_FIELDS = (
    "keyframes", "local_ba_attempts", "local_ba_accepted", "btc_descriptors",
    "btc_candidates", "loops_accepted", "loops_applied", "loops_graph_inliers",
    "hba_runs", "hba_accepted", "local_ba_time_ms", "btc_time_ms",
    "pose_graph_time_ms", "hba_time_ms",
)


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "m3dgr_backend_controlled_20260721",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "analysis" / "m3dgr_backend_controlled_20260721",
    )
    parser.add_argument(
        "--inventory", type=Path,
        default=repo.parents[1] / "_m3dgr_work" / "bench" / "inventory" / "bag_inventory.json",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--max-interpolation-gap", type=float, default=0.2)
    return parser.parse_args()


def parse_flat_yaml(path: Path) -> dict[str, int | float | str]:
    values: dict[str, int | float | str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, raw = (item.strip() for item in line.split(":", 1))
        try:
            values[key] = int(raw)
        except ValueError:
            try:
                values[key] = float(raw)
            except ValueError:
                values[key] = raw
    return values


def parse_legacy(path: Path) -> dict[str, int]:
    text = path.read_text(encoding="utf-8", errors="replace")
    loop_counts = [int(value) for value in re.findall(r"optimize finished, loops:\s*(\d+)", text)]
    return {
        "legacy_loop_success_events": len(re.findall(r"loop_closing\.cc:\d+\]\s+success:", text)),
        "legacy_pose_graph_optimization_events": len(loop_counts),
        "legacy_final_reported_loop_edges": loop_counts[-1] if loop_counts else 0,
    }


def load_loop_events(path: Path, sequence: str, method: str, repeat: int) -> list[dict[str, object]]:
    if not path.is_file():
        return []
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            if row.get("candidate") != "1":
                continue
            rows.append({
                "sequence": sequence, "method": method, "repeat": repeat,
                "current_keyframe": int(row["current_keyframe"]),
                "history_keyframe": int(row["history_keyframe"]),
                "candidate_source": row.get("candidate_source", "btc"),
                "accepted": int(row["accepted"]),
                "optimization_warranted": int(row["optimization_warranted"]),
                "score": float(row["score"]),
                "drift_translation_m": float(row["drift_translation"]),
                "drift_rotation_deg": float(row["drift_rotation_deg"]),
                "journey_span_m": float(row["journey_span"]),
                "drift_ratio": float(row["drift_ratio"]),
                "odom_revisit_distance_m": float(row["odom_revisit_distance"]),
                "plane_icp_observability": float(row["observability"]),
                "plane_icp_matches": int(row["plane_icp_matches"]),
                "confirmation_count": int(row["confirmation_count"]),
                "reason": row["reason"],
            })
    return rows


def main() -> int:
    args = parse_args()
    inventory = {row["sequence"]: row for row in json.loads(args.inventory.read_text(encoding="utf-8"))}
    failures: list[str] = []
    warnings: list[str] = []
    run_data: dict[tuple[str, str, int], tuple[Path, np.ndarray, dict[str, str]]] = {}
    evidence_files: list[Path] = []
    for sequence in SEQUENCES:
        for method in METHODS:
            for repeat in range(1, args.repeats + 1):
                run_dir = args.runs_root / sequence / method / f"repeat_{repeat:02d}"
                metadata_path = run_dir / "run_metadata.txt"
                trajectory_path = run_dir / "results" / "trajectory_slam_opt.tum"
                resource_path = run_dir / "resource_summary.json"
                if not all(path.is_file() for path in (metadata_path, trajectory_path, resource_path)):
                    failures.append(f"missing run artifact: {sequence}/{method}/repeat_{repeat:02d}")
                    continue
                metadata = parse_metadata(metadata_path)
                if not (
                    metadata.get("completion") == "reached_final_lidar"
                    and metadata.get("algorithm_rc") == "0"
                    and metadata.get("watchdog_status") == "completed"
                    and metadata.get("backend_evaluation_only") == "true"
                    and metadata.get("invalid_count") == "0"
                    and metadata.get("nonmonotonic_count") == "0"
                    and metadata.get("excessive_output_gap_count") == "0"
                ):
                    failures.append(f"failed run contract: {sequence}/{method}/repeat_{repeat:02d}")
                    continue
                run_data[(sequence, method, repeat)] = (run_dir, load_tum(trajectory_path), metadata)

    run_rows: list[dict[str, object]] = []
    aligned_rows: list[dict[str, object]] = []
    activity_rows: list[dict[str, object]] = []
    loop_rows: list[dict[str, object]] = []
    endpoint_rows: list[dict[str, object]] = []
    gt_preprocessing: dict[str, object] = {}
    for sequence in SEQUENCES:
        keys = [(sequence, method, repeat) for method in METHODS for repeat in range(1, args.repeats + 1)]
        if not all(key in run_data for key in keys):
            failures.append(f"incomplete sequence matrix: {sequence}")
            continue
        gt, diagnostics = normalize_gt(load_tum(GT_PATHS[sequence], require_monotonic=False))
        gt_preprocessing[sequence] = diagnostics
        reference_run_dir = run_data[(sequence, "frontend_only", 1)][0]
        raw_lio = load_tum(reference_run_dir / "results" / "trajectory_slam_keyframes_lio.tum")
        endpoint_rows.append({
            "sequence": sequence,
            "raw_lio_start_end_distance_m": float(np.linalg.norm(raw_lio[-1, 1:4] - raw_lio[0, 1:4])),
            "official_rtk_start_end_distance_m": float(np.linalg.norm(gt[-1, 1:4] - gt[0, 1:4])),
            "raw_lio_pose_count": len(raw_lio),
            "official_rtk_pose_count_after_normalization": len(gt),
            "raw_lio_sha256": sha256(reference_run_dir / "results" / "trajectory_slam_keyframes_lio.tum"),
        })
        common = np.ones(len(gt), dtype=bool)
        for key in keys:
            common &= support_mask(run_data[key][1], gt[:, 0], args.max_interpolation_gap)
        common_gt = gt[common]
        if len(common_gt) < 100:
            failures.append(f"insufficient common support: {sequence} ({len(common_gt)})")
            continue
        truth_position = common_gt[:, 1:4]
        truth_rotation = Rotation.from_quat(common_gt[:, 4:8])
        orientation_available = gt_orientation_available(common_gt)
        if not orientation_available:
            warnings.append(f"{sequence}: RTK quaternion is a placeholder; rotation metrics are N/A")
        for method in METHODS:
            for repeat in range(1, args.repeats + 1):
                run_dir, trajectory, metadata = run_data[(sequence, method, repeat)]
                estimate_position, estimate_rotation = interpolate(trajectory, common_gt[:, 0])
                alignment_rotation, alignment_translation = align_se3(estimate_position, truth_position)
                aligned_position = alignment_rotation.apply(estimate_position) + alignment_translation
                aligned_rotation = alignment_rotation * estimate_rotation
                translation_error = np.linalg.norm(aligned_position - truth_position, axis=1)
                rotation_error = (
                    np.degrees((truth_rotation.inv() * aligned_rotation).magnitude())
                    if orientation_available else None
                )
                ate = stats(translation_error, "m")
                ate_rotation = stats(rotation_error, "deg") if rotation_error is not None else unavailable_stats("deg")
                rpe_1m = rpe(aligned_position, aligned_rotation, truth_position, truth_rotation, 1.0, orientation_available)
                rpe_10m = rpe(aligned_position, aligned_rotation, truth_position, truth_rotation, 10.0, orientation_available)
                resources = read_resources(run_dir / "resource_summary.json")
                output_frames = int(metadata["trajectory_lines"])
                wall_time = float(metadata["wall_time_s"])
                optimized_path = run_dir / "results" / "trajectory_slam_keyframes_opt.tum"
                lio_path = run_dir / "results" / "trajectory_slam_keyframes_lio.tum"
                row = {
                    "sequence": sequence, "method": method, "repeat": repeat,
                    "common_gt_pose_count": len(common_gt),
                    "common_gt_coverage": len(common_gt) / len(gt),
                    "ate_rmse_m": ate["rmse_m"], "ate_median_m": ate["median_m"],
                    "ate_p95_m": ate["p95_m"], "ate_max_m": ate["max_m"],
                    "ate_rotation_rmse_deg": ate_rotation["rmse_deg"],
                    "rpe_1m_translation_rmse_m": rpe_1m["translation"]["rmse_m"],
                    "rpe_10m_translation_rmse_m": rpe_10m["translation"]["rmse_m"],
                    "wall_time_s": wall_time,
                    "wall_ms_per_output_frame": wall_time * 1000.0 / output_frames,
                    "output_frames": output_frames,
                    "expected_lidar_frames": int(inventory[sequence]["lidar_count"]),
                    "frontend_lio_sha256": sha256(lio_path),
                    "optimized_keyframes_sha256": sha256(optimized_path),
                    **resources,
                }
                run_rows.append(row)
                if repeat == 1:
                    aligned_rows.extend({
                        "sequence": sequence, "method": method, "timestamp_s": float(timestamp),
                        "truth_x_m": float(truth[0]), "truth_y_m": float(truth[1]), "truth_z_m": float(truth[2]),
                        "estimate_x_m": float(estimate[0]), "estimate_y_m": float(estimate[1]),
                        "estimate_z_m": float(estimate[2]), "translation_error_m": float(error),
                    } for timestamp, truth, estimate, error in zip(
                        common_gt[:, 0], truth_position, aligned_position, translation_error, strict=True
                    ))

                if method == "legacy_controlled":
                    log = run_dir / "logs" / "algorithm.stderr.log"
                    activity_rows.append({"sequence": sequence, "method": method, "repeat": repeat, **parse_legacy(log)})
                    evidence_files.append(log)
                elif method in {"local_ba_only", "loop_pgo", "full_backend"}:
                    summary_path = run_dir / "data" / "new_map" / "backend_diagnostics" / "backend_summary.yaml"
                    values = parse_flat_yaml(summary_path)
                    missing = [field for field in BACKEND_FIELDS if field not in values]
                    if missing:
                        failures.append(f"missing backend fields: {sequence}/{method}/{repeat}: {missing}")
                    activity_rows.append({
                        "sequence": sequence, "method": method, "repeat": repeat,
                        **{field: values.get(field) for field in BACKEND_FIELDS},
                    })
                    evidence_files.append(summary_path)
                    candidate_path = summary_path.parent / "btc_loop_candidates.csv"
                    loop_rows.extend(load_loop_events(candidate_path, sequence, method, repeat))
                    if candidate_path.is_file():
                        evidence_files.append(candidate_path)

    numeric = (
        "ate_rmse_m", "ate_median_m", "ate_p95_m", "ate_max_m",
        "rpe_1m_translation_rmse_m", "rpe_10m_translation_rmse_m",
        "wall_time_s", "wall_ms_per_output_frame", "mean_cpu_cores", "peak_cpu_cores",
        "p95_cpu_cores", "mean_rss_mb", "peak_rss_mb",
    )
    summary_rows: list[dict[str, object]] = []
    for sequence in SEQUENCES:
        for method in METHODS:
            group = [row for row in run_rows if row["sequence"] == sequence and row["method"] == method]
            out: dict[str, object] = {
                "sequence": sequence, "method": method, "successful_repeats": len(group),
                "unique_frontend_lio_hashes": len({str(row["frontend_lio_sha256"]) for row in group}),
                "unique_optimized_hashes": len({str(row["optimized_keyframes_sha256"]) for row in group}),
            }
            for field in numeric:
                mean, std = mean_std([float(row[field]) for row in group if row.get(field) is not None])
                out[f"{field}_mean"] = mean
                out[f"{field}_std"] = std
            summary_rows.append(out)

    effect_rows: list[dict[str, object]] = []
    comparisons = (
        ("legacy_controlled", "frontend_only", "legacy_minus_frontend"),
        ("local_ba_only", "frontend_only", "local_ba_minus_frontend"),
        ("loop_pgo", "local_ba_only", "loop_pgo_minus_local_ba"),
        ("full_backend", "loop_pgo", "hba_minus_loop_pgo"),
        ("full_backend", "frontend_only", "full_minus_frontend"),
    )
    for sequence in SEQUENCES:
        by_method = {row["method"]: row for row in summary_rows if row["sequence"] == sequence}
        for treatment, control, name in comparisons:
            effect_rows.append({
                "sequence": sequence, "comparison": name,
                "ate_rmse_delta_m": float(by_method[treatment]["ate_rmse_m_mean"]) - float(by_method[control]["ate_rmse_m_mean"]),
                "ate_rmse_delta_percent": 100.0 * (
                    float(by_method[treatment]["ate_rmse_m_mean"]) / float(by_method[control]["ate_rmse_m_mean"]) - 1.0
                ),
                "wall_ms_per_frame_delta": float(by_method[treatment]["wall_ms_per_output_frame_mean"]) - float(by_method[control]["wall_ms_per_output_frame_mean"]),
                "mean_cpu_cores_delta": float(by_method[treatment]["mean_cpu_cores_mean"]) - float(by_method[control]["mean_cpu_cores_mean"]),
                "peak_rss_mb_delta": float(by_method[treatment]["peak_rss_mb_mean"]) - float(by_method[control]["peak_rss_mb_mean"]),
            })

    frontend_hash_groups = {
        sequence: sorted({str(row["frontend_lio_sha256"]) for row in run_rows if row["sequence"] == sequence})
        for sequence in SEQUENCES
    }
    if any(len(values) != 1 for values in frontend_hash_groups.values()):
        failures.append(f"frontend LIO hash mismatch: {frontend_hash_groups}")
    for sequence in SEQUENCES:
        for method in METHODS:
            if len({str(row["optimized_keyframes_sha256"]) for row in run_rows if row["sequence"] == sequence and row["method"] == method}) != 1:
                warnings.append(
                    f"{sequence}/{method}: optimized trajectory differs across repeats; "
                    "report precision as mean +/- sample standard deviation"
                )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "run_metrics.csv", run_rows)
    write_csv(args.output_dir / "summary_metrics.csv", summary_rows)
    write_csv(args.output_dir / "aligned_samples.csv", aligned_rows)
    write_csv(args.output_dir / "backend_activity.csv", activity_rows)
    write_csv(args.output_dir / "loop_candidates.csv", loop_rows)
    write_csv(args.output_dir / "endpoint_loop_diagnostics.csv", endpoint_rows)
    write_csv(args.output_dir / "paired_effects.csv", effect_rows)
    validation = {
        "status": "failed" if failures else ("passed_with_warnings" if warnings else "passed"),
        "required_run_count": len(SEQUENCES) * len(METHODS) * args.repeats,
        "evaluated_run_count": len(run_rows), "failures": failures, "warnings": warnings,
        "frontend_lio_hashes": frontend_hash_groups,
        "alignment": "fixed-scale SE(3) on common official-RTK support",
        "timing_contract": "unpaced backend-evaluation-only wall ms/frame; no RTF",
    }
    (args.output_dir / "validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    manifest = {
        "schema_version": 1,
        "evaluator": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
        "run_manifest": {"path": str((args.runs_root / "_state" / "experiment_manifest.json").resolve()),
                         "sha256": sha256(args.runs_root / "_state" / "experiment_manifest.json")},
        "inventory": {"path": str(args.inventory.resolve()), "sha256": sha256(args.inventory)},
        "ground_truth": {sequence: {"path": str(GT_PATHS[sequence].resolve()), "sha256": sha256(GT_PATHS[sequence]),
                                     "preprocessing": gt_preprocessing.get(sequence)} for sequence in SEQUENCES},
        "backend_evidence": {str(path.resolve()): sha256(path) for path in sorted(set(evidence_files))},
    }
    (args.output_dir / "analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
