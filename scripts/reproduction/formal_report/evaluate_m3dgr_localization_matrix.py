#!/usr/bin/env python3
"""Evaluate M3DGR localization for Lightning- and Voxel-built maps."""

from __future__ import annotations

import argparse
import json
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


MAP_VARIANTS = ("lightning_new_backend", "voxel_slam_full")


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "m3dgr_localization",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "analysis" / "m3dgr_localization",
    )
    parser.add_argument("--sequences", default=",".join(SEQUENCES))
    parser.add_argument("--map-variants", default=",".join(MAP_VARIANTS))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--max-interpolation-gap", type=float, default=0.2)
    return parser.parse_args()


def run_passed(metadata: dict[str, str]) -> bool:
    return (
        metadata.get("method") == "lightning_lm_offline_localization"
        and metadata.get("completion") == "reached_final_lidar"
        and metadata.get("algorithm_rc") == "0"
        and metadata.get("watchdog_status") == "completed"
        and metadata.get("invalid_count") == "0"
        and metadata.get("nonmonotonic_count") == "0"
        and metadata.get("excessive_output_gap_count") == "0"
    )


def aligned_metrics(
    trajectory: np.ndarray,
    common_gt: np.ndarray,
    truth_position: np.ndarray,
    truth_rotation: Rotation,
    orientation_available: bool,
) -> tuple[dict[str, object], np.ndarray, np.ndarray, np.ndarray | None]:
    position, rotation = interpolate(trajectory, common_gt[:, 0])
    alignment_rotation, alignment_translation = align_se3(position, truth_position)
    position = alignment_rotation.apply(position) + alignment_translation
    rotation = alignment_rotation * rotation
    translation_error = np.linalg.norm(position - truth_position, axis=1)
    rotation_error = (
        np.degrees((truth_rotation.inv() * rotation).magnitude())
        if orientation_available else None
    )
    return {
        "ate_translation": stats(translation_error, "m"),
        "ate_rotation": stats(rotation_error, "deg") if rotation_error is not None else unavailable_stats("deg"),
        "rpe_1m": rpe(position, rotation, truth_position, truth_rotation, 1.0, orientation_available),
        "rpe_10m": rpe(position, rotation, truth_position, truth_rotation, 10.0, orientation_available),
    }, position, translation_error, rotation_error


def main() -> int:
    args = parse_args()
    sequences = tuple(item for item in args.sequences.split(",") if item)
    variants = tuple(item for item in args.map_variants.split(",") if item)
    if set(sequences) - set(SEQUENCES) or set(variants) - set(MAP_VARIANTS):
        raise SystemExit("unknown sequence or map variant")

    runs: dict[str, dict[tuple[str, int], tuple[Path, np.ndarray, np.ndarray, dict[str, str], dict]]] = {}
    failures: list[str] = []
    warnings: list[str] = []
    for sequence in sequences:
        sequence_runs = {}
        for variant in variants:
            for repeat in range(1, args.repeats + 1):
                run_dir = args.runs_root / sequence / variant / f"repeat_{repeat:02d}"
                required = [
                    run_dir / "run_metadata.txt",
                    run_dir / "resource_summary.json",
                    run_dir / "results" / "trajectory_loc.tum",
                    run_dir / "results" / "trajectory_lidar_loc.tum",
                    run_dir / "results" / "localization_summary.json",
                ]
                if not all(path.is_file() for path in required):
                    failures.append(f"missing artifacts: {sequence}/{variant}/repeat_{repeat:02d}")
                    continue
                metadata = parse_metadata(required[0])
                if not run_passed(metadata):
                    failures.append(f"failed run contract: {sequence}/{variant}/repeat_{repeat:02d}")
                    continue
                summary = json.loads(required[4].read_text(encoding="utf-8"))
                sequence_runs[(variant, repeat)] = (
                    run_dir, load_tum(required[2]), load_tum(required[3]), metadata, summary
                )
        runs[sequence] = sequence_runs

    run_rows: list[dict[str, object]] = []
    aligned_sample_rows: list[dict[str, object]] = []
    result_sequences: dict[str, object] = {}
    gt_preprocessing: dict[str, dict[str, int | str]] = {}
    for sequence in sequences:
        sequence_runs = runs[sequence]
        required_keys = [(variant, repeat) for variant in variants for repeat in range(1, args.repeats + 1)]
        if not all(key in sequence_runs for key in required_keys):
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
        for _, final_trajectory, raw_trajectory, _, _ in sequence_runs.values():
            common &= support_mask(final_trajectory, gt[:, 0], args.max_interpolation_gap)
            common &= support_mask(raw_trajectory, gt[:, 0], args.max_interpolation_gap)
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
            "map_variants": {},
        }
        for (variant, repeat), (run_dir, final_trajectory, raw_trajectory, metadata, summary) in sequence_runs.items():
            final_metrics, aligned_position, translation_error, rotation_error = aligned_metrics(
                final_trajectory, common_gt, truth_position, truth_rotation, orientation_available
            )
            raw_metrics, _, _, _ = aligned_metrics(
                raw_trajectory, common_gt, truth_position, truth_rotation, orientation_available
            )
            localization = summary["localization"]
            resources = read_resources(run_dir / "resource_summary.json")
            wall_time = float(metadata["wall_time_s"])
            sensor_duration = float(metadata["sensor_duration_s"])
            row = {
                "sequence": sequence,
                "map_variant": variant,
                "repeat": repeat,
                "common_gt_pose_count": len(common_gt),
                "common_gt_coverage": len(common_gt) / len(gt),
                "gt_orientation_available": orientation_available,
                "ate_rmse_m": final_metrics["ate_translation"]["rmse_m"],
                "ate_median_m": final_metrics["ate_translation"]["median_m"],
                "ate_p95_m": final_metrics["ate_translation"]["p95_m"],
                "ate_max_m": final_metrics["ate_translation"]["max_m"],
                "ate_rotation_rmse_deg": final_metrics["ate_rotation"]["rmse_deg"],
                "rpe_1m_translation_rmse_m": final_metrics["rpe_1m"]["translation"]["rmse_m"],
                "rpe_1m_rotation_rmse_deg": final_metrics["rpe_1m"]["rotation"]["rmse_deg"],
                "rpe_10m_translation_rmse_m": final_metrics["rpe_10m"]["translation"]["rmse_m"],
                "rpe_10m_rotation_rmse_deg": final_metrics["rpe_10m"]["rotation"]["rmse_deg"],
                "raw_lidar_loc_ate_rmse_m": raw_metrics["ate_translation"]["rmse_m"],
                "raw_lidar_loc_rpe_10m_translation_rmse_m": raw_metrics["rpe_10m"]["translation"]["rmse_m"],
                "localization_valid_ratio": float(localization["valid_ratio"]),
                "processing_ms_mean": float(localization["processing_ms"]["mean"]),
                "processing_ms_median": float(localization["processing_ms"]["median"]),
                "processing_ms_p95": float(localization["processing_ms"]["p95"]),
                "processing_ms_max": float(localization["processing_ms"]["max"]),
                "confidence_mean": float(localization["confidence"]["mean"]),
                "confidence_p95": float(localization["confidence"]["p95"]),
                "wall_time_s": wall_time,
                "realtime_factor": wall_time / sensor_duration,
                **resources,
            }
            run_rows.append(row)
            sequence_payload["map_variants"].setdefault(variant, []).append({
                **row,
                "final_metrics": final_metrics,
                "raw_lidar_localization_metrics": raw_metrics,
            })
            aligned_sample_rows.extend({
                "sequence": sequence,
                "map_variant": variant,
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

    numeric_fields = [
        "ate_rmse_m", "ate_median_m", "ate_p95_m", "ate_max_m", "ate_rotation_rmse_deg",
        "rpe_1m_translation_rmse_m", "rpe_1m_rotation_rmse_deg",
        "rpe_10m_translation_rmse_m", "rpe_10m_rotation_rmse_deg",
        "raw_lidar_loc_ate_rmse_m", "raw_lidar_loc_rpe_10m_translation_rmse_m",
        "localization_valid_ratio", "processing_ms_mean", "processing_ms_median", "processing_ms_p95",
        "processing_ms_max", "confidence_mean", "confidence_p95", "wall_time_s", "realtime_factor",
        "mean_cpu_cores", "peak_cpu_cores", "p95_cpu_cores", "mean_rss_mb", "peak_rss_mb",
    ]
    summary_rows: list[dict[str, object]] = []
    for sequence in sequences:
        for variant in variants:
            group = [row for row in run_rows if row["sequence"] == sequence and row["map_variant"] == variant]
            summary: dict[str, object] = {
                "sequence": sequence, "map_variant": variant, "successful_repeats": len(group)
            }
            for field in numeric_fields:
                values = [float(item[field]) for item in group if item.get(field) is not None]
                mean, std = mean_std(values)
                summary[f"{field}_mean"] = mean
                summary[f"{field}_std"] = std
                if field in {"localization_valid_ratio", "confidence_mean", "confidence_p95"}:
                    summary[f"{field}_worst"] = min(values, default=None)
                else:
                    summary[f"{field}_worst"] = max(values, default=None)
            summary_rows.append(summary)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "run_metrics.csv", run_rows)
    write_csv(args.output_dir / "summary_metrics.csv", summary_rows)
    write_csv(args.output_dir / "aligned_samples.csv", aligned_sample_rows)
    (args.output_dir / "metrics.json").write_text(
        json.dumps(
            {"sequences": result_sequences, "failures": failures, "warnings": warnings},
            ensure_ascii=False,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    run_manifest = args.runs_root / "_state" / "experiment_manifest.json"
    analysis_manifest = {
        "schema_version": 1,
        "evaluator": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
        "run_manifest": {
            "path": str(run_manifest.resolve()),
            "sha256": sha256(run_manifest) if run_manifest.is_file() else None,
        },
        "ground_truth": {
            sequence: {
                "path": str(GT_PATHS[sequence].resolve()),
                "sha256": sha256(GT_PATHS[sequence]),
                "preprocessing": gt_preprocessing.get(sequence),
            }
            for sequence in sequences
        },
        "parameters": {
            "sequences": sequences, "map_variants": variants, "repeats": args.repeats,
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
        "required_run_count": len(sequences) * len(variants) * args.repeats,
        "evaluated_run_count": len(run_rows),
        "failures": failures,
        "warnings": warnings,
        "alignment": "fixed-scale SE(3)",
        "ground_truth_orientation_contract": (
            "all-identity RTK quaternions are placeholders: rotation metrics N/A; "
            "translation RPE uses aligned global displacement vectors"
        ),
        "common_support": "intersection of final/raw localization support across both map sources and all repeats",
    }
    (args.output_dir / "validation.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(validation, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
