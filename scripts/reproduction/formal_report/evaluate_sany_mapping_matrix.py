#!/usr/bin/env python3
"""Evaluate formal SANY 20260701 single/four-LiDAR mapping trials."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from evaluate_m3dgr_frontend_matrix import (
    align_se3,
    interpolate,
    mean_std,
    normalize_gt,
    parse_metadata,
    read_resources,
    rpe,
    sha256,
    stats,
    support_mask,
    write_csv,
)


VARIANTS = ("single114", "four_lidar", "four_no_noise", "drop114", "drop127", "drop187", "drop195")
PIPELINES = ("frontend", "slam")
SLAM_VARIANTS = ("single114", "four_lidar")


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "sany_20260701_mapping",
    )
    parser.add_argument(
        "--truth", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "sany_voxel114_reference"
        / "data_20260701" / "results" / "trajectory_voxel_opt.tum",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "analysis" / "sany_20260701_mapping",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--max-interpolation-gap", type=float, default=0.2)
    return parser.parse_args()


def load_first8(path: Path) -> np.ndarray:
    rows = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 8:
            raise ValueError(f"{path}:{line_number}: fewer than 8 columns")
        values = [float(value) for value in fields[:8]]
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{path}:{line_number}: non-finite pose")
        rows.append(values)
    poses = np.asarray(rows, dtype=np.float64)
    if len(poses) < 2 or np.any(np.diff(poses[:, 0]) <= 0):
        raise ValueError(f"{path}: missing or non-monotonic poses")
    norms = np.linalg.norm(poses[:, 4:8], axis=1)
    if np.any(np.abs(norms - 1.0) > 1e-3):
        raise ValueError(f"{path}: invalid quaternion")
    poses[:, 4:8] /= norms[:, None]
    return poses


def trajectory_path(run_dir: Path, pipeline: str) -> Path:
    if pipeline == "frontend":
        return run_dir / "results" / "trajectory_lidar114.tum"
    return run_dir / "results" / "trajectory_slam_opt.tum"


def run_passed(metadata: dict[str, str], pipeline: str) -> bool:
    expected_method = "lightning_lm" if pipeline == "frontend" else "lightning_lm_offline_slam_map_export"
    return (
        metadata.get("method") == expected_method
        and metadata.get("completion") == "reached_final_lidar"
        and metadata.get("algorithm_rc") == "0"
        and metadata.get("watchdog_status") == "completed"
        and metadata.get("invalid_count") == "0"
        and metadata.get("nonmonotonic_count") == "0"
        and metadata.get("excessive_output_gap_count") == "0"
    )


def frame_contract(path: Path) -> dict[str, float | int | None]:
    if not path.is_file() or path.stat().st_size == 0:
        return {"assembled_frames": 0, "partial_frames": 0, "partial_ratio": None}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    partial = sum(int(row["partial"]) for row in rows)
    return {
        "assembled_frames": len(rows),
        "partial_frames": partial,
        "partial_ratio": partial / len(rows) if rows else None,
    }


def main() -> int:
    args = parse_args()
    truth, truth_normalization = normalize_gt(load_first8(args.truth))
    runs: dict[tuple[str, str, int], tuple[Path, np.ndarray, dict[str, str]]] = {}
    failures: list[str] = []
    cells = [(variant, "frontend") for variant in VARIANTS] + [(variant, "slam") for variant in SLAM_VARIANTS]
    for variant, pipeline in cells:
        for repeat in range(1, args.repeats + 1):
            run_dir = args.runs_root / variant / pipeline / f"repeat_{repeat:02d}"
            metadata_path = run_dir / "run_metadata.txt"
            trajectory = trajectory_path(run_dir, pipeline)
            resources = run_dir / "resource_summary.json"
            if not all(path.is_file() for path in (metadata_path, trajectory, resources)):
                failures.append(f"missing artifacts: {variant}/{pipeline}/repeat_{repeat:02d}")
                continue
            metadata = parse_metadata(metadata_path)
            if not run_passed(metadata, pipeline):
                failures.append(f"failed run contract: {variant}/{pipeline}/repeat_{repeat:02d}")
                continue
            runs[(variant, pipeline, repeat)] = (run_dir, load_first8(trajectory), metadata)

    required = [(variant, pipeline, repeat) for variant, pipeline in cells
                for repeat in range(1, args.repeats + 1)]
    common = np.ones(len(truth), dtype=bool)
    if all(key in runs for key in required):
        for _, trajectory, _ in runs.values():
            common &= support_mask(trajectory, truth[:, 0], args.max_interpolation_gap)
    else:
        failures.append("incomplete common-support matrix")
        common[:] = False
    common_truth = truth[common]
    if len(common_truth) < 100:
        failures.append(f"insufficient common proxy-truth support: {len(common_truth)}")

    run_rows: list[dict[str, object]] = []
    aligned_sample_rows: list[dict[str, object]] = []
    metrics_payload: dict[str, object] = {}
    if len(common_truth) >= 100:
        truth_position = common_truth[:, 1:4]
        truth_rotation = Rotation.from_quat(common_truth[:, 4:8])
        for (variant, pipeline, repeat), (run_dir, trajectory, metadata) in runs.items():
            position, rotation = interpolate(trajectory, common_truth[:, 0])
            alignment_rotation, alignment_translation = align_se3(position, truth_position)
            position = alignment_rotation.apply(position) + alignment_translation
            rotation = alignment_rotation * rotation
            translation_error = np.linalg.norm(position - truth_position, axis=1)
            rotation_error = np.degrees((truth_rotation.inv() * rotation).magnitude())
            ate_translation = stats(translation_error, "m")
            ate_rotation = stats(rotation_error, "deg")
            rpe_1m = rpe(position, rotation, truth_position, truth_rotation, 1.0)
            rpe_10m = rpe(position, rotation, truth_position, truth_rotation, 10.0)
            resources = read_resources(run_dir / "resource_summary.json")
            wall_time = float(metadata["wall_time_s"])
            sensor_duration = float(metadata["sensor_duration_s"])
            frame_stats = frame_contract(run_dir / "results" / "frame_stats.csv")
            row = {
                "variant": variant,
                "pipeline": pipeline,
                "repeat": repeat,
                "common_truth_pose_count": len(common_truth),
                "common_truth_coverage": len(common_truth) / len(truth),
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
                "output_frames": int(metadata["trajectory_lines"]),
                "wall_ms_per_output_frame": wall_time * 1000.0 / int(metadata["trajectory_lines"]),
                **frame_stats,
                **resources,
            }
            run_rows.append(row)
            metrics_payload[f"{variant}/{pipeline}/repeat_{repeat:02d}"] = {
                **row,
                "ate_translation": ate_translation,
                "ate_rotation": ate_rotation,
                "rpe_1m": rpe_1m,
                "rpe_10m": rpe_10m,
            }
            aligned_sample_rows.extend({
                "variant": variant,
                "pipeline": pipeline,
                "repeat": repeat,
                "timestamp_s": float(timestamp),
                "truth_x_m": float(truth_sample[0]),
                "truth_y_m": float(truth_sample[1]),
                "truth_z_m": float(truth_sample[2]),
                "estimate_x_m": float(estimate[0]),
                "estimate_y_m": float(estimate[1]),
                "estimate_z_m": float(estimate[2]),
                "translation_error_m": float(error),
                "rotation_error_deg": float(rotation_error_value),
            } for timestamp, truth_sample, estimate, error, rotation_error_value in zip(
                common_truth[:, 0], truth_position, position, translation_error, rotation_error, strict=True
            ))

    numeric_fields = [
        "ate_rmse_m", "ate_median_m", "ate_p95_m", "ate_max_m", "ate_rotation_rmse_deg",
        "rpe_1m_translation_rmse_m", "rpe_1m_rotation_rmse_deg",
        "rpe_10m_translation_rmse_m", "rpe_10m_rotation_rmse_deg",
        "wall_time_s", "realtime_factor", "wall_ms_per_output_frame", "partial_ratio",
        "mean_cpu_cores", "peak_cpu_cores", "p95_cpu_cores", "mean_rss_mb", "peak_rss_mb",
    ]
    summary_rows: list[dict[str, object]] = []
    for variant, pipeline in cells:
        group = [row for row in run_rows if row["variant"] == variant and row["pipeline"] == pipeline]
        summary: dict[str, object] = {
            "variant": variant, "pipeline": pipeline, "successful_repeats": len(group)
        }
        for field in numeric_fields:
            values = [float(item[field]) for item in group if item.get(field) is not None]
            mean, std = mean_std(values)
            summary[f"{field}_mean"] = mean
            summary[f"{field}_std"] = std
            summary[f"{field}_worst"] = max(values, default=None)
        summary_rows.append(summary)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "run_metrics.csv", run_rows)
    write_csv(args.output_dir / "summary_metrics.csv", summary_rows)
    write_csv(args.output_dir / "aligned_samples.csv", aligned_sample_rows)
    (args.output_dir / "metrics.json").write_text(
        json.dumps({"runs": metrics_payload, "failures": failures}, ensure_ascii=False, indent=2) + "\n",
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
        "proxy_reference": {"path": str(args.truth.resolve()), "sha256": sha256(args.truth)},
        "parameters": {
            "variants": VARIANTS, "slam_variants": SLAM_VARIANTS, "repeats": args.repeats,
            "max_interpolation_gap_s": args.max_interpolation_gap,
            "alignment": "fixed-scale SE(3)",
        },
    }
    (args.output_dir / "analysis_manifest.json").write_text(
        json.dumps(analysis_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    validation = {
        "status": "passed" if not failures else "failed",
        "required_run_count": len(required),
        "evaluated_run_count": len(run_rows),
        "failures": failures,
        "reference": "Voxel-SLAM 114 trajectory; cross-implementation proxy, not absolute ground truth",
        "reference_normalization": truth_normalization,
        "alignment": "fixed-scale SE(3)",
        "common_support": "intersection across both sensor variants, both pipelines and all repeats",
    }
    (args.output_dir / "validation.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(validation, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
