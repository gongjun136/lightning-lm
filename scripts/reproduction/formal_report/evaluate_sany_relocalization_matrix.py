#!/usr/bin/env python3
"""Evaluate SANY 20260716 startup and forced-loss relocalization trials."""

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
from evaluate_sany_mapping_matrix import load_first8


DATASETS = ("data1", "data2")
MODES = ("startup", "forced_loss")
FORCED_FRAMES = {"data1": 600, "data2": 350}


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "sany_20260716_relocalization",
    )
    parser.add_argument(
        "--reference-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "sany_voxel114_reference",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "analysis" / "sany_20260716_relocalization",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--max-interpolation-gap", type=float, default=0.2)
    return parser.parse_args()


def truth_path(root: Path, dataset: str) -> Path:
    return root / dataset / "results" / "trajectory_voxel_opt.tum"


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


def event_metrics(path: Path, dataset: str, mode: str) -> dict[str, float | int | None]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    attempts = [row for row in rows if int(row["relocalization_attempted"]) == 1]
    accepts = [row for row in rows if int(row["relocalization_accepted"]) == 1]
    first_timestamp = float(rows[0]["timestamp"])
    startup_accept = accepts[0] if accepts else None
    payload: dict[str, float | int | None] = {
        "relocalization_attempts": len(attempts),
        "relocalization_accepts": len(accepts),
        "startup_accept_frame": int(startup_accept["frame_index"]) if startup_accept else None,
        "startup_latency_s": float(startup_accept["timestamp"]) - first_timestamp if startup_accept else None,
        "recovery_accept_frame": None,
        "recovery_attempts": None,
        "recovery_latency_s": None,
    }
    if mode != "forced_loss":
        return payload
    forced_frame = FORCED_FRAMES[dataset]
    injection = next((row for row in rows if int(row["frame_index"]) == forced_frame), None)
    recovery_accepts = [row for row in accepts if int(row["frame_index"]) >= forced_frame]
    if injection and recovery_accepts:
        accepted = recovery_accepts[0]
        payload["recovery_accept_frame"] = int(accepted["frame_index"])
        payload["recovery_attempts"] = sum(
            int(row["relocalization_attempted"]) for row in rows
            if forced_frame <= int(row["frame_index"]) <= int(accepted["frame_index"])
        )
        payload["recovery_latency_s"] = float(accepted["timestamp"]) - float(injection["timestamp"])
    return payload


def main() -> int:
    args = parse_args()
    failures: list[str] = []
    runs: dict[str, dict[tuple[str, int], tuple[Path, np.ndarray, dict[str, str], dict]]] = {}
    for dataset in DATASETS:
        dataset_runs = {}
        for mode in MODES:
            for repeat in range(1, args.repeats + 1):
                run_dir = args.runs_root / dataset / mode / f"repeat_{repeat:02d}"
                required = [
                    run_dir / "run_metadata.txt",
                    run_dir / "resource_summary.json",
                    run_dir / "results" / "trajectory_loc.tum",
                    run_dir / "results" / "localization_summary.json",
                    run_dir / "results" / "localization_stats.csv",
                ]
                if not all(path.is_file() for path in required):
                    failures.append(f"missing artifacts: {dataset}/{mode}/repeat_{repeat:02d}")
                    continue
                metadata = parse_metadata(required[0])
                if not run_passed(metadata):
                    failures.append(f"failed run contract: {dataset}/{mode}/repeat_{repeat:02d}")
                    continue
                summary = json.loads(required[3].read_text(encoding="utf-8"))
                dataset_runs[(mode, repeat)] = (run_dir, load_first8(required[2]), metadata, summary)
        runs[dataset] = dataset_runs

    run_rows: list[dict[str, object]] = []
    aligned_sample_rows: list[dict[str, object]] = []
    payload: dict[str, object] = {}
    truth_normalization: dict[str, dict[str, int | str]] = {}
    for dataset in DATASETS:
        dataset_runs = runs[dataset]
        required_keys = [(mode, repeat) for mode in MODES for repeat in range(1, args.repeats + 1)]
        if not all(key in dataset_runs for key in required_keys):
            failures.append(f"incomplete common-support matrix: {dataset}")
            continue
        truth, truth_normalization[dataset] = normalize_gt(
            load_first8(truth_path(args.reference_root, dataset))
        )
        common = np.ones(len(truth), dtype=bool)
        for _, trajectory, _, _ in dataset_runs.values():
            common &= support_mask(trajectory, truth[:, 0], args.max_interpolation_gap)
        common_truth = truth[common]
        if len(common_truth) < 100:
            failures.append(f"insufficient common proxy-truth support: {dataset} ({len(common_truth)})")
            continue
        truth_position = common_truth[:, 1:4]
        truth_rotation = Rotation.from_quat(common_truth[:, 4:8])
        for (mode, repeat), (run_dir, trajectory, metadata, summary) in dataset_runs.items():
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
            early_mask = common_truth[:, 0] <= common_truth[0, 0] + 30.0
            early_error = np.linalg.norm(position[early_mask] - truth_position[early_mask], axis=1)
            localization = summary["localization"]
            events = event_metrics(run_dir / "results" / "localization_stats.csv", dataset, mode)
            resources = read_resources(run_dir / "resource_summary.json")
            wall_time = float(metadata["wall_time_s"])
            sensor_duration = float(metadata["sensor_duration_s"])
            row = {
                "dataset": dataset,
                "mode": mode,
                "repeat": repeat,
                "common_truth_pose_count": len(common_truth),
                "common_truth_coverage": len(common_truth) / len(truth),
                "ate_rmse_m": ate_translation["rmse_m"],
                "ate_median_m": ate_translation["median_m"],
                "ate_p95_m": ate_translation["p95_m"],
                "ate_max_m": ate_translation["max_m"],
                "ate_rotation_rmse_deg": ate_rotation["rmse_deg"],
                "early_30s_ate_rmse_m": float(np.sqrt(np.mean(early_error * early_error))),
                "rpe_1m_translation_rmse_m": rpe_1m["translation"]["rmse_m"],
                "rpe_1m_rotation_rmse_deg": rpe_1m["rotation"]["rmse_deg"],
                "rpe_10m_translation_rmse_m": rpe_10m["translation"]["rmse_m"],
                "rpe_10m_rotation_rmse_deg": rpe_10m["rotation"]["rmse_deg"],
                "localization_valid_ratio": float(localization["valid_ratio"]),
                "processing_ms_mean": float(localization["processing_ms"]["mean"]),
                "processing_ms_median": float(localization["processing_ms"]["median"]),
                "processing_ms_p95": float(localization["processing_ms"]["p95"]),
                "processing_ms_max": float(localization["processing_ms"]["max"]),
                "wall_time_s": wall_time,
                "realtime_factor": wall_time / sensor_duration,
                **events,
                **resources,
            }
            run_rows.append(row)
            payload[f"{dataset}/{mode}/repeat_{repeat:02d}"] = {
                **row,
                "ate_translation": ate_translation,
                "ate_rotation": ate_rotation,
                "rpe_1m": rpe_1m,
                "rpe_10m": rpe_10m,
            }
            aligned_sample_rows.extend({
                "dataset": dataset,
                "mode": mode,
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
        "early_30s_ate_rmse_m", "rpe_1m_translation_rmse_m", "rpe_1m_rotation_rmse_deg",
        "rpe_10m_translation_rmse_m", "rpe_10m_rotation_rmse_deg", "localization_valid_ratio",
        "processing_ms_mean", "processing_ms_median", "processing_ms_p95", "processing_ms_max",
        "startup_latency_s", "recovery_latency_s", "recovery_attempts", "wall_time_s", "realtime_factor",
        "mean_cpu_cores", "peak_cpu_cores", "p95_cpu_cores", "mean_rss_mb", "peak_rss_mb",
    ]
    summary_rows: list[dict[str, object]] = []
    for dataset in DATASETS:
        for mode in MODES:
            group = [row for row in run_rows if row["dataset"] == dataset and row["mode"] == mode]
            summary: dict[str, object] = {
                "dataset": dataset, "mode": mode, "successful_repeats": len(group)
            }
            for field in numeric_fields:
                values = [float(item[field]) for item in group if item.get(field) is not None]
                mean, std = mean_std(values)
                summary[f"{field}_mean"] = mean
                summary[f"{field}_std"] = std
                if field == "localization_valid_ratio":
                    summary[f"{field}_worst"] = min(values, default=None)
                else:
                    summary[f"{field}_worst"] = max(values, default=None)
            summary_rows.append(summary)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "run_metrics.csv", run_rows)
    write_csv(args.output_dir / "summary_metrics.csv", summary_rows)
    write_csv(args.output_dir / "aligned_samples.csv", aligned_sample_rows)
    (args.output_dir / "metrics.json").write_text(
        json.dumps({"runs": payload, "failures": failures}, ensure_ascii=False, indent=2) + "\n",
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
        "proxy_reference": {
            dataset: {
                "path": str(truth_path(args.reference_root, dataset).resolve()),
                "sha256": sha256(truth_path(args.reference_root, dataset)),
            }
            for dataset in DATASETS
        },
        "parameters": {
            "datasets": DATASETS, "modes": MODES, "forced_frames": FORCED_FRAMES,
            "repeats": args.repeats, "max_interpolation_gap_s": args.max_interpolation_gap,
            "alignment": "fixed-scale SE(3)",
        },
    }
    (args.output_dir / "analysis_manifest.json").write_text(
        json.dumps(analysis_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    validation = {
        "status": "passed" if not failures else "failed",
        "required_run_count": len(DATASETS) * len(MODES) * args.repeats,
        "evaluated_run_count": len(run_rows),
        "failures": failures,
        "reference": "corresponding Voxel-SLAM 114 optimized trajectory; proxy, not absolute ground truth",
        "reference_normalization": truth_normalization,
        "alignment": "fixed-scale SE(3)",
        "common_support": "intersection across startup/forced-loss trials and all repeats per dataset",
    }
    (args.output_dir / "validation.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(validation, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
