#!/usr/bin/env python3
"""Evaluate the four-sequence M3DGR frontend matrix on common GT support."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation, Slerp


# 正式实验矩阵及其官方 RTK 文件；命令行参数可选取其中的子集。
SEQUENCES = ("Grass02", "Outdoor04", "Z-Rough-Road01", "Dark01")
METHODS = ("lightning_lm", "fastlio", "fastlivo2_lio", "voxel_slam_frontend")
GT_PATHS = {
    "Grass02": Path(r"F:\datasets\M3DGR\GT\Wheel_Slippage\Outdoor\Grass02.txt"),
    "Outdoor04": Path(r"F:\datasets\M3DGR\GT\Standard\Outdoor04.txt"),
    "Z-Rough-Road01": Path(r"F:\datasets\M3DGR\GT\Wheel_Slippage\Outdoor\Z-Rough-Road01.txt"),
    "Dark01": Path(r"F:\datasets\M3DGR\GT\Visual_Challenge\Outdoor\Dark01.txt"),
}


# 命令行入口：定义运行目录、分析输出目录和评测口径。
def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "m3dgr_frontend",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "analysis" / "m3dgr_frontend",
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


# 轨迹读取与真值预处理：严格检查估计轨迹；官方 RTK 则允许先读入再排序去重。
def load_tum(path: Path, *, require_monotonic: bool = True) -> np.ndarray:
    rows = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        values = [float(value) for value in line.split()]
        if len(values) != 8 or not np.all(np.isfinite(values)):
            raise ValueError(f"{path}:{line_number}: invalid TUM row")
        rows.append(values)
    poses = np.asarray(rows, dtype=np.float64)
    if len(poses) < 2:
        raise ValueError(f"{path}: missing poses")
    if require_monotonic and np.any(np.diff(poses[:, 0]) <= 0):
        raise ValueError(f"{path}: non-monotonic poses")
    norms = np.linalg.norm(poses[:, 4:8], axis=1)
    if np.any(np.abs(norms - 1.0) > 1e-3):
        raise ValueError(f"{path}: invalid quaternion")
    poses[:, 4:8] /= norms[:, None]
    return poses


def normalize_gt(poses: np.ndarray) -> tuple[np.ndarray, dict[str, int | str]]:
    original_deltas = np.diff(poses[:, 0])
    order = np.argsort(poses[:, 0], kind="stable")
    sorted_poses = poses[order]
    sorted_deltas = np.diff(sorted_poses[:, 0])
    keep = np.concatenate(([True], sorted_deltas > 0))
    normalized = sorted_poses[keep]
    if len(normalized) < 2 or np.any(np.diff(normalized[:, 0]) <= 0):
        raise ValueError("ground truth normalization did not produce strictly increasing timestamps")
    diagnostics: dict[str, int | str] = {
        "raw_pose_count": int(len(poses)),
        "normalized_pose_count": int(len(normalized)),
        "raw_duplicate_transition_count": int(np.count_nonzero(original_deltas == 0)),
        "raw_negative_transition_count": int(np.count_nonzero(original_deltas < 0)),
        "stable_sort_moved_row_count": int(np.count_nonzero(order != np.arange(len(poses)))),
        "duplicate_timestamp_rows_removed": int(len(sorted_poses) - len(normalized)),
        "policy": "stable sort by timestamp; retain the first row at each exact timestamp",
    }
    return normalized, diagnostics


def gt_orientation_available(poses: np.ndarray) -> bool:
    """Return false for M3DGR RTK files whose quaternions are identity placeholders."""
    identity = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float64)
    return not bool(np.all(np.abs(poses[:, 4:8] @ identity) >= 1.0 - 1e-12))


# 公共时间支持、位姿插值与固定尺度 SE(3) 对齐。
def support_mask(source: np.ndarray, times: np.ndarray, max_gap: float) -> np.ndarray:
    right = np.searchsorted(source[:, 0], times, side="left")
    valid = (right > 0) & (right < len(source))
    left = np.clip(right - 1, 0, len(source) - 1)
    right = np.clip(right, 0, len(source) - 1)
    valid &= source[right, 0] - source[left, 0] <= max_gap
    return valid


def interpolate(source: np.ndarray, times: np.ndarray) -> tuple[np.ndarray, Rotation]:
    positions = np.column_stack([
        np.interp(times, source[:, 0], source[:, column]) for column in range(1, 4)
    ])
    rotations = Slerp(source[:, 0], Rotation.from_quat(source[:, 4:8]))(times)
    return positions, rotations


def align_se3(estimate: np.ndarray, truth: np.ndarray) -> tuple[Rotation, np.ndarray]:
    estimate_center = estimate.mean(axis=0)
    truth_center = truth.mean(axis=0)
    u, _, vt = np.linalg.svd((estimate - estimate_center).T @ (truth - truth_center))
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(vt.T @ u.T)
    matrix = vt.T @ correction @ u.T
    return Rotation.from_matrix(matrix), truth_center - matrix @ estimate_center


# 精度统计：ATE 使用逐姿态误差，RPE 使用 1 m/10 m 路径间隔的相对运动误差。
def stats(values: np.ndarray, suffix: str) -> dict[str, float | int]:
    return {
        "count": int(len(values)),
        f"rmse_{suffix}": float(np.sqrt(np.mean(values * values))),
        f"mean_{suffix}": float(np.mean(values)),
        f"median_{suffix}": float(np.median(values)),
        f"p95_{suffix}": float(np.quantile(values, 0.95)),
        f"max_{suffix}": float(np.max(values)),
    }


def unavailable_stats(suffix: str) -> dict[str, float | int | None]:
    return {
        "count": 0,
        f"rmse_{suffix}": None,
        f"mean_{suffix}": None,
        f"median_{suffix}": None,
        f"p95_{suffix}": None,
        f"max_{suffix}": None,
    }


def rpe(
    estimate_position: np.ndarray,
    estimate_rotation: Rotation,
    truth_position: np.ndarray,
    truth_rotation: Rotation,
    distance: float,
    orientation_available: bool = True,
) -> dict[str, object]:
    cumulative = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(truth_position, axis=0), axis=1))))
    end = np.searchsorted(cumulative, cumulative + distance, side="left")
    start = np.arange(len(cumulative))
    valid = end < len(cumulative)
    start, end = start[valid], end[valid]
    if orientation_available:
        truth_rel_rotation = truth_rotation[start].inv() * truth_rotation[end]
        estimate_rel_rotation = estimate_rotation[start].inv() * estimate_rotation[end]
        truth_rel_translation = truth_rotation[start].inv().apply(truth_position[end] - truth_position[start])
        estimate_rel_translation = estimate_rotation[start].inv().apply(estimate_position[end] - estimate_position[start])
        translation_error = np.linalg.norm(estimate_rel_translation - truth_rel_translation, axis=1)
        rotation = stats(
            np.degrees((truth_rel_rotation.inv() * estimate_rel_rotation).magnitude()), "deg"
        )
    else:
        # Positions are already rigidly aligned.  Without GT attitude, compare
        # global displacement vectors instead of inventing body frames.
        translation_error = np.linalg.norm(
            (estimate_position[end] - estimate_position[start])
            - (truth_position[end] - truth_position[start]),
            axis=1,
        )
        rotation = unavailable_stats("deg")
    return {
        "distance_m": distance,
        "translation": stats(translation_error, "m"),
        "rotation": rotation,
    }


# 运行证据与通用输出：读取控制器元数据/资源监控，写 CSV，并计算材料哈希。
def parse_metadata(path: Path) -> dict[str, str]:
    values = {}
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    return values


def read_resources(path: Path) -> dict[str, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {
        key: float(payload[key])
        for key in ("mean_cpu_cores", "peak_cpu_cores", "p95_cpu_cores", "mean_rss_mb", "peak_rss_mb")
        if payload.get(key) is not None
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def mean_std(values: list[float]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    array = np.asarray(values, dtype=np.float64)
    return float(np.mean(array)), float(np.std(array, ddof=1)) if len(array) > 1 else 0.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# 完成性合同：Lightning 与 ROS1 基线采用各自的正常结束判据，并保留诊断警告。
def run_contract(metadata: dict[str, str], method: str) -> tuple[bool, list[str]]:
    base = (
        metadata.get("completion") == "reached_final_lidar"
    )
    if method == "lightning_lm":
        passed = (
            base
            and metadata.get("algorithm_rc") == "0"
            and metadata.get("watchdog_status") == "completed"
            and metadata.get("invalid_count") == "0"
            and metadata.get("nonmonotonic_count") == "0"
            and metadata.get("excessive_output_gap_count") == "0"
        )
        return passed, []
    passed = (
        base
        and metadata.get("rosbag_play_rc") == "0"
        and metadata.get("algorithm_contract") == "passed"
        and metadata.get("subscriber_contract") == "passed"
    )
    warnings = []
    warning_fields = {
        "invalid_count": "invalid recorder observations",
        "nonmonotonic_count": "non-monotonic recorder observations",
        "excessive_output_gap_count": "output gaps above the 0.2 s diagnostic threshold",
        "shutdown_abnormal_exit_count": "abnormal shutdown log observations",
    }
    for key, label in warning_fields.items():
        if int(metadata.get(key, "0")) > 0:
            warnings.append(f"{label}: {metadata[key]}")
    if metadata.get("algorithm_launch_rc") != "0":
        warnings.append(
            f"roslaunch return code after deliberate process-group shutdown: {metadata.get('algorithm_launch_rc')}"
        )
    return passed, warnings


def main() -> int:
    args = parse_args()
    sequences = tuple(item for item in args.sequences.split(",") if item)
    methods = tuple(item for item in args.methods.split(",") if item)
    if set(sequences) - set(SEQUENCES) or set(methods) - set(METHODS):
        raise SystemExit("unknown sequence or method")
    inventory_rows = json.loads(args.inventory.read_text(encoding="utf-8"))
    inventory = {row["sequence"]: row for row in inventory_rows}
    missing_inventory = sorted(set(sequences) - set(inventory))
    if missing_inventory:
        raise SystemExit(f"missing inventory rows: {missing_inventory}")

    # 阶段 1：收集每个序列/方法/重复的轨迹，先执行文件和完成性合同检查。
    all_runs: dict[str, dict[tuple[str, int], tuple[Path, np.ndarray, dict[str, str], list[str]]]] = {}
    failures = []
    warnings = []
    provenance_correction = args.runs_root / "_state" / "provenance_correction.json"
    if provenance_correction.is_file():
        warnings.append(
            "provenance correction: the four M3DGR sequences are development-exposed; "
            "this is a frozen confirmatory rerun, not an unseen holdout evaluation"
        )
    for sequence in sequences:
        sequence_runs = {}
        for method in methods:
            for repeat in range(1, args.repeats + 1):
                run_dir = args.runs_root / sequence / method / f"repeat_{repeat:02d}"
                metadata_path = run_dir / "run_metadata.txt"
                trajectory = run_dir / "trajectory_mid360.tum"
                resources = run_dir / "resource_summary.json"
                if not all(path.is_file() for path in (metadata_path, trajectory, resources)):
                    failures.append(f"missing artifacts: {sequence}/{method}/repeat_{repeat:02d}")
                    continue
                metadata = parse_metadata(metadata_path)
                passed, run_warnings = run_contract(metadata, method)
                if not passed:
                    failures.append(f"failed run contract: {sequence}/{method}/repeat_{repeat:02d}")
                    continue
                for warning in run_warnings:
                    warnings.append(f"{sequence}/{method}/repeat_{repeat:02d}: {warning}")
                sequence_runs[(method, repeat)] = (run_dir, load_tum(trajectory), metadata, run_warnings)
        all_runs[sequence] = sequence_runs

    # 阶段 2：在同一序列所有正式运行的共同 GT 时间支持上计算逐次和逐帧指标。
    run_rows = []
    aligned_sample_rows: list[dict[str, object]] = []
    result_sequences = {}
    gt_preprocessing: dict[str, dict[str, int | str]] = {}
    for sequence in sequences:
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
        sequence_runs = all_runs[sequence]
        if not all((method, repeat) in sequence_runs for method in methods for repeat in range(1, args.repeats + 1)):
            failures.append(f"incomplete common-support matrix: {sequence}")
            continue
        common = np.ones(len(gt), dtype=bool)
        for _, trajectory, _, _ in sequence_runs.values():
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
        sequence_payload = {
            "common_gt_pose_count": len(common_gt),
            "gt_orientation_available": orientation_available,
            "methods": {},
        }
        for (method, repeat), (run_dir, trajectory, metadata, run_warnings) in sequence_runs.items():
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
            expected_lidar_frames = int(inventory[sequence]["lidar_count"])
            output_deficit_ratio = max(0.0, 1.0 - output_frames / expected_lidar_frames)
            payload = {
                "sequence": sequence, "method": method, "repeat": repeat,
                "common_gt_pose_count": len(common_gt),
                "common_gt_coverage": len(common_gt) / len(gt),
                "gt_orientation_available": orientation_available,
                "ate_translation": ate_translation,
                "ate_rotation": ate_rotation,
                "rpe_1m": rpe_1m, "rpe_10m": rpe_10m,
                "wall_time_s": wall_time,
                "sensor_duration_s": sensor_duration,
                "realtime_factor": wall_time / sensor_duration,
                "output_frames": output_frames,
                "expected_lidar_frames": expected_lidar_frames,
                "output_deficit_ratio": output_deficit_ratio,
                "wall_ms_per_output_frame": wall_time * 1000.0 / output_frames,
                "maximum_output_gap_s": float(metadata["maximum_output_gap_s"]),
                "contract_warnings": run_warnings,
                "resources": resources,
            }
            sequence_payload["methods"].setdefault(method, []).append(payload)
            run_rows.append({
                "sequence": sequence, "method": method, "repeat": repeat,
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
                "wall_time_s": wall_time, "realtime_factor": wall_time / sensor_duration,
                "output_frames": output_frames,
                "expected_lidar_frames": expected_lidar_frames,
                "output_deficit_ratio": output_deficit_ratio,
                "wall_ms_per_output_frame": wall_time * 1000.0 / output_frames,
                "maximum_output_gap_s": float(metadata["maximum_output_gap_s"]),
                "invalid_recorder_observations": int(metadata.get("invalid_count", "0")),
                "nonmonotonic_recorder_observations": int(metadata.get("nonmonotonic_count", "0")),
                "excessive_output_gap_count": int(metadata.get("excessive_output_gap_count", "0")),
                "shutdown_abnormal_exit_count": int(metadata.get("shutdown_abnormal_exit_count", "0")),
                "algorithm_launch_rc": int(metadata.get("algorithm_launch_rc", metadata.get("algorithm_rc", "0"))),
                "contract_warning_count": len(run_warnings),
                **resources,
            })
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

    # 阶段 3：按“序列 × 方法”汇总三次重复的均值、样本标准差和最差值。
    summary_rows = []
    numeric_fields = [
        "ate_rmse_m", "ate_median_m", "ate_p95_m", "ate_max_m", "ate_rotation_rmse_deg",
        "rpe_1m_translation_rmse_m", "rpe_1m_rotation_rmse_deg",
        "rpe_10m_translation_rmse_m", "rpe_10m_rotation_rmse_deg",
        "wall_time_s", "realtime_factor", "output_deficit_ratio", "wall_ms_per_output_frame",
        "maximum_output_gap_s", "excessive_output_gap_count", "nonmonotonic_recorder_observations",
        "shutdown_abnormal_exit_count", "contract_warning_count", "mean_cpu_cores", "peak_cpu_cores",
        "p95_cpu_cores", "mean_rss_mb", "peak_rss_mb",
    ]
    for sequence in sequences:
        for method in methods:
            group = [row for row in run_rows if row["sequence"] == sequence and row["method"] == method]
            row: dict[str, object] = {"sequence": sequence, "method": method, "successful_repeats": len(group)}
            for field in numeric_fields:
                mean, std = mean_std([float(item[field]) for item in group if item.get(field) is not None])
                row[f"{field}_mean"] = mean
                row[f"{field}_std"] = std
                row[f"{field}_worst"] = max(
                    (float(item[field]) for item in group if item.get(field) is not None), default=None
                )
            summary_rows.append(row)

    # 阶段 4：写出分析表、机器可读指标、来源清单和最终验证状态。
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "run_metrics.csv", run_rows)
    write_csv(args.output_dir / "summary_metrics.csv", summary_rows)
    write_csv(args.output_dir / "aligned_samples.csv", aligned_sample_rows)
    (args.output_dir / "metrics.json").write_text(
        json.dumps({"sequences": result_sequences, "failures": failures, "warnings": warnings}, ensure_ascii=False, indent=2) + "\n",
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
        "provenance_correction": {
            "path": str(provenance_correction.resolve()),
            "sha256": sha256(provenance_correction) if provenance_correction.is_file() else None,
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
        "common_support": "intersection of interpolable official-GT timestamps across every required run per sequence",
    }
    (args.output_dir / "validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
