#!/usr/bin/env python3
"""Compare the beginning of a localization trajectory with a SLAM trajectory."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True, help="SLAM TUM trajectory")
    parser.add_argument("--estimate", type=Path, required=True, help="localization TUM trajectory")
    parser.add_argument("--output", type=Path, required=True, help="output JSON report")
    parser.add_argument("--window-seconds", type=float, default=30.0)
    parser.add_argument("--maximum-window-seconds", type=float, default=45.0)
    parser.add_argument("--minimum-window-distance", type=float, default=10.0)
    parser.add_argument("--motion-onset-window-seconds", type=float, default=2.0)
    parser.add_argument("--motion-onset-distance", type=float, default=0.5)
    parser.add_argument("--maximum-interpolation-gap", type=float, default=0.2)
    parser.add_argument("--minimum-associated-poses", type=int, default=30)
    parser.add_argument("--minimum-coverage", type=float, default=0.90)
    parser.add_argument("--ate-rmse-threshold", type=float, default=0.50)
    parser.add_argument("--rpe-1m-rmse-threshold", type=float, default=0.30)
    parser.add_argument("--rpe-10m-rmse-threshold", type=float, default=0.50)
    parser.add_argument("--no-plot", action="store_true", help="skip PNG generation")
    return parser.parse_args()


def load_tum(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8-sig") as stream:
        for line_number, raw in enumerate(stream, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 8:
                raise ValueError(f"{path}:{line_number}: expected 8 TUM fields")
            row = [float(value) for value in fields]
            if not np.all(np.isfinite(row)):
                raise ValueError(f"{path}:{line_number}: non-finite pose")
            rows.append(row)
    if not rows:
        raise ValueError(f"{path}: no poses")
    poses = np.asarray(rows, dtype=np.float64)
    poses = poses[np.argsort(poses[:, 0], kind="stable")]
    _, unique = np.unique(poses[:, 0], return_index=True)
    poses = poses[np.sort(unique)]
    if len(poses) < 3 or np.any(np.diff(poses[:, 0]) <= 0.0):
        raise ValueError(f"{path}: trajectory timestamps are invalid")
    return poses


def interpolate_positions(
    source: np.ndarray, timestamps: np.ndarray, maximum_gap: float
) -> tuple[np.ndarray, np.ndarray]:
    right = np.searchsorted(source[:, 0], timestamps, side="left")
    exact = (right < len(source)) & (source[np.minimum(right, len(source) - 1), 0] == timestamps)
    valid = exact | ((right > 0) & (right < len(source)))
    left = np.clip(right - 1, 0, len(source) - 1)
    right = np.clip(right, 0, len(source) - 1)
    gaps = source[right, 0] - source[left, 0]
    valid &= exact | (gaps <= maximum_gap)
    alpha = np.zeros_like(timestamps)
    interpolate = valid & ~exact & (gaps > 0.0)
    alpha[interpolate] = (
        timestamps[interpolate] - source[left[interpolate], 0]
    ) / gaps[interpolate]
    positions = source[left, 1:4] * (1.0 - alpha[:, None]) + source[right, 1:4] * alpha[:, None]
    positions[exact] = source[np.minimum(np.searchsorted(source[:, 0], timestamps[exact]), len(source) - 1), 1:4]
    positions[~valid] = np.nan
    return positions, valid


def align_se3(estimate: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    estimate_center = estimate.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (estimate - estimate_center).T @ (reference - reference_center)
    u, _, vt = np.linalg.svd(covariance)
    handedness = np.eye(3)
    handedness[2, 2] = np.linalg.det(vt.T @ u.T)
    rotation = vt.T @ handedness @ u.T
    translation = reference_center - rotation @ estimate_center
    return rotation, translation


def statistics(errors: np.ndarray) -> dict[str, float | int]:
    return {
        "count": int(len(errors)),
        "rmse_m": float(np.sqrt(np.mean(errors * errors))),
        "mean_m": float(np.mean(errors)),
        "median_m": float(np.median(errors)),
        "p95_m": float(np.percentile(errors, 95.0)),
        "max_m": float(np.max(errors)),
    }


def rpe_at_distance(
    estimate: np.ndarray, reference: np.ndarray, distance: float
) -> dict[str, float | int] | None:
    cumulative = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(reference, axis=0), axis=1))))
    end = np.searchsorted(cumulative, cumulative + distance, side="left")
    start = np.arange(len(reference))
    valid = end < len(reference)
    if not np.any(valid):
        return None
    error = np.linalg.norm(
        (estimate[end[valid]] - estimate[start[valid]])
        - (reference[end[valid]] - reference[start[valid]]),
        axis=1,
    )
    return statistics(error)


def select_window(
    reference: np.ndarray, estimate: np.ndarray, args: argparse.Namespace
) -> tuple[np.ndarray, bool]:
    start = max(reference[0, 0], estimate[0, 0])
    overlap_end = min(reference[-1, 0], estimate[-1, 0])
    candidates = reference[(reference[:, 0] >= start) & (reference[:, 0] <= overlap_end)]
    if len(candidates) < 2:
        return candidates, False
    motion_onset = None
    for index in range(len(candidates) - 1):
        window_end = np.searchsorted(
            candidates[:, 0], candidates[index, 0] + args.motion_onset_window_seconds, side="right"
        )
        if window_end <= index + 1:
            continue
        distance = float(
            np.sum(np.linalg.norm(np.diff(candidates[index:window_end, 1:4], axis=0), axis=1))
        )
        if distance >= args.motion_onset_distance:
            motion_onset = index
            break
    motion_detected = motion_onset is not None
    if motion_detected:
        candidates = candidates[motion_onset:]
        start = float(candidates[0, 0])
    hard_end = min(overlap_end, start + args.maximum_window_seconds)
    candidates = candidates[candidates[:, 0] <= hard_end]
    cumulative = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(candidates[:, 1:4], axis=0), axis=1))))
    nominal_end = min(hard_end, start + args.window_seconds)
    time_index = int(np.searchsorted(candidates[:, 0], nominal_end, side="right"))
    distance_index = int(np.searchsorted(cumulative, args.minimum_window_distance, side="left")) + 1
    count = min(len(candidates), max(time_index, distance_index))
    return candidates[:count], motion_detected


def save_plot(path: Path, timestamps: np.ndarray, reference: np.ndarray, aligned: np.ndarray, errors: np.ndarray) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    figure, axes = plt.subplots(1, 2, figsize=(13, 6))
    axes[0].plot(reference[:, 0], reference[:, 1], color="#222222", linewidth=2.2, label="SLAM reference")
    axes[0].plot(aligned[:, 0], aligned[:, 1], color="#e6550d", linewidth=1.5, label="localization aligned")
    axes[0].scatter(reference[0, 0], reference[0, 1], color="#238b45", s=45, label="window start")
    axes[0].set_aspect("equal", adjustable="datalim")
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("y (m)")
    axes[0].set_title("Early trajectory alignment (fixed-scale SE(3))")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend()
    axes[1].plot(timestamps - timestamps[0], errors, color="#756bb1", linewidth=1.5)
    axes[1].set_xlabel("time since relocalization window start (s)")
    axes[1].set_ylabel("translation error (m)")
    axes[1].set_title("Aligned translation error")
    axes[1].grid(True, alpha=0.25)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    positive = (
        args.window_seconds,
        args.maximum_window_seconds,
        args.minimum_window_distance,
        args.motion_onset_window_seconds,
        args.motion_onset_distance,
        args.maximum_interpolation_gap,
        args.minimum_coverage,
        args.ate_rmse_threshold,
        args.rpe_1m_rmse_threshold,
        args.rpe_10m_rmse_threshold,
    )
    if not all(math.isfinite(value) and value > 0.0 for value in positive):
        raise ValueError("window and threshold options must be finite and positive")
    if args.maximum_window_seconds < args.window_seconds:
        raise ValueError("maximum-window-seconds must be >= window-seconds")

    reference_all = load_tum(args.reference)
    estimate_all = load_tum(args.estimate)
    reference_window, motion_onset_detected = select_window(reference_all, estimate_all, args)
    estimate_interpolated, valid = interpolate_positions(
        estimate_all, reference_window[:, 0], args.maximum_interpolation_gap
    )
    if np.count_nonzero(valid) < 3:
        raise ValueError("fewer than three time-associated poses in the early window")
    reference = reference_window[valid, 1:4]
    estimate = estimate_interpolated[valid]
    timestamps = reference_window[valid, 0]
    rotation, translation = align_se3(estimate, reference)
    aligned = (rotation @ estimate.T).T + translation
    errors = np.linalg.norm(aligned - reference, axis=1)
    ate = statistics(errors)
    rpe_1m = rpe_at_distance(aligned, reference, 1.0)
    rpe_10m = rpe_at_distance(aligned, reference, 10.0)
    travelled_distance = float(np.sum(np.linalg.norm(np.diff(reference, axis=0), axis=1)))
    coverage = float(np.count_nonzero(valid) / len(reference_window))
    checks = {
        "minimum_associated_poses": int(np.count_nonzero(valid)) >= args.minimum_associated_poses,
        "motion_onset_detected": motion_onset_detected,
        "minimum_window_distance": travelled_distance >= args.minimum_window_distance,
        "coverage": coverage >= args.minimum_coverage,
        "ate_rmse": ate["rmse_m"] <= args.ate_rmse_threshold,
        "rpe_1m_rmse": rpe_1m is not None and rpe_1m["rmse_m"] <= args.rpe_1m_rmse_threshold,
        "rpe_10m_rmse": rpe_10m is not None and rpe_10m["rmse_m"] <= args.rpe_10m_rmse_threshold,
    }
    report = {
        "reference": str(args.reference.resolve()),
        "estimate": str(args.estimate.resolve()),
        "alignment": "fixed-scale SE(3) Horn/SVD; scale fixed to 1.0",
        "window_start_s": float(timestamps[0]),
        "window_end_s": float(timestamps[-1]),
        "window_duration_s": float(timestamps[-1] - timestamps[0]),
        "window_reference_distance_m": travelled_distance,
        "motion_onset_detected": motion_onset_detected,
        "associated_poses": int(np.count_nonzero(valid)),
        "coverage": coverage,
        "se3_rotation": rotation.tolist(),
        "se3_translation_m": translation.tolist(),
        "ate_translation": ate,
        "rpe_translation": {"1m": rpe_1m, "10m": rpe_10m},
        "thresholds": {
            "minimum_associated_poses": args.minimum_associated_poses,
            "minimum_window_distance_m": args.minimum_window_distance,
            "motion_onset_window_s": args.motion_onset_window_seconds,
            "motion_onset_distance_m": args.motion_onset_distance,
            "minimum_coverage": args.minimum_coverage,
            "ate_rmse_m": args.ate_rmse_threshold,
            "rpe_1m_rmse_m": args.rpe_1m_rmse_threshold,
            "rpe_10m_rmse_m": args.rpe_10m_rmse_threshold,
        },
        "checks": checks,
        "passed": all(checks.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not args.no_plot:
        save_plot(args.output.with_suffix(".png"), timestamps, reference, aligned, errors)
    print(args.output.resolve())
    print(f"passed={report['passed']} ate_rmse={ate['rmse_m']:.6f} distance={travelled_distance:.3f}")


if __name__ == "__main__":
    main()
