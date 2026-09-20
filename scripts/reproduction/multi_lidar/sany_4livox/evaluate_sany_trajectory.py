#!/usr/bin/env python3
"""Evaluate SANY trajectories against alidarState with fixed-scale SE(3)."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--run", action="append", required=True, metavar="LABEL=TUM")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-interpolation-gap", type=float, default=0.2)
    parser.add_argument("--ate-rmse-threshold", type=float, default=1.0)
    parser.add_argument("--rpe-10m-rmse-threshold", type=float, default=0.5)
    parser.add_argument("--minimum-coverage", type=float, default=0.95)
    return parser.parse_args()


def load_poses(path: Path, *, allow_extra_fields: bool) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8-sig") as stream:
        for line_number, raw in enumerate(stream, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            invalid_field_count = len(fields) < 8 if allow_extra_fields else len(fields) != 8
            if invalid_field_count:
                expected = "at least 8" if allow_extra_fields else "8"
                raise ValueError(f"{path}:{line_number}: expected {expected} pose fields")
            values = [float(value) for value in fields[:8]]
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{path}:{line_number}: non-finite pose")
            rows.append(values)
    if not rows:
        raise ValueError(f"{path}: no poses")
    poses = np.asarray(rows, dtype=np.float64)
    poses = poses[np.argsort(poses[:, 0], kind="stable")]
    _, unique_indices = np.unique(poses[:, 0], return_index=True)
    poses = poses[np.sort(unique_indices)]
    if np.any(np.diff(poses[:, 0]) <= 0.0):
        raise ValueError(f"{path}: timestamps are not strictly increasing")
    quaternion_norms = np.linalg.norm(poses[:, 4:8], axis=1)
    if np.any(np.abs(quaternion_norms - 1.0) > 1e-3):
        raise ValueError(f"{path}: invalid quaternion norm")
    return poses


def interpolate_positions(
    source: np.ndarray, query_times: np.ndarray, maximum_gap: float
) -> tuple[np.ndarray, np.ndarray]:
    right = np.searchsorted(source[:, 0], query_times, side="left")
    valid = (right > 0) & (right < len(source))
    left = np.clip(right - 1, 0, len(source) - 1)
    right = np.clip(right, 0, len(source) - 1)
    gaps = source[right, 0] - source[left, 0]
    valid &= gaps <= maximum_gap
    alpha = np.zeros_like(query_times)
    nonzero = gaps > 0.0
    alpha[nonzero] = (
        query_times[nonzero] - source[left[nonzero], 0]
    ) / gaps[nonzero]
    positions = (
        source[left, 1:4] * (1.0 - alpha[:, None])
        + source[right, 1:4] * alpha[:, None]
    )
    positions[~valid] = np.nan
    return positions, valid


def align_se3(estimate: np.ndarray, truth: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    estimate_center = estimate.mean(axis=0)
    truth_center = truth.mean(axis=0)
    covariance = (estimate - estimate_center).T @ (truth - truth_center)
    u, _, vt = np.linalg.svd(covariance)
    handedness = np.eye(3)
    handedness[2, 2] = np.linalg.det(vt.T @ u.T)
    rotation = vt.T @ handedness @ u.T
    translation = truth_center - rotation @ estimate_center
    return rotation, translation


def error_statistics(errors: np.ndarray) -> dict[str, float | int]:
    return {
        "count": int(len(errors)),
        "rmse_m": float(np.sqrt(np.mean(errors * errors))),
        "mean_m": float(np.mean(errors)),
        "median_m": float(np.median(errors)),
        "p95_m": float(np.percentile(errors, 95.0)),
        "max_m": float(np.max(errors)),
    }


def rpe_at_distance(
    estimate: np.ndarray, truth: np.ndarray, distance: float
) -> dict[str, float | int]:
    cumulative = np.concatenate(
        ([0.0], np.cumsum(np.linalg.norm(np.diff(truth, axis=0), axis=1)))
    )
    end = np.searchsorted(cumulative, cumulative + distance, side="left")
    start = np.arange(len(truth))
    valid = end < len(truth)
    errors = np.linalg.norm(
        (estimate[end[valid]] - estimate[start[valid]])
        - (truth[end[valid]] - truth[start[valid]]),
        axis=1,
    )
    if not len(errors):
        raise ValueError(f"trajectory is too short for {distance:g} m RPE")
    return error_statistics(errors)


def evaluate(path: Path, truth: np.ndarray, args: argparse.Namespace) -> dict[str, object]:
    estimate = load_poses(path, allow_extra_fields=False)
    interpolated, valid = interpolate_positions(
        estimate, truth[:, 0], args.max_interpolation_gap
    )
    if np.count_nonzero(valid) < 3:
        raise ValueError(f"{path}: fewer than three associated poses")
    associated_truth = truth[valid, 1:4]
    associated_estimate = interpolated[valid]
    rotation, translation = align_se3(associated_estimate, associated_truth)
    aligned = (rotation @ associated_estimate.T).T + translation
    ate = error_statistics(np.linalg.norm(aligned - associated_truth, axis=1))
    rpe = {
        f"{distance:g}m": rpe_at_distance(aligned, associated_truth, distance)
        for distance in (1.0, 10.0, 50.0)
    }
    coverage = float(np.count_nonzero(valid) / len(truth))
    checks = {
        "coverage": coverage >= args.minimum_coverage,
        "ate_rmse": ate["rmse_m"] <= args.ate_rmse_threshold,
        "rpe_10m_rmse": rpe["10m"]["rmse_m"] <= args.rpe_10m_rmse_threshold,
    }
    return {
        "trajectory": str(path.resolve()),
        "estimated_poses": int(len(estimate)),
        "associated_ground_truth_poses": int(np.count_nonzero(valid)),
        "coverage": coverage,
        "associated_time_span_s": float(truth[valid, 0][-1] - truth[valid, 0][0]),
        "alignment": "fixed-scale SE(3) Horn/SVD; scale fixed to 1.0",
        "se3_rotation": rotation.tolist(),
        "se3_translation_m": translation.tolist(),
        "ate_translation": ate,
        "rpe_translation": rpe,
        "checks": checks,
        "passed": all(checks.values()),
    }


def markdown(report: dict[str, object]) -> str:
    lines = [
        "# SANY trajectory evaluation",
        "",
        "Fixed-scale SE(3) alignment is used; scale estimation is prohibited.",
        "",
        "| Run | Coverage | ATE RMSE (m) | RPE 10 m RMSE (m) | Result |",
        "|---|---:|---:|---:|---|",
    ]
    for label, run in report["runs"].items():
        lines.append(
            f"| {label} | {run['coverage']:.3%} | "
            f"{run['ate_translation']['rmse_m']:.4f} | "
            f"{run['rpe_translation']['10m']['rmse_m']:.4f} | "
            f"{'PASS' if run['passed'] else 'FAIL'} |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    if not all(
        math.isfinite(value) and value > 0.0
        for value in (
            args.max_interpolation_gap,
            args.ate_rmse_threshold,
            args.rpe_10m_rmse_threshold,
            args.minimum_coverage,
        )
    ):
        raise ValueError("evaluation thresholds must be finite and positive")
    truth = load_poses(args.ground_truth, allow_extra_fields=True)
    report: dict[str, object] = {
        "ground_truth": str(args.ground_truth.resolve()),
        "ground_truth_poses": int(len(truth)),
        "protocol": {
            "alignment": "fixed-scale SE(3); Sim(3) prohibited",
            "maximum_interpolation_gap_s": args.max_interpolation_gap,
            "minimum_coverage": args.minimum_coverage,
            "ate_rmse_threshold_m": args.ate_rmse_threshold,
            "rpe_10m_rmse_threshold_m": args.rpe_10m_rmse_threshold,
        },
        "runs": {},
    }
    for assignment in args.run:
        label, separator, raw_path = assignment.partition("=")
        if not separator or not label or not raw_path:
            raise ValueError(f"expected LABEL=TUM, got {assignment!r}")
        if label in report["runs"]:
            raise ValueError(f"duplicate run label: {label}")
        report["runs"][label] = evaluate(Path(raw_path), truth, args)
    report["passed"] = all(run["passed"] for run in report["runs"].values())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    args.output.with_suffix(".md").write_text(markdown(report), encoding="utf-8")
    print(args.output.resolve())
    print(f"passed={report['passed']}")


if __name__ == "__main__":
    main()
