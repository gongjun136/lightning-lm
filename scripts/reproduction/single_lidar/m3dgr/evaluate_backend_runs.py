#!/usr/bin/env python3
"""Evaluate M3DGR backend runs with fixed-scale SE(3) and RTK loop truth."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument(
        "--run", action="append", required=True, metavar="LABEL=TUM",
        help="repeat for every trajectory being compared",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-interpolation-gap", type=float, default=0.2)
    parser.add_argument("--rpe-distance", type=float, action="append", default=[])
    parser.add_argument("--loop-min-time", type=float, default=30.0)
    parser.add_argument("--loop-max-distance", type=float, default=5.0)
    return parser.parse_args()


def split_assignment(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"expected LABEL=PATH, got {value!r}")
    label, path = value.split("=", 1)
    if not label or not path:
        raise ValueError(f"expected LABEL=PATH, got {value!r}")
    return label, Path(path)


def load_tum(path: Path, *, ground_truth: bool) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8-sig") as stream:
        for line_number, raw in enumerate(stream, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 8:
                raise ValueError(f"{path}:{line_number}: expected 8 TUM fields")
            values = [float(value) for value in fields]
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{path}:{line_number}: non-finite pose")
            rows.append(values)
    if not rows:
        raise ValueError(f"{path}: no poses")
    poses = np.asarray(rows, dtype=np.float64)
    norms = np.linalg.norm(poses[:, 4:8], axis=1)
    if np.any(np.abs(norms - 1.0) > 1e-3):
        raise ValueError(f"{path}: invalid quaternion norm")
    poses[:, 4:8] /= norms[:, None]

    if ground_truth:
        poses = poses[np.argsort(poses[:, 0], kind="stable")]
        retained: list[np.ndarray] = []
        begin = 0
        while begin < len(poses):
            end = begin + 1
            while end < len(poses) and poses[end, 0] == poses[begin, 0]:
                end += 1
            group = poses[begin:end]
            if np.max(np.linalg.norm(group[:, 1:4] - group[0, 1:4], axis=1)) > 1e-8:
                raise ValueError(f"{path}: duplicate timestamp has inconsistent positions")
            retained.append(group[0])
            begin = end
        poses = np.asarray(retained)
    if np.any(np.diff(poses[:, 0]) <= 0.0):
        raise ValueError(f"{path}: timestamps are not strictly increasing")
    return poses


def interpolate_positions(
    source: np.ndarray, query_times: np.ndarray, max_gap: float
) -> tuple[np.ndarray, np.ndarray]:
    right = np.searchsorted(source[:, 0], query_times, side="left")
    valid = (right > 0) & (right < len(source))
    left = np.clip(right - 1, 0, len(source) - 1)
    right = np.clip(right, 0, len(source) - 1)
    gaps = source[right, 0] - source[left, 0]
    valid &= gaps <= max_gap
    alpha = np.zeros_like(query_times, dtype=np.float64)
    nonzero = gaps > 0.0
    alpha[nonzero] = (
        (query_times[nonzero] - source[left[nonzero], 0]) / gaps[nonzero]
    )
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
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(vt.T @ u.T)
    rotation = vt.T @ correction @ u.T
    translation = truth_center - rotation @ estimate_center
    return rotation, translation


def statistics(errors: np.ndarray) -> dict[str, float | int | None]:
    if len(errors) == 0:
        return {"count": 0, "rmse_m": None, "mean_m": None, "median_m": None,
                "std_m": None, "max_m": None}
    return {
        "count": int(len(errors)),
        "rmse_m": float(np.sqrt(np.mean(errors * errors))),
        "mean_m": float(np.mean(errors)),
        "median_m": float(np.median(errors)),
        "std_m": float(np.std(errors)),
        "max_m": float(np.max(errors)),
    }


def rpe_by_distance(
    estimate: np.ndarray, truth: np.ndarray, distances: list[float]
) -> dict[str, dict[str, float | int | None]]:
    cumulative = np.concatenate(
        ([0.0], np.cumsum(np.linalg.norm(np.diff(truth, axis=0), axis=1)))
    )
    result: dict[str, dict[str, float | int | None]] = {}
    for distance in distances:
        end = np.searchsorted(cumulative, cumulative + distance, side="left")
        start = np.arange(len(truth))
        valid = end < len(truth)
        delta_estimate = estimate[end[valid]] - estimate[start[valid]]
        delta_truth = truth[end[valid]] - truth[start[valid]]
        errors = np.linalg.norm(delta_estimate - delta_truth, axis=1)
        result[f"{distance:g}m"] = statistics(errors)
    return result


def find_run_root(trajectory: Path) -> Path | None:
    for parent in trajectory.resolve().parents:
        if (parent / "run_metadata.txt").is_file():
            return parent
    return None


def read_optional_yaml(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    with path.open("r", encoding="utf-8") as stream:
        value = yaml.safe_load(stream)
    return value if isinstance(value, dict) else None


def read_optional_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    return value if isinstance(value, dict) else None


def read_run_metadata(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    values: dict[str, Any] = {}
    with path.open("r", encoding="utf-8-sig") as stream:
        for line in stream:
            key, separator, value = line.strip().partition("=")
            if not separator or not key:
                continue
            try:
                values[key] = float(value)
            except ValueError:
                values[key] = value
    return values


def load_loop_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def loop_metrics(
    rows: list[dict[str, str]], truth: np.ndarray, min_time: float, max_distance: float
) -> dict[str, Any] | None:
    if not rows:
        return None
    current_times = np.asarray([float(row["current_timestamp"]) for row in rows])
    current_positions, current_valid = interpolate_positions(truth, current_times, 0.5)
    opportunities = np.zeros(len(rows), dtype=bool)
    opportunity_pairs: list[dict[str, Any]] = []
    for current in range(len(rows)):
        if not current_valid[current]:
            continue
        prior = np.flatnonzero(
            current_valid[:current] & ((current_times[current] - current_times[:current]) >= min_time)
        )
        if len(prior):
            distances = np.linalg.norm(
                current_positions[prior] - current_positions[current], axis=1
            )
            nearest = int(np.argmin(distances))
            if distances[nearest] <= max_distance:
                history = int(prior[nearest])
                opportunities[current] = True
                opportunity_pairs.append(
                    {
                        "current_descriptor": int(rows[current]["current_descriptor"]),
                        "nearest_history_descriptor": int(rows[history]["current_descriptor"]),
                        "time_separation_s": float(current_times[current] - current_times[history]),
                        "rtk_distance_m": float(distances[nearest]),
                    }
                )

    candidate = np.asarray([int(row.get("candidate", "0")) != 0 for row in rows])
    accepted = np.asarray([int(row.get("accepted", "0")) != 0 for row in rows])
    candidate_true = np.zeros(len(rows), dtype=bool)
    accepted_true = np.zeros(len(rows), dtype=bool)
    candidate_details: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        if not candidate[index] or not current_valid[index]:
            continue
        history_time = float(row.get("history_timestamp", "0"))
        history_position, valid = interpolate_positions(
            truth, np.asarray([history_time]), 0.5
        )
        if valid[0]:
            time_separation = current_times[index] - history_time
            rtk_distance = float(
                np.linalg.norm(current_positions[index] - history_position[0])
            )
            candidate_true[index] = (
                time_separation >= min_time and rtk_distance <= max_distance
            )
            accepted_true[index] = accepted[index] and candidate_true[index]
            candidate_details.append(
                {
                    "current_descriptor": int(row["current_descriptor"]),
                    "history_descriptor": int(row["history_descriptor"]),
                    "time_separation_s": time_separation,
                    "rtk_distance_m": rtk_distance,
                    "score": float(row.get("score", "0")),
                    "accepted": bool(accepted[index]),
                    "rtk_true": bool(candidate_true[index]),
                    "reason": row.get("reason", ""),
                }
            )

    descriptor_interval = float(np.median(np.diff(current_times))) if len(rows) > 1 else 1.0
    episode_gap = max(2.5, 3.0 * descriptor_interval)
    episode_indices: list[list[int]] = []
    for index in np.flatnonzero(opportunities):
        if not episode_indices or current_times[index] - current_times[episode_indices[-1][-1]] > episode_gap:
            episode_indices.append([int(index)])
        else:
            episode_indices[-1].append(int(index))
    covered_episodes = sum(bool(np.any(accepted_true[episode])) for episode in episode_indices)
    accepted_count = int(np.count_nonzero(accepted))
    true_count = int(np.count_nonzero(accepted_true))
    candidate_count = int(np.count_nonzero(candidate))
    true_candidate_count = int(np.count_nonzero(candidate_true))
    opportunity_count = int(np.count_nonzero(opportunities))
    return {
        "truth_definition": {
            "minimum_time_separation_s": min_time,
            "maximum_rtk_distance_m": max_distance,
            "unit": "descriptor-query loop episode",
        },
        "retrieved_candidates": candidate_count,
        "true_retrieved_candidates": true_candidate_count,
        "candidate_precision": (
            true_candidate_count / candidate_count if candidate_count else None
        ),
        "candidates": candidate_details,
        "accepted_detections": accepted_count,
        "true_accepted_detections": true_count,
        "false_accepted_detections": accepted_count - true_count,
        "precision": true_count / accepted_count if accepted_count else None,
        "opportunity_queries": opportunity_count,
        "opportunities": opportunity_pairs,
        "true_opportunity_queries_detected": int(np.count_nonzero(accepted_true & opportunities)),
        "query_recall": (
            float(np.count_nonzero(accepted_true & opportunities)) / opportunity_count
            if opportunity_count else None
        ),
        "truth_episodes": len(episode_indices),
        "detected_truth_episodes": covered_episodes,
        "recall": covered_episodes / len(episode_indices) if episode_indices else None,
    }


def evaluate_run(
    trajectory_path: Path,
    truth: np.ndarray,
    interpolation_gap: float,
    rpe_distances: list[float],
    loop_min_time: float,
    loop_max_distance: float,
) -> dict[str, Any]:
    estimate = load_tum(trajectory_path, ground_truth=False)
    estimate_positions, valid = interpolate_positions(
        estimate, truth[:, 0], interpolation_gap
    )
    if np.count_nonzero(valid) < 3:
        raise ValueError(f"{trajectory_path}: fewer than three GT poses can be associated")
    associated_truth = truth[valid, 1:4]
    associated_estimate = estimate_positions[valid]
    rotation, translation = align_se3(associated_estimate, associated_truth)
    aligned_estimate = (rotation @ associated_estimate.T).T + translation
    ate = np.linalg.norm(aligned_estimate - associated_truth, axis=1)

    result: dict[str, Any] = {
        "trajectory": str(trajectory_path.resolve()),
        "alignment": "fixed-scale SE(3) Horn/SVD; scale fixed to 1.0",
        "estimated_poses": int(len(estimate)),
        "associated_gt_poses": int(np.count_nonzero(valid)),
        "gt_pose_coverage": float(np.count_nonzero(valid) / len(truth)),
        "associated_time_span_s": float(truth[valid, 0][-1] - truth[valid, 0][0]),
        "se3_rotation": rotation.tolist(),
        "se3_translation_m": translation.tolist(),
        "ate_translation": statistics(ate),
        "rpe_translation": rpe_by_distance(aligned_estimate, associated_truth, rpe_distances),
        "orientation_metrics": None,
        "orientation_note": "M3DGR outdoor RTK files contain identity placeholder quaternions",
    }

    run_root = find_run_root(trajectory_path)
    if run_root is not None:
        diagnostics = run_root / "data" / "new_map" / "backend_diagnostics"
        result["run_root"] = str(run_root)
        result["run_metadata"] = read_run_metadata(run_root / "run_metadata.txt")
        result["resource"] = read_optional_json(run_root / "resource_summary.json")
        result["backend"] = read_optional_yaml(diagnostics / "backend_summary.yaml")
        result["loop"] = loop_metrics(
            load_loop_rows(diagnostics / "btc_loop_candidates.csv"),
            truth,
            loop_min_time,
            loop_max_distance,
        )
    return result


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# M3DGR backend evaluation",
        "",
        "Fixed-scale SE(3) alignment is used; no scale is estimated. RTK is used only here as ground truth.",
        "",
        "| Run | ATE RMSE (m) | ATE median (m) | RPE 10 m RMSE (m) | Loop P | Loop R | Wall (s) | Peak RSS (MB) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for label, run in report["runs"].items():
        ate = run["ate_translation"]
        rpe = run["rpe_translation"].get("10m", {})
        loop = run.get("loop") or {}
        resource = run.get("resource") or {}
        metadata = run.get("run_metadata") or {}
        def value(number: Any, digits: int = 4) -> str:
            return "N/A" if number is None else f"{float(number):.{digits}f}"
        lines.append(
            f"| {label} | {value(ate['rmse_m'])} | {value(ate['median_m'])} | "
            f"{value(rpe.get('rmse_m'))} | {value(loop.get('precision'), 3)} | "
            f"{value(loop.get('recall'), 3)} | "
            f"{value(metadata.get('wall_time_s', resource.get('duration_s')), 2)} | "
            f"{value(resource.get('peak_rss_mb'), 1)} |"
        )
    lines.extend(
        [
            "",
            "Loop truth: RTK separation <= configured distance and temporal separation >= configured time. "
            "Recall is computed over contiguous descriptor-query loop episodes; query-level recall is also retained in JSON.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    truth = load_tum(args.ground_truth, ground_truth=True)
    distances = args.rpe_distance or [1.0, 10.0, 100.0]
    report: dict[str, Any] = {
        "ground_truth": str(args.ground_truth.resolve()),
        "ground_truth_poses": int(len(truth)),
        "protocol": {
            "alignment": "fixed-scale SE(3); Sim(3) prohibited",
            "maximum_interpolation_gap_s": args.max_interpolation_gap,
            "rpe_distances_m": distances,
            "loop_minimum_time_separation_s": args.loop_min_time,
            "loop_maximum_rtk_distance_m": args.loop_max_distance,
            "rtk_optimizer_input": False,
        },
        "runs": {},
    }
    for assignment in args.run:
        label, path = split_assignment(assignment)
        if label in report["runs"]:
            raise ValueError(f"duplicate run label: {label}")
        report["runs"][label] = evaluate_run(
            path,
            truth,
            args.max_interpolation_gap,
            distances,
            args.loop_min_time,
            args.loop_max_distance,
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    markdown_path = args.output.with_suffix(".md")
    markdown_path.write_text(markdown_report(report), encoding="utf-8")
    print(args.output.resolve())
    print(markdown_path.resolve())


if __name__ == "__main__":
    main()
