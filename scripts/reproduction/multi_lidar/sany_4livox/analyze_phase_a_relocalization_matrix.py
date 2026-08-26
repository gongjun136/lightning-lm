#!/usr/bin/env python3
"""Summarize the SANY Phase-A multi-offset cold-start matrix."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


EXPECTED_DATASETS = ("data1", "data2", "data3", "data4", "data5")
EXPECTED_OFFSETS_S = (0, 3, 5, 7, 10, 15, 20, 25, 30, 35)
SUCCESS_RATE_MINIMUM = 1.0


def unit_interval_float(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"expected a float, got {value!r}") from error
    if not math.isfinite(parsed) or not 0.0 <= parsed <= 1.0:
        raise argparse.ArgumentTypeError("expected a finite float in [0, 1]")
    return parsed


def percentile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def yaw_deg(row: dict[str, str]) -> float:
    qx = float(row["lidar_loc_qx"])
    qy = float(row["lidar_loc_qy"])
    qz = float(row["lidar_loc_qz"])
    qw = float(row["lidar_loc_qw"])
    return math.degrees(math.atan2(2.0 * (qw * qz + qx * qy), 1.0 - 2.0 * (qy * qy + qz * qz)))


def angle_difference_deg(lhs: float, rhs: float) -> float:
    return abs((lhs - rhs + 180.0) % 360.0 - 180.0)


def read_trial(path: Path, dataset: str, offset_s: int) -> dict[str, object]:
    stats = path / "results" / "localization_stats.csv"
    if not stats.is_file():
        return {"dataset": dataset, "offset_s": offset_s, "success": False, "reason": "missing_stats"}
    with stats.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return {"dataset": dataset, "offset_s": offset_s, "success": False, "reason": "empty_stats"}
    accepted = next((row for row in rows if row["relocalization_accepted"] == "1"), None)
    if accepted is None:
        attempts = [row for row in rows if row["relocalization_attempted"] == "1"]
        candidates = sum(row["relocalization_candidate_found"] == "1" for row in rows)
        return {
            "dataset": dataset,
            "offset_s": offset_s,
            "success": False,
            "reason": "no_accepted_pose",
            "frames": len(rows),
            "attempts": len(attempts),
            "candidate_frames": candidates,
            "observed_processing_total_ms": sum(float(row["processing_ms"]) for row in rows),
            "observed_search_total_ms": sum(
                float(row["relocalization_search_time_ms"]) for row in attempts
            ),
        }
    accepted_index = rows.index(accepted)
    startup_rows = rows[: accepted_index + 1]
    attempts = [row for row in startup_rows if row["relocalization_attempted"] == "1"]
    first_attempt = attempts[0]
    return {
        "dataset": dataset,
        "offset_s": offset_s,
        "success": True,
        "reason": "accepted",
        "frames": len(rows),
        "accepted_frame": int(accepted["frame_index"]),
        "delay_s": float(accepted["timestamp"]) - float(rows[0]["timestamp"]),
        "attempts": len(attempts),
        "startup_processing_total_ms": sum(
            float(row["processing_ms"]) for row in startup_rows
        ),
        "relocalization_search_total_ms": sum(
            float(row["relocalization_search_time_ms"]) for row in attempts
        ),
        "first_attempt_processing_ms": float(first_attempt["processing_ms"]),
        "first_attempt_search_ms": float(
            first_attempt["relocalization_search_time_ms"]
        ),
        "score": float(accepted["relocalization_score"]),
        "overlap": float(accepted["map_overlap_ratio"]),
        "confirmations": int(accepted["relocalization_confirmation_count"]),
        "x": float(accepted["lidar_loc_x"]),
        "y": float(accepted["lidar_loc_y"]),
        "z": float(accepted["lidar_loc_z"]),
        "yaw_deg": yaw_deg(accepted),
    }


def add_reference_error(
    trial: dict[str, object], reference_path: Path | None
) -> None:
    if not trial["success"] or reference_path is None:
        return
    reference = read_trial(
        reference_path, str(trial["dataset"]), int(trial["offset_s"])
    )
    if not reference["success"]:
        return
    trial["reference_position_error_m"] = math.dist(
        [float(trial[key]) for key in ("x", "y", "z")],
        [float(reference[key]) for key in ("x", "y", "z")],
    )
    trial["reference_yaw_error_deg"] = angle_difference_deg(
        float(trial["yaw_deg"]), float(reference["yaw_deg"])
    )


def summarize(
    trials: list[dict[str, object]], map_overlap_minimum: float
) -> dict[str, object]:
    accepted = [trial for trial in trials if trial["success"]]
    delays = [float(trial["delay_s"]) for trial in accepted]
    attempts = [int(trial["attempts"]) for trial in accepted]
    startup_processing = [
        float(trial["startup_processing_total_ms"]) for trial in accepted
    ]
    relocalization_search = [
        float(trial["relocalization_search_total_ms"]) for trial in accepted
    ]
    first_attempt_processing = [
        float(trial["first_attempt_processing_ms"]) for trial in accepted
    ]
    first_attempt_search = [
        float(trial["first_attempt_search_ms"]) for trial in accepted
    ]
    reference_trials = [
        trial for trial in accepted if "reference_position_error_m" in trial
    ]
    summary: dict[str, object] = {
        "trials": len(trials),
        "accepted": len(accepted),
        "success_rate": len(accepted) / len(trials),
        "within_5s": sum(delay <= 5.0 for delay in delays),
        "within_5s_rate": sum(delay <= 5.0 for delay in delays) / len(trials),
        "delay_p95_s": percentile(delays, 0.95),
        "delay_max_s": max(delays) if delays else None,
        "attempts_p95": percentile(attempts, 0.95),
        "attempts_max": max(attempts) if attempts else None,
        "startup_processing_total_p50_ms": percentile(startup_processing, 0.50),
        "startup_processing_total_p95_ms": percentile(startup_processing, 0.95),
        "startup_processing_total_max_ms": (
            max(startup_processing) if startup_processing else None
        ),
        "relocalization_search_total_p50_ms": percentile(
            relocalization_search, 0.50
        ),
        "relocalization_search_total_p95_ms": percentile(
            relocalization_search, 0.95
        ),
        "relocalization_search_total_max_ms": (
            max(relocalization_search) if relocalization_search else None
        ),
        "first_attempt_processing_p50_ms": percentile(
            first_attempt_processing, 0.50
        ),
        "first_attempt_processing_p95_ms": percentile(
            first_attempt_processing, 0.95
        ),
        "first_attempt_search_p50_ms": percentile(first_attempt_search, 0.50),
        "first_attempt_search_p95_ms": percentile(first_attempt_search, 0.95),
        "minimum_overlap": min(
            (float(trial["overlap"]) for trial in accepted), default=None
        ),
        "minimum_confirmations": min(
            (int(trial["confirmations"]) for trial in accepted), default=None
        ),
        "unsafe_overlap_accepts": sum(
            float(trial["overlap"]) < map_overlap_minimum for trial in accepted
        ),
        "reference_comparisons": len(reference_trials),
        "reference_position_error_p95_m": percentile(
            [float(trial["reference_position_error_m"]) for trial in reference_trials],
            0.95,
        ),
        "reference_position_error_max_m": max(
            (float(trial["reference_position_error_m"]) for trial in reference_trials),
            default=None,
        ),
        "reference_yaw_error_max_deg": max(
            (float(trial["reference_yaw_error_deg"]) for trial in reference_trials),
            default=None,
        ),
    }
    if accepted:
        center = [statistics.median(float(trial[key]) for trial in accepted) for key in ("x", "y", "z")]
        center_yaw = statistics.median(float(trial["yaw_deg"]) for trial in accepted)
        summary["accepted_pose_median_xyz"] = center
        summary["accepted_pose_max_spread_m"] = max(
            math.dist(center, [float(trial[key]) for key in ("x", "y", "z")]) for trial in accepted
        )
        summary["accepted_yaw_max_spread_deg"] = max(
            angle_difference_deg(float(trial["yaw_deg"]), center_yaw) for trial in accepted
        )
    else:
        summary["accepted_pose_median_xyz"] = None
        summary["accepted_pose_max_spread_m"] = None
        summary["accepted_yaw_max_spread_deg"] = None
    reference_safe = all(
        float(trial["reference_position_error_m"]) <= 1.0
        and float(trial["reference_yaw_error_deg"]) <= 5.0
        for trial in reference_trials
    )
    summary["passed"] = bool(
        summary["success_rate"] >= SUCCESS_RATE_MINIMUM
        and summary["delay_p95_s"] is not None
        and float(summary["delay_p95_s"]) <= 5.0
        and summary["minimum_confirmations"] is not None
        and int(summary["minimum_confirmations"]) >= 2
        and summary["minimum_overlap"] is not None
        and float(summary["minimum_overlap"]) >= map_overlap_minimum
        and reference_safe
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-root", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument(
        "--map-overlap-minimum",
        type=unit_interval_float,
        default=0.95,
        help="Minimum accepted map overlap ratio (default: 0.95)",
    )
    parser.add_argument(
        "--reference-runs-root",
        type=Path,
        help="Optional same-offset reference runs used for accepted-pose checks",
    )
    args = parser.parse_args()

    trials: list[dict[str, object]] = []
    for dataset in EXPECTED_DATASETS:
        for offset_s in EXPECTED_OFFSETS_S:
            trial_dir = args.runs_root / dataset / f"offset_{offset_s:02d}s"
            trial = read_trial(trial_dir, dataset, offset_s)
            reference_path = None
            if args.reference_runs_root is not None:
                reference_path = (
                    args.reference_runs_root / dataset / f"offset_{offset_s:02d}s"
                )
            add_reference_error(trial, reference_path)
            trials.append(trial)

    datasets = {
        dataset: summarize(
            [trial for trial in trials if trial["dataset"] == dataset],
            args.map_overlap_minimum,
        )
        for dataset in EXPECTED_DATASETS
    }
    overall = summarize(trials, args.map_overlap_minimum)
    overall["datasets_passed"] = sum(bool(summary["passed"]) for summary in datasets.values())
    overall["all_datasets_passed"] = all(bool(summary["passed"]) for summary in datasets.values())
    overall["passed"] = bool(overall["passed"] and overall["all_datasets_passed"])
    report = {
        "schema_version": 1,
        "runs_root": str(args.runs_root.resolve()),
        "criteria": {
            "success_rate_minimum": SUCCESS_RATE_MINIMUM,
            "delay_p95_maximum_s": 5.0,
            "confirmation_minimum": 2,
            "map_overlap_minimum": args.map_overlap_minimum,
            "reference_position_error_maximum_m": 1.0,
            "reference_yaw_error_maximum_deg": 5.0,
            "reference_scope": "only trials accepted by the optional reference runs",
            "accepted_pose_spread": "diagnostic_only_to_allow_dynamic_initialization",
        },
        "datasets": datasets,
        "overall": overall,
        "trials": trials,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    fields = sorted({key for trial in trials for key in trial})
    with args.output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(trials)
    print(args.output_json.resolve())
    print(f"datasets={len(datasets)} trials={len(trials)} passed={report['overall']['passed']}")


if __name__ == "__main__":
    main()
