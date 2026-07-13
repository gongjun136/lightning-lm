#!/usr/bin/env python3
"""Extract structured Lightning-LM timing statistics from the offline log."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import OrderedDict
from pathlib import Path


NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
TIMER_LINE = re.compile(
    rf"> \[(?P<stage>[^\]]+)\]\s+average time usage:\s*"
    rf"(?P<average>{NUMBER})\s*ms,\s*med:\s*(?P<median>{NUMBER})\s+"
    rf"95%:\s*(?P<p95>{NUMBER}),\s*called times\s*:\s*(?P<samples>\d+)"
)
TIMER_SAMPLE_WINDOW_LIMIT = 2000


def positive_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a finite number greater than zero")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < 0:
        raise argparse.ArgumentTypeError("value must be a finite number zero or greater")
    return parsed


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return parsed


def boolean(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise argparse.ArgumentTypeError("value must be true/false or 1/0")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--wall-time-s", required=True, type=positive_float)
    parser.add_argument("--sensor-duration-s", required=True, type=positive_float)
    parser.add_argument("--trajectory-frames", required=True, type=nonnegative_int)
    parser.add_argument("--completion", required=True)
    parser.add_argument("--max-lidar-frames", required=True, type=nonnegative_int)
    parser.add_argument("--algorithm-rc", required=True, type=nonnegative_int)
    parser.add_argument("--watchdog-status", required=True)
    parser.add_argument("--playback-rate", required=True, type=nonnegative_float)
    parser.add_argument("--wait-ui", required=True, type=boolean)
    return parser.parse_args()


def parse_timer_rows(log_path: Path) -> list[dict[str, float | int | str]]:
    stages: OrderedDict[str, dict[str, float | int | str]] = OrderedDict()
    with log_path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = TIMER_LINE.search(line)
            if not match:
                continue
            stage = match.group("stage").strip()
            retained_samples = int(match.group("samples"))
            average_ms = float(match.group("average"))
            stages[stage] = {
                "stage": stage,
                "average_ms": average_ms,
                "median_ms": float(match.group("median")),
                "p95_ms": float(match.group("p95")),
                "retained_samples": retained_samples,
                "retained_window_estimated_total_ms": round(average_ms * retained_samples, 6),
                "sample_window_limit_reached": retained_samples >= TIMER_SAMPLE_WINDOW_LIMIT,
            }
    return list(stages.values())


def write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    fields = [
        "stage",
        "average_ms",
        "median_ms",
        "p95_ms",
        "retained_samples",
        "retained_window_estimated_total_ms",
        "sample_window_limit_reached",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, args: argparse.Namespace, rows: list[dict[str, float | int | str]]) -> None:
    ratio_suppressed_reasons = []
    if args.max_lidar_frames != 0:
        ratio_suppressed_reasons.append("frame_limited_run")
    elif args.completion != "reached_final_lidar":
        ratio_suppressed_reasons.append("input_not_fully_processed")
    if args.algorithm_rc != 0:
        ratio_suppressed_reasons.append("algorithm_failed")
    if args.watchdog_status != "completed":
        ratio_suppressed_reasons.append("watchdog_not_completed")
    if args.wait_ui:
        ratio_suppressed_reasons.append("ui_wait_included_in_wall_time")
    full_bag_performance_valid = not ratio_suppressed_reasons
    wall_throughput_valid = not args.wait_ui
    end_to_end = {
        "wall_time_s": args.wall_time_s,
        "input_sensor_duration_s": args.sensor_duration_s,
        "trajectory_frames": args.trajectory_frames,
        "trajectory_frames_per_wall_s": (
            args.trajectory_frames / args.wall_time_s
            if args.trajectory_frames and wall_throughput_valid
            else None
        ),
        "wall_ms_per_trajectory_frame": (
            args.wall_time_s * 1000.0 / args.trajectory_frames
            if args.trajectory_frames and wall_throughput_valid
            else None
        ),
        "realtime_factor": args.wall_time_s / args.sensor_duration_s if full_bag_performance_valid else None,
        "processing_speed_x": args.sensor_duration_s / args.wall_time_s if full_bag_performance_valid else None,
    }
    payload = {
        "schema_version": 1,
        "status": "ok" if rows else "no_timer_records",
        "source_log": str(args.log.resolve()),
        "completion": args.completion,
        "max_lidar_frames": args.max_lidar_frames,
        "algorithm_rc": args.algorithm_rc,
        "watchdog_status": args.watchdog_status,
        "playback_rate": args.playback_rate,
        "wait_ui": args.wait_ui,
        "performance_ratio_status": "available" if full_bag_performance_valid else "suppressed",
        "performance_ratio_suppressed_reasons": ratio_suppressed_reasons,
        "stage_count": len(rows),
        "timer_sample_window_limit": TIMER_SAMPLE_WINDOW_LIMIT,
        "timer_semantics": (
            "Each stage keeps at most the latest 2000 measurements. Stage timers may overlap or be nested; "
            "retained_window_estimated_total_ms values must not be summed as end-to-end wall time."
        ),
        "end_to_end": end_to_end,
        "stages": rows,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def main() -> int:
    args = parse_args()
    rows = parse_timer_rows(args.log)
    write_csv(args.csv, rows)
    write_summary(args.summary, args, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
