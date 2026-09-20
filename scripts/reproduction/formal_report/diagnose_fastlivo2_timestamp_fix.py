#!/usr/bin/env python3
"""Validate the FAST-LIVO2 sensor-time timestamp fix against archived runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


SEQUENCES = ("Grass02", "Outdoor04", "Z-Rough-Road01", "Dark01")
METHOD = "fastlivo2_lio"


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    formal = repo / "runs" / "formal_report_20260717"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--old-analysis",
        type=Path,
        default=formal / "analysis" / "_superseded_m3dgr_frontend_ros_now_20260720",
    )
    parser.add_argument(
        "--new-analysis", type=Path, default=formal / "analysis" / "m3dgr_frontend",
    )
    parser.add_argument(
        "--old-runs",
        type=Path,
        default=formal / "m3dgr_frontend_v2" / "_superseded_fastlivo_ros_now_20260720",
    )
    parser.add_argument(
        "--new-runs", type=Path, default=formal / "m3dgr_frontend_v2",
    )
    parser.add_argument("--spike-excess-m", type=float, default=0.5)
    parser.add_argument("--gap-threshold-s", type=float, default=0.2)
    parser.add_argument("--post-gap-window-s", type=float, default=2.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=formal / "analysis" / "m3dgr_frontend" / "fastlivo2_timestamp_fix_validation.json",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def parse_metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    return values


def summary_lookup(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    return {(row["sequence"], row["method"]): row for row in read_csv(path)}


def old_spike_diagnostics(
    rows: list[dict[str, str]],
    sequence: str,
    old_tum: Path,
    spike_excess_m: float,
    gap_threshold_s: float,
    post_gap_window_s: float,
) -> dict[str, float | int]:
    repeat_rows = [
        row for row in rows if row["sequence"] == sequence and int(row["repeat"]) == 1
    ]
    errors: dict[tuple[str, float], float] = {
        (row["method"], float(row["timestamp_s"])): float(row["translation_error_m"])
        for row in repeat_rows
    }
    fast_times = sorted(time for method, time in errors if method == METHOD)
    excess = np.asarray([
        errors[(METHOD, time)]
        - np.median([value for (method, stamp), value in errors.items() if stamp == time and method != METHOD])
        for time in fast_times
    ])
    times = np.asarray(fast_times, dtype=np.float64)
    trajectory = np.loadtxt(old_tum, dtype=np.float64)
    gaps = np.diff(trajectory[:, 0])
    gap_right_indices = np.flatnonzero(gaps > gap_threshold_s) + 1
    gap_right_times = trajectory[gap_right_indices, 0]
    latest_gap = np.searchsorted(gap_right_times, times, side="right") - 1
    since_gap = np.full(len(times), np.inf)
    valid = latest_gap >= 0
    since_gap[valid] = times[valid] - gap_right_times[latest_gap[valid]]
    spikes = excess > spike_excess_m
    inside_window = spikes & (since_gap >= 0.0) & (since_gap <= post_gap_window_s)
    spike_count = int(np.count_nonzero(spikes))
    captured = int(np.count_nonzero(inside_window))
    return {
        "representative_repeat": 1,
        "old_output_gap_count": int(np.count_nonzero(gaps > gap_threshold_s)),
        "old_max_output_gap_s": float(np.max(gaps)),
        "fastlivo_excess_spike_count": spike_count,
        "spikes_within_post_gap_window": captured,
        "spike_capture_ratio": captured / spike_count if spike_count else 1.0,
        "maximum_excess_error_m": float(np.max(excess)),
    }


def excess_spike_summary(
    rows: list[dict[str, str]], sequence: str, spike_excess_m: float,
) -> dict[str, float | int]:
    repeat_rows = [
        row for row in rows if row["sequence"] == sequence and int(row["repeat"]) == 1
    ]
    errors: dict[tuple[str, float], float] = {
        (row["method"], float(row["timestamp_s"])): float(row["translation_error_m"])
        for row in repeat_rows
    }
    fast_times = sorted(time for method, time in errors if method == METHOD)
    excess = np.asarray([
        errors[(METHOD, time)]
        - np.median([value for (method, stamp), value in errors.items() if stamp == time and method != METHOD])
        for time in fast_times
    ])
    return {
        "spike_count": int(np.count_nonzero(excess > spike_excess_m)),
        "maximum_excess_error_m": float(np.max(excess)),
    }


def main() -> int:
    args = parse_args()
    old_aligned_path = args.old_analysis / "aligned_samples.csv"
    old_rows = read_csv(old_aligned_path)
    new_aligned_path = args.new_analysis / "aligned_samples.csv"
    new_rows = read_csv(new_aligned_path)
    old_summary = summary_lookup(args.old_analysis / "summary_metrics.csv")
    new_summary = summary_lookup(args.new_analysis / "summary_metrics.csv")
    repo = Path(__file__).resolve().parents[3]
    source = (
        repo.parents[1]
        / "WSL_Ubuntu_20.04"
        / "ros1_ws"
        / "ws_fastlivo2"
        / "src"
        / "FAST-LIVO2"
        / "src"
        / "LIVMapper.cpp"
    )

    sequences: dict[str, object] = {}
    new_contract_rows: list[dict[str, object]] = []
    for sequence in SEQUENCES:
        old_tum = args.old_runs / sequence / METHOD / "repeat_01" / "trajectory_mid360.tum"
        diagnostics = old_spike_diagnostics(
            old_rows,
            sequence,
            old_tum,
            args.spike_excess_m,
            args.gap_threshold_s,
            args.post_gap_window_s,
        )
        old_metrics = old_summary[(sequence, METHOD)]
        new_metrics = new_summary[(sequence, METHOD)]
        diagnostics["old_ate_rmse_m"] = float(old_metrics["ate_rmse_m_mean"])
        diagnostics["new_ate_rmse_m"] = float(new_metrics["ate_rmse_m_mean"])
        diagnostics["old_ate_max_m"] = float(old_metrics["ate_max_m_mean"])
        diagnostics["new_ate_max_m"] = float(new_metrics["ate_max_m_mean"])
        new_spikes = excess_spike_summary(new_rows, sequence, args.spike_excess_m)
        diagnostics["new_excess_spike_count"] = new_spikes["spike_count"]
        diagnostics["new_maximum_excess_error_m"] = new_spikes["maximum_excess_error_m"]
        sequences[sequence] = diagnostics
        for repeat in range(1, 4):
            metadata_path = (
                args.new_runs / sequence / METHOD / f"repeat_{repeat:02d}" / "run_metadata.txt"
            )
            metadata = parse_metadata(metadata_path)
            new_contract_rows.append({
                "sequence": sequence,
                "repeat": repeat,
                "completion": metadata.get("completion"),
                "trajectory_lines": int(metadata["trajectory_lines"]),
                "invalid_count": int(metadata["invalid_count"]),
                "nonmonotonic_count": int(metadata["nonmonotonic_count"]),
                "maximum_output_gap_s": float(metadata["maximum_output_gap_s"]),
                "excessive_output_gap_count": int(metadata["excessive_output_gap_count"]),
                "algorithm_binary_sha256": metadata["algorithm_binary_sha256"],
            })

    all_spikes = sum(int(row["fastlivo_excess_spike_count"]) for row in sequences.values())
    captured_spikes = sum(int(row["spikes_within_post_gap_window"]) for row in sequences.values())
    new_spikes = sum(int(row["new_excess_spike_count"]) for row in sequences.values())
    payload = {
        "schema_version": 1,
        "diagnosis": (
            "FAST-LIVO2 used processing-completion ros::Time::now() for a LiDAR state. "
            "Offline playback stalls therefore assigned late timestamps and caused GT-time mismatch spikes."
        ),
        "fix": "publish /aft_mapped_to_init with LidarMeasures.last_lio_update_time",
        "parameters": {
            "spike_definition": "FAST-LIVO2 error minus median error of the other three methods",
            "spike_excess_m": args.spike_excess_m,
            "gap_threshold_s": args.gap_threshold_s,
            "post_gap_window_s": args.post_gap_window_s,
        },
        "archived_representative_run_result": {
            "spike_count": all_spikes,
            "spikes_within_post_gap_window": captured_spikes,
            "capture_ratio": captured_spikes / all_spikes if all_spikes else 1.0,
        },
        "new_run_contract": {
            "run_count": len(new_contract_rows),
            "total_excessive_output_gaps": sum(
                int(row["excessive_output_gap_count"]) for row in new_contract_rows
            ),
            "maximum_output_gap_s": max(
                float(row["maximum_output_gap_s"]) for row in new_contract_rows
            ),
            "invalid_count": sum(int(row["invalid_count"]) for row in new_contract_rows),
            "nonmonotonic_count": sum(int(row["nonmonotonic_count"]) for row in new_contract_rows),
            "runs": new_contract_rows,
        },
        "new_representative_run_result": {
            "spike_count": new_spikes,
        },
        "sequences": sequences,
        "provenance": {
            "old_aligned_samples": {"path": str(old_aligned_path.resolve()), "sha256": sha256(old_aligned_path)},
            "new_aligned_samples": {
                "path": str(new_aligned_path.resolve()),
                "sha256": sha256(new_aligned_path),
            },
            "fixed_source": {"path": str(source.resolve()), "sha256": sha256(source)},
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(args.output.resolve()),
        "old_spikes": all_spikes,
        "captured_spikes": captured_spikes,
        "new_spikes": new_spikes,
        "new_total_gaps": payload["new_run_contract"]["total_excessive_output_gaps"],
        "new_max_gap_s": payload["new_run_contract"]["maximum_output_gap_s"],
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
