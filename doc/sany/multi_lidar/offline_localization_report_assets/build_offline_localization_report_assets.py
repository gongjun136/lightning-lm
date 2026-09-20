#!/usr/bin/env python3
"""Build figures and summary metrics for the SANY offline localization report."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[4]
ASSET_DIR = Path(__file__).resolve().parent
RUN_ROOT = REPO / "runs" / "sany_offline_loc_full_20260713"
SLAM_RUN = RUN_ROOT / "full_20260713_slam_map"
LOC_RUN = RUN_ROOT / "full_20260713_localization_final"


def read_metadata(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def read_tum(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        values = [float(item) for item in line.split()]
        if len(values) != 8:
            continue
        rows.append(values)
    return pd.DataFrame(rows, columns=["timestamp", "x", "y", "z", "qx", "qy", "qz", "qw"])


def numeric_stats(values: pd.Series) -> dict[str, float]:
    return {
        "count": int(values.count()),
        "mean": float(values.mean()),
        "median": float(values.median()),
        "p95": float(values.quantile(0.95)),
        "max": float(values.max()),
    }


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 180,
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#D7DEE8",
            "axes.labelcolor": "#253041",
            "axes.titlecolor": "#172033",
            "xtick.color": "#44546A",
            "ytick.color": "#44546A",
            "grid.color": "#E8EDF3",
            "grid.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def savefig(name: str) -> None:
    plt.tight_layout()
    plt.savefig(ASSET_DIR / name, bbox_inches="tight")
    plt.close()


def plot_trajectory(reference: pd.DataFrame, estimate: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ref = reference.copy()
    est = estimate.copy()
    x0, y0 = ref.iloc[0][["x", "y"]]
    ref["x_rel"] = ref["x"] - x0
    ref["y_rel"] = ref["y"] - y0
    est["x_rel"] = est["x"] - x0
    est["y_rel"] = est["y"] - y0
    ax.plot(ref["x_rel"], ref["y_rel"], color="#2563EB", linewidth=2.0, label="Offline SLAM reference")
    ax.plot(est["x_rel"], est["y_rel"], color="#F59E0B", linewidth=1.7, linestyle="--", label="Offline localization")
    ax.scatter(ref["x_rel"].iloc[0], ref["y_rel"].iloc[0], color="#16A34A", s=38, zorder=3, label="Start")
    ax.scatter(ref["x_rel"].iloc[-1], ref["y_rel"].iloc[-1], color="#DC2626", s=38, zorder=3, label="End")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)
    ax.set_title("XY trajectory overlay")
    ax.set_xlabel("Relative X (m)")
    ax.set_ylabel("Relative Y (m)")
    ax.legend(loc="best", frameon=False)
    savefig("trajectory_xy_overlay.png")


def plot_errors(errors: pd.DataFrame) -> None:
    elapsed = errors["timestamp"] - errors["timestamp"].iloc[0]
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ax.plot(elapsed, errors["error_xy"], color="#2563EB", linewidth=1.5, label="XY error")
    ax.plot(elapsed, errors["error_3d"], color="#F59E0B", linewidth=1.3, alpha=0.85, label="3D error")
    ax.axhline(errors["error_xy"].quantile(0.95), color="#2563EB", linestyle=":", linewidth=1.2, label="XY P95")
    ax.set_title("Trajectory error over sensor time")
    ax.set_xlabel("Elapsed sensor time (s)")
    ax.set_ylabel("Error (m)")
    ax.grid(True)
    ax.legend(loc="upper right", frameon=False)
    savefig("trajectory_error_timeseries.png")

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.hist(errors["error_xy"], bins=32, color="#60A5FA", edgecolor="#1E3A8A", alpha=0.88)
    ax.axvline(errors["error_xy"].mean(), color="#0F172A", linewidth=1.4, label="Mean")
    ax.axvline(errors["error_xy"].quantile(0.95), color="#F59E0B", linewidth=1.4, label="P95")
    ax.set_title("XY error distribution")
    ax.set_xlabel("XY error (m)")
    ax.set_ylabel("Associated poses")
    ax.grid(axis="y")
    ax.legend(frameon=False)
    savefig("xy_error_distribution.png")


def plot_localization_quality(localization: pd.DataFrame) -> None:
    window = 25
    frame = localization["frame_index"]
    confidence = localization["confidence"].rolling(window, min_periods=1).median()
    processing = localization["processing_ms"].rolling(window, min_periods=1).median()
    fig, ax1 = plt.subplots(figsize=(8.8, 4.8))
    ax1.plot(frame, confidence, color="#2563EB", linewidth=1.7, label="Confidence, rolling median")
    ax1.set_xlabel("Localization frame")
    ax1.set_ylabel("NDT confidence", color="#2563EB")
    ax1.tick_params(axis="y", labelcolor="#2563EB")
    ax1.grid(True)
    ax2 = ax1.twinx()
    ax2.plot(frame, processing, color="#F59E0B", linewidth=1.3, label="Processing ms, rolling median")
    ax2.set_ylabel("Processing time (ms)", color="#B45309")
    ax2.tick_params(axis="y", labelcolor="#B45309")
    fig.legend(loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    ax1.set_title("Localization quality and per-frame cost")
    savefig("localization_quality_timeseries.png")


def plot_lidar_points(frame_stats: pd.DataFrame) -> None:
    lidar_cols = [column for column in frame_stats.columns if column.startswith("points_lidar_")]
    totals = frame_stats[lidar_cols].sum()
    means = frame_stats[lidar_cols].mean()
    labels = [column.replace("points_lidar_", "LiDAR ") for column in lidar_cols]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    bars = ax.bar(x, totals.values / 1e6, color=["#2563EB", "#F59E0B", "#64748B", "#EC4899"])
    ax.set_xticks(x, labels)
    ax.set_ylabel("Total input points (million)")
    ax.set_title("Input point contribution by LiDAR")
    ax.grid(axis="y")
    for bar, mean in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{mean:,.0f}/frame",
                ha="center", va="bottom", fontsize=9, color="#334155")
    savefig("lidar_point_contribution.png")

    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    elapsed = frame_stats["end_time"] - frame_stats["end_time"].iloc[0]
    ax.plot(elapsed, frame_stats["merged_points"].rolling(25, min_periods=1).median(),
            color="#2563EB", linewidth=1.6)
    ax.set_title("Merged point count stability")
    ax.set_xlabel("Elapsed sensor time (s)")
    ax.set_ylabel("Merged points per fused frame")
    ax.grid(True)
    savefig("merged_points_timeseries.png")


def plot_timing_and_resources(loc_timing: dict, slam_timing: dict, loc_res: dict, slam_res: dict) -> None:
    stage_order = [
        "Preprocess (Standard)",
        "Undistort Pcl",
        "ObsModel (Lidar Match)",
        "Incremental Mapping",
        "IVox Add Points",
    ]
    rows = []
    for stage in stage_order:
        loc_stage = next((item for item in loc_timing["stages"] if item["stage"] == stage), None)
        slam_stage = next((item for item in slam_timing["stages"] if item["stage"] == stage), None)
        if loc_stage and slam_stage:
            rows.append((stage, loc_stage["average_ms"], slam_stage["average_ms"]))
    labels = [row[0].replace(" (Standard)", "") for row in rows]
    y = np.arange(len(labels))
    height = 0.36
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ax.barh(y - height / 2, [row[1] for row in rows], height, color="#2563EB", label="Localization run")
    ax.barh(y + height / 2, [row[2] for row in rows], height, color="#F59E0B", label="SLAM map run")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Average stage time (ms)")
    ax.set_title("Core per-frame stage cost")
    ax.grid(axis="x")
    ax.legend(frameon=False)
    savefig("core_stage_timing.png")

    labels = ["SLAM map export", "Localization"]
    cpu_mean = [slam_res["mean_cpu_cores"], loc_res["mean_cpu_cores"]]
    cpu_p95 = [slam_res["p95_cpu_cores"], loc_res["p95_cpu_cores"]]
    rss_peak = [slam_res["peak_rss_mb"], loc_res["peak_rss_mb"]]
    x = np.arange(len(labels))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.5))
    ax1.bar(x - 0.18, cpu_mean, 0.36, color="#2563EB", label="Mean")
    ax1.bar(x + 0.18, cpu_p95, 0.36, color="#60A5FA", label="P95")
    ax1.set_xticks(x, labels, rotation=10)
    ax1.set_ylabel("CPU cores")
    ax1.set_title("CPU usage")
    ax1.grid(axis="y")
    ax1.legend(frameon=False)
    ax2.bar(x, rss_peak, color=["#F59E0B", "#EC4899"])
    ax2.set_xticks(x, labels, rotation=10)
    ax2.set_ylabel("Peak RSS (MB)")
    ax2.set_title("Peak memory")
    ax2.grid(axis="y")
    savefig("resource_usage.png")


def map_chunk_count(index_path: Path) -> int:
    count = 0
    for line in index_path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[0].lstrip("-").isdigit():
            count += 1
    return count


def main() -> int:
    setup_style()
    slam_meta = read_metadata(SLAM_RUN / "run_metadata.txt")
    loc_meta = read_metadata(LOC_RUN / "run_metadata.txt")
    loc_summary = json.loads((LOC_RUN / "results" / "localization_summary.json").read_text(encoding="utf-8"))
    loc_timing = json.loads((LOC_RUN / "results" / "processing_timing_summary.json").read_text(encoding="utf-8"))
    slam_timing = json.loads((SLAM_RUN / "results" / "processing_timing_summary.json").read_text(encoding="utf-8"))
    loc_resource = json.loads((LOC_RUN / "resource_summary.json").read_text(encoding="utf-8"))
    slam_resource = json.loads((SLAM_RUN / "resource_summary.json").read_text(encoding="utf-8"))
    localization = pd.read_csv(LOC_RUN / "results" / "localization_stats.csv")
    frame_stats = pd.read_csv(LOC_RUN / "results" / "frame_stats.csv")
    errors = pd.read_csv(LOC_RUN / "results" / "trajectory_reference_errors.csv")
    reference = read_tum(SLAM_RUN / "results" / "trajectory_slam.tum")
    estimate = read_tum(LOC_RUN / "results" / "trajectory_loc.tum")

    plot_trajectory(reference, estimate)
    plot_errors(errors)
    plot_localization_quality(localization)
    plot_lidar_points(frame_stats)
    plot_timing_and_resources(loc_timing, slam_timing, loc_resource, slam_resource)

    lidar_cols = [column for column in frame_stats.columns if column.startswith("points_lidar_")]
    completed_frames = int((frame_stats["partial"] == 0).sum())
    partial_frames = int((frame_stats["partial"] != 0).sum())
    summary = {
        "generated_assets": [
            "trajectory_xy_overlay.png",
            "trajectory_error_timeseries.png",
            "xy_error_distribution.png",
            "localization_quality_timeseries.png",
            "lidar_point_contribution.png",
            "merged_points_timeseries.png",
            "core_stage_timing.png",
            "resource_usage.png",
        ],
        "dataset": {
            "bag": loc_meta["bag"],
            "sensor_duration_s": float(loc_meta["sensor_duration_s"]),
            "primary_lidar_topic": loc_meta["primary_lidar_topic"],
            "payload_sha256": loc_meta["bag_payload_sha256"],
        },
        "map_export": {
            "algorithm_rc": int(slam_meta["algorithm_rc"]),
            "trajectory_lines": int(slam_meta["trajectory_lines"]),
            "last_stamp": float(slam_meta["last_stamp"]),
            "expected_last_lidar_end_s": float(slam_meta["expected_last_lidar_end_s"]),
            "tail_delta_s": float(slam_meta["expected_last_lidar_end_s"]) - float(slam_meta["last_stamp"]),
            "map_chunk_count_index": map_chunk_count(SLAM_RUN / "data" / "new_map" / "index.txt"),
            "map_chunk_count_metadata": int(slam_meta["map_chunk_count"]),
            "global_map_points": int(slam_meta["global_map_points"]),
            "wall_time_s": float(slam_meta["wall_time_s"]),
            "max_output_gap_s": float(slam_meta["maximum_output_gap_s"]),
        },
        "localization": {
            "algorithm_rc": int(loc_meta["algorithm_rc"]),
            "trajectory_lines": int(loc_meta["trajectory_lines"]),
            "completion": loc_meta["completion"],
            "valid_frames": int(loc_summary["localization"]["valid_frames"]),
            "valid_ratio": float(loc_summary["localization"]["valid_ratio"]),
            "associated_pairs": int(loc_summary["trajectory"]["associated_pairs"]),
            "xy_error": loc_summary["trajectory"]["error_xy"],
            "error_3d": loc_summary["trajectory"]["error_3d"],
            "confidence": loc_summary["localization"]["confidence"],
            "processing_ms": loc_summary["localization"]["processing_ms"],
            "loc_odom_delta": loc_summary["localization"]["loc_odom_delta"],
            "match_iterations": loc_summary["localization"]["match_iterations"],
            "active_map_chunks": loc_summary["localization"]["active_map_chunks"],
            "wall_time_s": float(loc_meta["wall_time_s"]),
            "realtime_factor": float(loc_meta["realtime_factor"]),
            "processing_speed_x": float(loc_meta["processing_speed_x"]),
            "max_output_gap_s": float(loc_meta["maximum_output_gap_s"]),
            "nonmonotonic_count": int(loc_meta["nonmonotonic_count"]),
            "invalid_count": int(loc_meta["invalid_count"]),
        },
        "multi_lidar_sync": {
            "fused_frames": int(len(frame_stats)),
            "complete_frames": completed_frames,
            "partial_frames": partial_frames,
            "complete_ratio": completed_frames / len(frame_stats),
            "mean_merged_points": float(frame_stats["merged_points"].mean()),
            "p95_merged_points": float(frame_stats["merged_points"].quantile(0.95)),
            "total_points_by_lidar": {column.replace("points_lidar_", "lidar_"): int(frame_stats[column].sum()) for column in lidar_cols},
            "mean_points_by_lidar": {column.replace("points_lidar_", "lidar_"): float(frame_stats[column].mean()) for column in lidar_cols},
        },
        "resources": {
            "slam": slam_resource,
            "localization": loc_resource,
        },
        "timing": {
            "slam": slam_timing["end_to_end"],
            "localization": loc_timing["end_to_end"],
        },
    }
    (ASSET_DIR / "summary_metrics.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
