#!/usr/bin/env python3
"""Build frontend pure-compute and Lightning-LM detailed timing figures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SEQUENCES = ("Grass02", "Outdoor04", "Z-Rough-Road01", "Dark01")
METHODS = ("lightning_lm", "fastlio", "fastlivo2_lio", "voxel_slam_frontend")
DISPLAY = {
    "lightning_lm": "Lightning-LM",
    "fastlio": "FAST-LIO",
    "fastlivo2_lio": "FAST-LIVO2 LIO",
    "voxel_slam_frontend": "Voxel-SLAM LIO",
    "preprocess_ms": "LiDAR preprocessing",
    "imu_undistort_ms": "IMU propagation / undistortion",
    "downsample_ms": "Downsampling",
    "match_setup_ms": "Matching setup",
    "scan_match_ms": "Iterated scan matching",
    "map_update_ms": "Map update",
}
COLORS = {
    "lightning_lm": "#0072B2",
    "fastlio": "#E69F00",
    "fastlivo2_lio": "#009E73",
    "voxel_slam_frontend": "#CC79A7",
}
STAGE_COLORS = {
    "preprocess_ms": "#D7EAF5",
    "imu_undistort_ms": "#A9D2E8",
    "downsample_ms": "#78B7D8",
    "match_setup_ms": "#4899C5",
    "scan_match_ms": "#1479AD",
    "map_update_ms": "#004E7A",
}


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    formal = repo / "runs" / "formal_report_20260717"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis-dir", type=Path,
        default=formal / "analysis" / "m3dgr_frontend_compute",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=repo / "doc" / "assets" / "formal_report_20260717",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row: dict[str, str], key: str) -> float:
    return float(row[key])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 180,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.6,
    })


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def frontend_summary(
    sequence_rows: list[dict[str, str]],
    overall_rows: list[dict[str, str]],
    output: Path,
) -> None:
    sequence_lookup = {(row["sequence"], row["method"]): row for row in sequence_rows}
    overall_lookup = {row["method"]: row for row in overall_rows}
    fig, axes = plt.subplots(2, 2, figsize=(14.2, 8.2))

    x = np.arange(len(SEQUENCES), dtype=float)
    width = 0.8 / len(METHODS)
    for index, method in enumerate(METHODS):
        values = [number(sequence_lookup[(sequence, method)], "total_ms_mean_repeat_mean") for sequence in SEQUENCES]
        errors = [number(sequence_lookup[(sequence, method)], "total_ms_mean_repeat_std") for sequence in SEQUENCES]
        offset = (index - (len(METHODS) - 1) / 2) * width
        axes[0, 0].bar(x + offset, values, width * 0.92, yerr=errors, capsize=2,
                       color=COLORS[method], label=DISPLAY[method])
    axes[0, 0].set_xticks(x, SEQUENCES, rotation=15, ha="right")
    axes[0, 0].set_ylabel("Pure compute time per LiDAR frame (ms)")
    axes[0, 0].set_title("Mean compute time by sequence")
    axes[0, 0].set_ylim(bottom=0)
    axes[0, 0].legend(frameon=False, ncol=2)

    quantiles = ("total_ms_median", "total_ms_p95", "total_ms_p99")
    quantile_labels = ("Median", "P95", "P99")
    x_methods = np.arange(len(METHODS), dtype=float)
    q_width = 0.23
    for index, (key, label) in enumerate(zip(quantiles, quantile_labels, strict=True)):
        values = [number(overall_lookup[method], key) for method in METHODS]
        axes[0, 1].bar(x_methods + (index - 1) * q_width, values, q_width * 0.94,
                       label=label, color=("#A6CEE3", "#4C9DCC", "#0B5D8C")[index])
    axes[0, 1].axhline(100.0, color="#D55E00", linestyle="--", linewidth=1.2, label="10 Hz deadline")
    axes[0, 1].set_xticks(x_methods, [DISPLAY[method] for method in METHODS], rotation=15, ha="right")
    axes[0, 1].set_ylabel("Compute time (ms)")
    axes[0, 1].set_title("Pooled latency distribution")
    axes[0, 1].set_ylim(bottom=0)
    axes[0, 1].legend(frameon=False, ncol=2)

    core = np.asarray([number(overall_lookup[method], "core_update_ms_mean") for method in METHODS])
    preprocess = np.asarray([number(overall_lookup[method], "preprocess_ms_mean") for method in METHODS])
    axes[1, 0].bar(x_methods, preprocess, color="#9ECAE1", label="LiDAR preprocessing")
    axes[1, 0].bar(x_methods, core, bottom=preprocess, color="#2171B5", label="Synchronized LIO core")
    axes[1, 0].set_xticks(x_methods, [DISPLAY[method] for method in METHODS], rotation=15, ha="right")
    axes[1, 0].set_ylabel("Mean compute time (ms)")
    axes[1, 0].set_title("Mean time composition")
    axes[1, 0].set_ylim(bottom=0)
    axes[1, 0].legend(frameon=False)

    throughput = np.asarray([number(overall_lookup[method], "effective_throughput_hz") for method in METHODS])
    utilization = np.asarray([number(overall_lookup[method], "compute_utilization") * 100.0 for method in METHODS])
    bars = axes[1, 1].bar(x_methods, throughput, color=[COLORS[method] for method in METHODS])
    axes[1, 1].set_xticks(x_methods, [DISPLAY[method] for method in METHODS], rotation=15, ha="right")
    axes[1, 1].set_ylabel("Equivalent compute throughput (frames/s)")
    axes[1, 1].set_title("Compute capacity and 10 Hz budget use")
    axes[1, 1].set_ylim(bottom=0)
    for bar, value in zip(bars, utilization, strict=True):
        axes[1, 1].text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                        f"{value:.1f}% budget", ha="center", va="bottom", fontsize=8)

    fig.suptitle("M3DGR frontend LIO: pure algorithm compute performance", fontsize=13, fontweight="bold")
    fig.text(
        0.5, 0.006,
        "Boundary: LiDAR preprocessing + synchronized LIO update; excludes playback wait, ROS output, file I/O and backend.",
        ha="center", fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.96), h_pad=2.2, w_pad=2.0)
    save(fig, output)


def lightning_stage_composition(rows: list[dict[str, str]], output: Path) -> None:
    rows = [row for row in rows if row["sequence"] in SEQUENCES]
    lookup = {(row["sequence"], row["stage"]): row for row in rows}
    stages = tuple(STAGE_COLORS)
    x = np.arange(len(SEQUENCES), dtype=float)
    bottom = np.zeros(len(SEQUENCES))
    fig, axis = plt.subplots(figsize=(10.8, 5.2))
    for stage in stages:
        values = np.asarray([number(lookup[(sequence, stage)], "stage_ms_mean") for sequence in SEQUENCES])
        axis.bar(x, values, bottom=bottom, color=STAGE_COLORS[stage], label=DISPLAY[stage])
        bottom += values
    axis.set_xticks(x, SEQUENCES)
    axis.set_ylabel("Mean compute time per LiDAR frame (ms)")
    axis.set_title("Lightning-LM detailed compute-time composition", fontweight="bold", fontsize=12)
    axis.set_ylim(bottom=0)
    axis.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    for position, total in enumerate(bottom):
        axis.text(position, total, f"{total:.2f} ms", ha="center", va="bottom", fontsize=8)
    fig.text(0.5, 0.008, "Pooled included frames from 3 runs; first 10 tracking frames per run excluded.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 0.82, 1))
    save(fig, output)


def lightning_distribution(frame_rows: list[dict[str, str]], output: Path) -> None:
    values = [
        np.asarray([
            number(row, "total_ms") for row in frame_rows
            if row["method"] == "lightning_lm" and row["sequence"] == sequence and row["included"] == "True"
        ])
        for sequence in SEQUENCES
    ]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
    violin = axes[0].violinplot(values, showmeans=False, showmedians=True, showextrema=False)
    for body in violin["bodies"]:
        body.set_facecolor(COLORS["lightning_lm"])
        body.set_edgecolor("#004E7A")
        body.set_alpha(0.65)
    violin["cmedians"].set_color("#D55E00")
    axes[0].boxplot(values, widths=0.13, showfliers=False, patch_artist=True,
                    boxprops={"facecolor": "white", "alpha": 0.8},
                    medianprops={"color": "#D55E00"})
    axes[0].set_xticks(np.arange(1, len(SEQUENCES) + 1), SEQUENCES, rotation=12, ha="right")
    axes[0].set_ylabel("Total pure compute time (ms)")
    zoom_limit = max(float(np.percentile(item, 99.5)) for item in values) * 1.15
    axes[0].set_title("Per-frame distribution (zoomed to P99.5)")
    axes[0].set_ylim(0, zoom_limit)

    keys = ("Median", "P95", "P99", "Maximum")
    probabilities = (50, 95, 99, 100)
    x = np.arange(len(SEQUENCES), dtype=float)
    width = 0.19
    colors = ("#B3D8EA", "#65AED2", "#1F7EAF", "#004E7A")
    for index, (label, probability) in enumerate(zip(keys, probabilities, strict=True)):
        quantile_values = [np.percentile(item, probability) for item in values]
        axes[1].bar(x + (index - 1.5) * width, quantile_values, width * 0.94,
                    color=colors[index], label=label)
    axes[1].axhline(100.0, color="#D55E00", linestyle="--", linewidth=1.2, label="10 Hz deadline")
    axes[1].set_xticks(x, SEQUENCES, rotation=12, ha="right")
    axes[1].set_ylabel("Total pure compute time (ms)")
    axes[1].set_title("Latency quantiles and worst case")
    axes[1].set_ylim(bottom=0)
    axes[1].legend(frameon=False, ncol=2)

    fig.suptitle("Lightning-LM frontend: detailed latency distribution", fontsize=13, fontweight="bold")
    fig.text(
        0.5, 0.006,
        "Left: distribution zoomed to P99.5 for readability. Right: P50/P95/P99 and the untrimmed maximum.",
        ha="center", fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95), w_pad=2.5)
    save(fig, output)


def rolling_median(values: np.ndarray, window: int) -> np.ndarray:
    if len(values) < window:
        return values.copy()
    padded = np.pad(values, (window - 1, 0), mode="edge")
    return np.asarray([np.median(padded[index:index + window]) for index in range(len(values))])


def lightning_timeseries(frame_rows: list[dict[str, str]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 7.6), sharey=False)
    for axis, sequence in zip(axes.flat, SEQUENCES, strict=True):
        selected = [
            row for row in frame_rows
            if row["method"] == "lightning_lm" and row["sequence"] == sequence
            and row["repeat"] == "1" and row["included"] == "True"
        ]
        indices = np.asarray([number(row, "tracking_index") for row in selected])
        total = np.asarray([number(row, "total_ms") for row in selected])
        axis.plot(indices, total, color="#8FC4DF", linewidth=0.55, alpha=0.65, label="Per frame")
        axis.plot(indices, rolling_median(total, 50), color=COLORS["lightning_lm"], linewidth=1.25,
                  label="50-frame rolling median")
        axis.axhline(100.0, color="#D55E00", linestyle="--", linewidth=1.0, label="10 Hz deadline")
        axis.set_title(sequence)
        axis.set_xlabel("Tracking-frame index")
        axis.set_ylabel("Pure compute time (ms)")
        axis.set_ylim(bottom=0)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle("Lightning-LM frontend: compute-time stability over each sequence", fontsize=13, fontweight="bold")
    fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.945))
    fig.text(0.5, 0.006, "Predeclared representative run: repeat 1; raw timings retained behind the rolling median.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 0.9), h_pad=2.0, w_pad=2.0)
    save(fig, output)


def main() -> int:
    args = parse_args()
    validation_path = args.analysis_dir / "validation.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if validation.get("status") not in {"passed", "passed_with_warnings"}:
        raise SystemExit(f"analysis validation did not pass: {validation.get('status')}")

    inputs = {
        "sequence_summary": args.analysis_dir / "sequence_method_summary.csv",
        "overall_summary": args.analysis_dir / "overall_method_summary.csv",
        "lightning_stages": args.analysis_dir / "lightning_stage_summary.csv",
        "frame_timing": args.analysis_dir / "frame_timing.csv",
        "validation": validation_path,
        "analysis_manifest": args.analysis_dir / "analysis_manifest.json",
    }
    missing = [str(path) for path in inputs.values() if not path.is_file()]
    if missing:
        raise SystemExit(f"missing validated analysis inputs: {missing}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    configure_style()
    sequence_rows = read_csv(inputs["sequence_summary"])
    overall_rows = read_csv(inputs["overall_summary"])
    stage_rows = read_csv(inputs["lightning_stages"])
    frame_rows = read_csv(inputs["frame_timing"])

    figures = {
        "frontend_compute_summary": args.output_dir / "m3dgr_frontend_compute_summary.png",
        "lightning_stage_composition": args.output_dir / "lightning_lio_timing_stage_summary.png",
        "lightning_latency_distribution": args.output_dir / "lightning_lio_timing_distribution.png",
        "lightning_timing_timeseries": args.output_dir / "lightning_lio_timing_timeseries.png",
    }
    frontend_summary(sequence_rows, overall_rows, figures["frontend_compute_summary"])
    lightning_stage_composition(stage_rows, figures["lightning_stage_composition"])
    lightning_distribution(frame_rows, figures["lightning_latency_distribution"])
    lightning_timeseries(frame_rows, figures["lightning_timing_timeseries"])

    script = Path(__file__).resolve()
    manifest = {
        "schema_version": 1,
        "generator": {"path": str(script), "sha256": sha256(script)},
        "inputs": {name: {"path": str(path.resolve()), "sha256": sha256(path)} for name, path in inputs.items()},
        "figures": {name: {"path": str(path.resolve()), "sha256": sha256(path)} for name, path in figures.items()},
    }
    manifest_path = args.output_dir / "frontend_compute_figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"figures": len(figures), "manifest": str(manifest_path)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
