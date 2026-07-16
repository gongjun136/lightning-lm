#!/usr/bin/env python3
"""Build reproducible figures for the multi-LiDAR relocalization report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BLUE = "#2471a3"
GOLD = "#d68910"
INK = "#263238"
GREY = "#7f8c8d"
LIGHT_GREY = "#d9e1e5"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def load_json(path: str | Path) -> dict:
    with Path(path).open("r", encoding="utf-8") as stream:
        return json.load(stream)


def load_tum(path: str | Path) -> np.ndarray:
    rows = []
    with Path(path).open("r", encoding="utf-8-sig") as stream:
        for line in stream:
            line = line.strip()
            if line and not line.startswith("#"):
                values = [float(value) for value in line.split()]
                if len(values) == 8 and np.all(np.isfinite(values)):
                    rows.append(values)
    if not rows:
        raise ValueError(f"no TUM poses in {path}")
    return np.asarray(rows, dtype=np.float64)


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "grid.color": LIGHT_GREY,
            "grid.linewidth": 0.7,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
        }
    )


def finish(figure: plt.Figure, path: Path) -> None:
    figure.tight_layout()
    figure.savefig(path, dpi=190, bbox_inches="tight")
    plt.close(figure)


def metric(evaluation: dict, name: str) -> float:
    if name == "ate":
        return float(evaluation["ate_translation"]["rmse_m"])
    return float(evaluation["rpe_translation"]["10m"]["rmse_m"])


def label_bars(axis: plt.Axes, bars, decimals: int = 3) -> None:
    for bar in bars:
        height = float(bar.get_height())
        axis.annotate(
            f"{height:.{decimals}f}",
            (bar.get_x() + bar.get_width() / 2.0, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            color=INK,
        )


def plot_20260701(manifest: dict, output_dir: Path) -> None:
    evaluation = load_json(manifest["sany_20260701_evaluation"])
    runs = evaluation["runs"]
    names = ["Map LIO", "Map backend", "Localization"]
    keys = ["mapping_frontend", "mapping_backend", "localization_final"]
    colors = [GREY, BLUE, GOLD]
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    for axis, field, title in zip(
        axes, ["ate", "rpe"], ["ATE RMSE", "10 m RPE RMSE"]
    ):
        values = [metric(runs[key], field) for key in keys]
        bars = axis.bar(names, values, color=colors, edgecolor=INK, linewidth=0.7)
        label_bars(axis, bars, decimals=4)
        axis.set_ylabel("translation error (m)")
        axis.set_title(title)
        axis.grid(axis="y")
        axis.set_axisbelow(True)
        axis.set_ylim(0.0, max(values) * 1.28)
    figure.suptitle("SANY 20260701 mapping and localization accuracy", fontsize=14, fontweight="bold")
    finish(figure, output_dir / "sany_20260701_accuracy.png")


def plot_reference_consistency(manifest: dict, output_dir: Path) -> None:
    comparisons = [
        ("LIO vs backend", "frontend_backend_evaluation", GREY),
        ("Voxel114 vs backend", "voxel_backend_evaluation", GOLD),
        ("Localization vs backend", "startup_evaluation", BLUE),
    ]
    datasets = ["data1", "data2"]
    figure, axes = plt.subplots(1, 2, figsize=(12.3, 4.5))
    x = np.arange(len(datasets), dtype=np.float64)
    width = 0.23
    for axis, field, title in zip(
        axes, ["ate", "rpe"], ["Fixed-scale SE(3) ATE RMSE", "10 m RPE RMSE"]
    ):
        for index, (label, key, color) in enumerate(comparisons):
            values = [
                metric(load_json(manifest["sany_20260716"][dataset][key]), field)
                for dataset in datasets
            ]
            bars = axis.bar(
                x + (index - 1) * width,
                values,
                width,
                label=label,
                color=color,
                edgecolor=INK,
                linewidth=0.6,
            )
            label_bars(axis, bars, decimals=4)
        axis.set_xticks(x, datasets)
        axis.set_ylabel("translation error (m)")
        axis.set_title(title)
        axis.grid(axis="y")
        axis.set_axisbelow(True)
    axes[0].legend(loc="upper left", frameon=False)
    figure.suptitle("SANY 20260716 three-trajectory consistency", fontsize=14, fontweight="bold")
    finish(figure, output_dir / "sany_three_trajectory_consistency.png")


def plot_trajectory_overlays(manifest: dict, output_dir: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12.0, 5.1))
    for axis, dataset in zip(axes, ["data1", "data2"]):
        config = manifest["sany_20260716"][dataset]
        evaluation = load_json(config["startup_evaluation"])
        reference = load_tum(config["reference_trajectory"])
        estimate = load_tum(config["startup_trajectory"])
        rotation = np.asarray(evaluation["se3_rotation"], dtype=np.float64)
        translation = np.asarray(evaluation["se3_translation_m"], dtype=np.float64)
        aligned = (rotation @ estimate[:, 1:4].T).T + translation
        axis.plot(reference[:, 1], reference[:, 2], color=INK, linewidth=2.0, label="SLAM backend")
        axis.plot(aligned[:, 0], aligned[:, 1], color=BLUE, linewidth=1.3, label="Localization aligned")
        axis.scatter(reference[0, 1], reference[0, 2], s=45, color=GOLD, edgecolor=INK, zorder=4, label="Start")
        axis.set_aspect("equal", adjustable="datalim")
        axis.set_xlabel("x (m)")
        axis.set_ylabel("y (m)")
        axis.set_title(dataset)
        axis.grid()
    axes[0].legend(loc="best", frameon=False)
    figure.suptitle("SANY 20260716 startup relocalization trajectory overlays", fontsize=14, fontweight="bold")
    finish(figure, output_dir / "sany_startup_trajectory_overlays.png")


def plot_relocalization_modes(manifest: dict, output_dir: Path) -> None:
    datasets = ["data1", "data2"]
    modes = [
        ("Startup", "startup_evaluation", BLUE),
        ("Forced loss", "recovery_evaluation", GOLD),
    ]
    figure, axes = plt.subplots(1, 3, figsize=(14.2, 4.4))
    x = np.arange(len(datasets), dtype=np.float64)
    width = 0.32
    for axis, field, title in zip(
        axes[:2], ["ate", "rpe"], ["ATE RMSE", "10 m RPE RMSE"]
    ):
        for index, (label, key, color) in enumerate(modes):
            values = [
                metric(load_json(manifest["sany_20260716"][dataset][key]), field)
                for dataset in datasets
            ]
            bars = axis.bar(
                x + (index - 0.5) * width,
                values,
                width,
                label=label,
                color=color,
                edgecolor=INK,
                linewidth=0.6,
            )
            label_bars(axis, bars, decimals=3)
        axis.set_xticks(x, datasets)
        axis.set_ylabel("translation error (m)")
        axis.set_title(title)
        axis.grid(axis="y")
        axis.set_axisbelow(True)
    latencies = [
        float(manifest["sany_20260716"][dataset]["recovery_latency_s"])
        for dataset in datasets
    ]
    bars = axes[2].bar(datasets, latencies, color=[BLUE, GOLD], edgecolor=INK, linewidth=0.7)
    label_bars(axes[2], bars, decimals=1)
    axes[2].set_ylabel("sensor time (s)")
    axes[2].set_title("Forced-loss recovery latency")
    axes[2].grid(axis="y")
    axes[2].set_axisbelow(True)
    axes[0].legend(loc="upper left", frameon=False)
    figure.suptitle("SANY 20260716 relocalization validation", fontsize=14, fontweight="bold")
    finish(figure, output_dir / "sany_relocalization_modes.png")


def plot_gravity_gate(manifest: dict, output_dir: Path) -> None:
    candidates = manifest["gravity_gate_candidates"]
    labels = [candidate["label"] for candidate in candidates]
    colors = [GREY if not candidate["accepted"] else BLUE for candidate in candidates]
    figure, axes = plt.subplots(1, 2, figsize=(10.7, 4.3))
    overlap = [float(candidate["map_overlap_ratio"]) for candidate in candidates]
    gravity = [float(candidate["gravity_alignment_cos"]) for candidate in candidates]
    bars = axes[0].bar(labels, overlap, color=colors, edgecolor=INK, linewidth=0.7)
    label_bars(axes[0], bars, decimals=3)
    axes[0].axhline(float(manifest["map_overlap_threshold"]), color=GOLD, linestyle="--", label="threshold")
    axes[0].set_ylim(0.0, 1.08)
    axes[0].set_ylabel("ratio")
    axes[0].set_title("Nearest-neighbor map overlap")
    axes[0].grid(axis="y")
    axes[0].legend(frameon=False)
    bars = axes[1].bar(labels, gravity, color=colors, edgecolor=INK, linewidth=0.7)
    for bar, value in zip(bars, gravity):
        axes[1].annotate(f"{value:.3f}", (bar.get_x() + bar.get_width() / 2.0, value),
                         xytext=(0, 4 if value >= 0 else -14), textcoords="offset points",
                         ha="center", va="bottom", fontsize=8)
    axes[1].axhline(float(manifest["gravity_alignment_threshold"]), color=GOLD, linestyle="--", label="threshold")
    axes[1].axhline(0.0, color=INK, linewidth=0.8)
    axes[1].set_ylim(-1.12, 1.12)
    axes[1].set_ylabel("cosine")
    axes[1].set_title("Map/odometry gravity alignment")
    axes[1].grid(axis="y")
    axes[1].legend(frameon=False)
    figure.suptitle("Data2 startup candidate gate evidence", fontsize=14, fontweight="bold")
    finish(figure, output_dir / "sany_gravity_gate_evidence.png")


def plot_m3dgr(manifest: dict, output_dir: Path) -> None:
    sequences = ["Grass02", "Outdoor04", "Dark01", "Z-Rough-Road01"]
    methods = [
        ("Optimized 20260715", "optimized_20260715", GREY),
        ("Regression 20260716", "regression_20260716", BLUE),
        ("ws_voxel_slam", "ws_voxel_slam", GOLD),
    ]
    evaluations = {name: load_json(manifest["m3dgr_evaluations"][name]) for name in sequences}
    figure, axes = plt.subplots(2, 1, figsize=(12.2, 8.0), sharex=True)
    x = np.arange(len(sequences), dtype=np.float64)
    width = 0.24
    for axis, field, title in zip(
        axes, ["ate", "rpe"], ["ATE RMSE", "10 m RPE RMSE"]
    ):
        for index, (label, run_key, color) in enumerate(methods):
            values = [metric(evaluations[sequence]["runs"][run_key], field) for sequence in sequences]
            bars = axis.bar(
                x + (index - 1) * width,
                values,
                width,
                label=label,
                color=color,
                edgecolor=INK,
                linewidth=0.6,
            )
            label_bars(axis, bars, decimals=3)
        axis.set_ylabel("translation error (m)")
        axis.set_title(title)
        axis.grid(axis="y")
        axis.set_axisbelow(True)
    axes[1].set_xticks(x, sequences)
    figure.suptitle(
        "M3DGR single-LiDAR regression (fixed-scale SE(3))",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.958),
        frameon=False,
        ncol=3,
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.91))
    figure.savefig(
        output_dir / "m3dgr_single_lidar_regression.png",
        dpi=190,
        bbox_inches="tight",
    )
    plt.close(figure)


def main() -> None:
    args = parse_args()
    manifest = load_json(args.manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    configure_style()
    plot_20260701(manifest, args.output_dir)
    plot_reference_consistency(manifest, args.output_dir)
    plot_trajectory_overlays(manifest, args.output_dir)
    plot_relocalization_modes(manifest, args.output_dir)
    plot_gravity_gate(manifest, args.output_dir)
    plot_m3dgr(manifest, args.output_dir)
    chart_map = {
        "surface": "Markdown technical report with static PNG figures",
        "palette": "blue/gold roots plus neutral grey; edge and labels remain legible without color",
        "figures": [
            "sany_20260701_accuracy.png",
            "sany_three_trajectory_consistency.png",
            "sany_startup_trajectory_overlays.png",
            "sany_relocalization_modes.png",
            "sany_gravity_gate_evidence.png",
            "m3dgr_single_lidar_regression.png",
        ],
        "manifest": str(args.manifest),
    }
    with (args.output_dir / "chart_map.json").open("w", encoding="utf-8") as stream:
        json.dump(chart_map, stream, ensure_ascii=False, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
