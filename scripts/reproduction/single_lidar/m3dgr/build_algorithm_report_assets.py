#!/usr/bin/env python3
"""Build reproducible figures and metric tables for the M3DGR backend report."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from evaluate_backend_runs import align_se3, interpolate_positions, load_tum


COLORS = {
    "legacy": "#7f8c8d",
    "lightning_final": "#2471a3",
    "ws_voxel_slam": "#d68910",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def load_entry(entry: dict[str, Any]) -> dict[str, Any]:
    report = load_json(Path(entry["evaluation"]))
    return report["runs"][entry["run_label"]]


def method_labels(manifest: dict[str, Any]) -> tuple[list[str], dict[str, str]]:
    order = manifest.get("methods_order", ["legacy", "lightning_final", "ws_voxel_slam"])
    display = manifest.get(
        "methods_display",
        {"legacy": "Legacy", "lightning_final": "Lightning BA+BTC+HBA", "ws_voxel_slam": "ws_voxel_slam"},
    )
    return order, display


def collect_rows(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for sequence, sequence_cfg in manifest["sequences"].items():
        for method, entry in sequence_cfg["methods"].items():
            run = load_entry(entry)
            loop = run.get("loop") or {}
            resource = run.get("resource") or {}
            metadata = run.get("run_metadata") or {}
            backend = run.get("backend") or {}
            rows.append(
                {
                    "sequence": sequence,
                    "method": method,
                    "valid_for_primary": bool(entry.get("valid_for_primary", True)),
                    "validity_note": entry.get("validity_note", ""),
                    "ate_rmse_m": run["ate_translation"]["rmse_m"],
                    "ate_median_m": run["ate_translation"]["median_m"],
                    "rpe_10m_rmse_m": (run.get("rpe_translation") or {}).get("10m", {}).get("rmse_m"),
                    "gt_coverage": run.get("gt_pose_coverage"),
                    "associated_span_s": run.get("associated_time_span_s"),
                    "wall_time_s": metadata.get("wall_time_s"),
                    "peak_rss_mb": resource.get("peak_rss_mb"),
                    "loop_precision": loop.get("precision"),
                    "loop_recall": loop.get("recall"),
                    "loops_accepted": backend.get("loops_accepted"),
                    "loops_applied": backend.get("loops_applied"),
                }
            )
    return rows


def write_tables(rows: list[dict[str, Any]], output_dir: Path) -> None:
    with (output_dir / "metrics_summary.json").open("w", encoding="utf-8") as stream:
        json.dump(rows, stream, ensure_ascii=False, indent=2)
    with (output_dir / "metrics_summary.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def row_index(rows: list[dict[str, Any]]) -> dict[tuple[str, str], dict[str, Any]]:
    return {(row["sequence"], row["method"]): row for row in rows}


def plot_accuracy(manifest: dict[str, Any], rows: list[dict[str, Any]], output_dir: Path) -> None:
    sequences = list(manifest["sequences"])
    methods, display = method_labels(manifest)
    lookup = row_index(rows)
    x = np.arange(len(sequences), dtype=float)
    width = 0.24
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8), constrained_layout=True)
    for method_i, method in enumerate(methods):
        offset = (method_i - (len(methods) - 1) / 2.0) * width
        for axis, field, title in (
            (axes[0], "ate_rmse_m", "ATE RMSE (fixed-scale SE(3))"),
            (axes[1], "rpe_10m_rmse_m", "10 m translational RPE RMSE"),
        ):
            values = []
            for sequence in sequences:
                row = lookup.get((sequence, method))
                values.append(
                    np.nan
                    if row is None or not row["valid_for_primary"] or row[field] is None
                    else row[field]
                )
            axis.bar(x + offset, values, width, color=COLORS.get(method), label=display[method])
            axis.set_title(title)
            axis.set_ylabel("metres (log scale)")
            axis.set_yscale("log")
            axis.grid(axis="y", which="both", alpha=0.25)
    for axis in axes:
        axis.set_xticks(x, sequences, rotation=18, ha="right")
    axes[0].legend(frameon=False, fontsize=9)
    fig.savefig(output_dir / "accuracy_ate_rpe.png", dpi=180)
    plt.close(fig)


def aligned_series(gt_path: Path, trajectory_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gt = load_tum(gt_path, ground_truth=True)
    estimate = load_tum(trajectory_path, ground_truth=False)
    interpolated, valid = interpolate_positions(estimate, gt[:, 0], 0.2)
    valid &= np.all(np.isfinite(interpolated), axis=1)
    rotation, translation = align_se3(interpolated[valid], gt[valid, 1:4])
    aligned_all = (rotation @ estimate[:, 1:4].T).T + translation
    aligned_associated = (rotation @ interpolated[valid].T).T + translation
    errors = np.linalg.norm(aligned_associated - gt[valid, 1:4], axis=1)
    return aligned_all, gt[valid, 1:4], errors


def plot_trajectories(manifest: dict[str, Any], output_dir: Path) -> None:
    methods, display = method_labels(manifest)
    sequences = list(manifest["sequences"])
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for axis, sequence in zip(axes.flat, sequences):
        cfg = manifest["sequences"][sequence]
        gt = load_tum(Path(cfg["gt"]), ground_truth=True)
        axis.plot(gt[:, 1], gt[:, 2], color="#202020", linewidth=1.4, label="RTK")
        for method in methods:
            entry = cfg["methods"].get(method)
            if entry is None:
                continue
            aligned, _, _ = aligned_series(Path(cfg["gt"]), Path(entry["trajectory"]))
            suffix = " (partial)" if not entry.get("valid_for_primary", True) else ""
            axis.plot(
                aligned[:, 0], aligned[:, 1], linewidth=1.0,
                color=COLORS.get(method), alpha=0.85,
                linestyle="--" if suffix else "-", label=display[method] + suffix,
            )
        axis.set_title(sequence)
        axis.set_aspect("equal", adjustable="datalim")
        axis.grid(alpha=0.2)
        axis.set_xlabel("x / m")
        axis.set_ylabel("y / m")
    for axis in axes.flat[len(sequences):]:
        axis.axis("off")
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.savefig(output_dir / "trajectory_xy_overlays.png", dpi=180)
    plt.close(fig)


def plot_error_curves(manifest: dict[str, Any], output_dir: Path) -> None:
    methods, display = method_labels(manifest)
    sequences = list(manifest["sequences"])
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for axis, sequence in zip(axes.flat, sequences):
        cfg = manifest["sequences"][sequence]
        for method in methods:
            entry = cfg["methods"].get(method)
            if entry is None or not entry.get("valid_for_primary", True):
                continue
            _, _, errors = aligned_series(Path(cfg["gt"]), Path(entry["trajectory"]))
            progress = np.linspace(0.0, 100.0, len(errors))
            axis.plot(progress, errors, color=COLORS.get(method), linewidth=1.0, label=display[method])
        axis.set_title(sequence)
        axis.set_xlabel("associated trajectory progress / %")
        axis.set_ylabel("translation error / m")
        axis.grid(alpha=0.2)
    for axis in axes.flat[len(sequences):]:
        axis.axis("off")
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.savefig(output_dir / "translation_error_curves.png", dpi=180)
    plt.close(fig)


def plot_resources_and_loops(manifest: dict[str, Any], rows: list[dict[str, Any]], output_dir: Path) -> None:
    sequences = list(manifest["sequences"])
    methods, display = method_labels(manifest)
    lookup = row_index(rows)
    x = np.arange(len(sequences), dtype=float)
    width = 0.24
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8), constrained_layout=True)
    for method_i, method in enumerate(methods):
        offset = (method_i - (len(methods) - 1) / 2.0) * width
        wall = [lookup.get((seq, method), {}).get("wall_time_s") for seq in sequences]
        memory = [lookup.get((seq, method), {}).get("peak_rss_mb") for seq in sequences]
        wall = [np.nan if value is None else value for value in wall]
        memory = [np.nan if value is None else value for value in memory]
        axes[0].bar(x + offset, wall, width, color=COLORS.get(method), label=display[method])
        axes[1].bar(x + offset, memory, width, color=COLORS.get(method), label=display[method])
    axes[0].set_title("Wall time (different playback contracts)")
    axes[0].set_ylabel("seconds")
    axes[1].set_title("Peak resident memory")
    axes[1].set_ylabel("MiB")
    for axis in axes:
        axis.set_xticks(x, sequences, rotation=18, ha="right")
        axis.grid(axis="y", alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    fig.savefig(output_dir / "runtime_memory.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(10.5, 5.2), constrained_layout=True)
    loop_rows = [lookup.get((sequence, "lightning_final"), {}) for sequence in sequences]
    precision = [0.0 if row.get("loop_precision") is None else row["loop_precision"] for row in loop_rows]
    recall = [0.0 if row.get("loop_recall") is None else row["loop_recall"] for row in loop_rows]
    x_loop = np.arange(len(sequences), dtype=float)
    axis.bar(x_loop - 0.18, precision, 0.36, color="#2471a3", label="accepted-loop precision")
    axis.bar(x_loop + 0.18, recall, 0.36, color="#d68910", label="RTK episode recall")
    for index, row in enumerate(loop_rows):
        if row.get("loop_precision") is None:
            axis.scatter(index - 0.18, 0.015, color="#2471a3", marker="x", s=55)
            axis.text(index - 0.18, 0.035, "N/A", ha="center", va="bottom", fontsize=7)
    axis.set_ylim(0.0, 1.0)
    axis.set_xticks(x_loop, sequences, rotation=18, ha="right")
    axis.set_ylabel("score")
    axis.set_title("Lightning BTC accepted-loop precision / recall")
    axis.legend(frameon=False)
    axis.grid(axis="y", alpha=0.25)
    fig.savefig(output_dir / "loop_precision_recall.png", dpi=180)
    plt.close(fig)


def plot_ablations(manifest: dict[str, Any], output_dir: Path) -> None:
    ablations = manifest.get("ablations", [])
    if not ablations:
        return
    labels = [item["label"] for item in ablations]
    values = [item["ate_rmse_m"] for item in ablations]
    colors = [item.get("color", "#5d6d7e") for item in ablations]
    fig, axis = plt.subplots(figsize=(12.5, 5.2), constrained_layout=True)
    bars = axis.bar(np.arange(len(labels)), values, color=colors)
    axis.set_yscale("log")
    axis.set_ylabel("ATE RMSE / m (log scale)")
    axis.set_title("Key frontend/backend ablations and rejected variants")
    axis.set_xticks(np.arange(len(labels)), labels, rotation=24, ha="right")
    axis.grid(axis="y", which="both", alpha=0.25)
    for bar, value in zip(bars, values):
        axis.text(bar.get_x() + bar.get_width() / 2, value * 1.08, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
    fig.savefig(output_dir / "ablation_ate.png", dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    manifest = load_json(args.manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.titleweight": "bold"})
    rows = collect_rows(manifest)
    if not rows:
        raise SystemExit("manifest contains no runs")
    write_tables(rows, args.output_dir)
    plot_accuracy(manifest, rows, args.output_dir)
    plot_trajectories(manifest, args.output_dir)
    plot_error_curves(manifest, args.output_dir)
    plot_resources_and_loops(manifest, rows, args.output_dir)
    plot_ablations(manifest, args.output_dir)
    print(args.output_dir)


if __name__ == "__main__":
    main()
