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
    "lightning_previous": "#85c1e9",
    "lightning_current": "#2471a3",
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


def subplot_grid(item_count: int) -> tuple[int, int]:
    columns = 2 if item_count > 1 else 1
    rows = (item_count + columns - 1) // columns
    return rows, columns


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
                    "hba_runs": backend.get("hba_runs"),
                    "hba_accepted": backend.get("hba_accepted"),
                    "config_sha256": metadata.get("config_sha256"),
                    "binary_sha256": metadata.get("algorithm_binary_sha256"),
                    "completion": metadata.get("completion"),
                    "algorithm_rc": metadata.get("algorithm_rc"),
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


def validate_and_summarize(
    manifest: dict[str, Any], rows: list[dict[str, Any]], output_dir: Path
) -> None:
    validation = manifest.get("validation", {})
    primary_method = manifest.get("primary_method")
    primary_rows = [row for row in rows if row["method"] == primary_method]
    if primary_method and len(primary_rows) != len(manifest["sequences"]):
        raise ValueError(f"primary method {primary_method!r} is missing one or more sequences")

    minimum_coverage = validation.get("minimum_gt_coverage")
    if minimum_coverage is not None:
        failed = [
            row["sequence"] for row in primary_rows
            if row["gt_coverage"] is None or row["gt_coverage"] < minimum_coverage
        ]
        if failed:
            raise ValueError(f"primary RTK coverage below {minimum_coverage}: {failed}")

    required_completion = validation.get("require_primary_completion")
    if required_completion is not None:
        failed = [row["sequence"] for row in primary_rows if row["completion"] != required_completion]
        if failed:
            raise ValueError(f"primary runs are incomplete: {failed}")

    required_rc = validation.get("require_primary_algorithm_rc")
    if required_rc is not None:
        failed = [row["sequence"] for row in primary_rows if float(row["algorithm_rc"]) != float(required_rc)]
        if failed:
            raise ValueError(f"primary runs returned a non-zero status: {failed}")

    config_hashes = sorted({row["config_sha256"] for row in primary_rows if row["config_sha256"]})
    binary_hashes = sorted({row["binary_sha256"] for row in primary_rows if row["binary_sha256"]})
    if validation.get("require_shared_primary_config") and len(config_hashes) != 1:
        raise ValueError(f"primary runs do not share one config hash: {config_hashes}")
    if validation.get("require_shared_primary_binary") and len(binary_hashes) != 1:
        raise ValueError(f"primary runs do not share one binary hash: {binary_hashes}")

    methods, _ = method_labels(manifest)
    aggregates: dict[str, dict[str, float]] = {}
    for method in methods:
        method_rows = [row for row in rows if row["method"] == method and row["valid_for_primary"]]
        aggregates[method] = {
            "mean_ate_rmse_m": float(np.mean([row["ate_rmse_m"] for row in method_rows])),
            "mean_rpe_10m_rmse_m": float(np.mean([row["rpe_10m_rmse_m"] for row in method_rows])),
            "mean_peak_rss_mb": float(np.mean([row["peak_rss_mb"] for row in method_rows])),
        }

    summary: dict[str, Any] = {
        "status": "passed",
        "sequence_count": len(manifest["sequences"]),
        "row_count": len(rows),
        "primary_method": primary_method,
        "primary_config_sha256": config_hashes,
        "primary_binary_sha256": binary_hashes,
        "aggregates": aggregates,
    }
    baseline_method = manifest.get("improvement_baseline_method")
    if primary_method and baseline_method:
        current = aggregates[primary_method]
        baseline = aggregates[baseline_method]
        summary["improvement_over_baseline_pct"] = {
            "mean_ate_rmse": 100.0 * (
                baseline["mean_ate_rmse_m"] - current["mean_ate_rmse_m"]
            ) / baseline["mean_ate_rmse_m"],
            "mean_rpe_10m_rmse": 100.0 * (
                baseline["mean_rpe_10m_rmse_m"] - current["mean_rpe_10m_rmse_m"]
            ) / baseline["mean_rpe_10m_rmse_m"],
            "mean_peak_rss": 100.0 * (
                baseline["mean_peak_rss_mb"] - current["mean_peak_rss_mb"]
            ) / baseline["mean_peak_rss_mb"],
        }
    with (output_dir / "validation_summary.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2)


def row_index(rows: list[dict[str, Any]]) -> dict[tuple[str, str], dict[str, Any]]:
    return {(row["sequence"], row["method"]): row for row in rows}


def plot_accuracy(manifest: dict[str, Any], rows: list[dict[str, Any]], output_dir: Path) -> None:
    sequences = list(manifest["sequences"])
    methods, display = method_labels(manifest)
    lookup = row_index(rows)
    x = np.arange(len(sequences), dtype=float)
    width = min(0.8 / max(len(methods), 1), 0.24)
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
    rows, columns = subplot_grid(len(sequences))
    fig, axes = plt.subplots(rows, columns, figsize=(7.2 * columns, 4.8 * rows), constrained_layout=True)
    axes_array = np.atleast_1d(axes).flat
    for axis, sequence in zip(axes_array, sequences):
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
    axes_array = np.atleast_1d(axes).flat
    for axis in list(axes_array)[len(sequences):]:
        axis.axis("off")
    np.atleast_1d(axes).flat[0].legend(frameon=False, fontsize=8)
    fig.savefig(output_dir / "trajectory_xy_overlays.png", dpi=180)
    plt.close(fig)


def plot_error_curves(manifest: dict[str, Any], output_dir: Path) -> None:
    methods, display = method_labels(manifest)
    sequences = list(manifest["sequences"])
    rows, columns = subplot_grid(len(sequences))
    fig, axes = plt.subplots(rows, columns, figsize=(7.2 * columns, 4.4 * rows), constrained_layout=True)
    axes_array = np.atleast_1d(axes).flat
    for axis, sequence in zip(axes_array, sequences):
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
    axes_array = np.atleast_1d(axes).flat
    for axis in list(axes_array)[len(sequences):]:
        axis.axis("off")
    np.atleast_1d(axes).flat[0].legend(frameon=False, fontsize=8)
    fig.savefig(output_dir / "translation_error_curves.png", dpi=180)
    plt.close(fig)


def plot_resources_and_loops(manifest: dict[str, Any], rows: list[dict[str, Any]], output_dir: Path) -> None:
    sequences = list(manifest["sequences"])
    methods, display = method_labels(manifest)
    lookup = row_index(rows)
    x = np.arange(len(sequences), dtype=float)
    width = min(0.8 / max(len(methods), 1), 0.24)
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
    loop_method = manifest.get("primary_method", "lightning_final")
    loop_rows = [lookup.get((sequence, loop_method), {}) for sequence in sequences]
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
    axis.set_title(f"{display[loop_method]} BTC accepted-loop precision / recall")
    axis.legend(frameon=False)
    axis.grid(axis="y", alpha=0.25)
    fig.savefig(output_dir / "loop_precision_recall.png", dpi=180)
    plt.close(fig)


def plot_improvements(manifest: dict[str, Any], rows: list[dict[str, Any]], output_dir: Path) -> None:
    current_method = manifest.get("primary_method")
    baseline_method = manifest.get("improvement_baseline_method")
    if not current_method or not baseline_method:
        return
    sequences = list(manifest["sequences"])
    lookup = row_index(rows)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8), constrained_layout=True)
    for axis, field, title in (
        (axes[0], "ate_rmse_m", "ATE improvement over previous backend"),
        (axes[1], "rpe_10m_rmse_m", "10 m RPE improvement over previous backend"),
    ):
        values = []
        for sequence in sequences:
            current = lookup[(sequence, current_method)][field]
            baseline = lookup[(sequence, baseline_method)][field]
            values.append(100.0 * (baseline - current) / baseline)
        colors = ["#229954" if value >= 0.0 else "#c0392b" for value in values]
        bars = axis.bar(np.arange(len(sequences)), values, color=colors)
        axis.axhline(0.0, color="#202020", linewidth=0.8)
        axis.set_title(title)
        axis.set_ylabel("improvement / %")
        axis.set_xticks(np.arange(len(sequences)), sequences, rotation=18, ha="right")
        axis.grid(axis="y", alpha=0.25)
        for bar, value in zip(bars, values):
            vertical = "bottom" if value >= 0.0 else "top"
            offset = 1.2 if value >= 0.0 else -1.2
            axis.text(bar.get_x() + bar.get_width() / 2, value + offset, f"{value:.1f}%", ha="center", va=vertical, fontsize=8)
    fig.savefig(output_dir / "improvement_vs_previous.png", dpi=180)
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
    validate_and_summarize(manifest, rows, args.output_dir)
    plot_accuracy(manifest, rows, args.output_dir)
    plot_trajectories(manifest, args.output_dir)
    plot_error_curves(manifest, args.output_dir)
    plot_resources_and_loops(manifest, rows, args.output_dir)
    plot_improvements(manifest, rows, args.output_dir)
    plot_ablations(manifest, args.output_dir)
    print(args.output_dir)


if __name__ == "__main__":
    main()
