#!/usr/bin/env python3
"""Build report figures for the controlled M3DGR backend experiment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


METHODS = ("frontend_only", "legacy_controlled", "local_ba_only", "loop_pgo", "full_backend")
LABELS = ("Frontend", "Legacy", "Local BA", "Loop + PGO", "Full backend")
COLORS = ("#4C78A8", "#E45756", "#72B7B2", "#F2CF5B", "#7A5195")
SEQUENCES = ("Grass02", "Outdoor04", "Z-Rough-Road01", "Dark01")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def style() -> None:
    plt.rcParams.update({
        "font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9,
        "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 160,
    })


def main() -> None:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis-dir", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "analysis" / "m3dgr_backend_controlled_20260721",
    )
    parser.add_argument("--output-dir", type=Path, default=repo / "doc" / "assets")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    style()

    summary = read_csv(args.analysis_dir / "summary_metrics.csv")
    effects = read_csv(args.analysis_dir / "paired_effects.csv")
    activity = read_csv(args.analysis_dir / "backend_activity.csv")
    lookup = {(row["sequence"], row["method"]): row for row in summary}
    outputs: list[Path] = []

    x = np.arange(len(SEQUENCES))
    width = 0.16
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.8), constrained_layout=True,
                             gridspec_kw={"width_ratios": (0.9, 1.35)})

    # Keep the failed legacy backend visible without allowing its metre-scale
    # errors to flatten the controlled new-backend comparison on the right.
    legacy_methods = ("frontend_only", "legacy_controlled")
    legacy_labels = ("Frontend", "Legacy")
    legacy_colors = (COLORS[0], COLORS[1])
    legacy_width = 0.34
    for index, (method, label, color) in enumerate(zip(legacy_methods, legacy_labels, legacy_colors, strict=True)):
        means = [float(lookup[(sequence, method)]["ate_rmse_m_mean"]) for sequence in SEQUENCES]
        stds = [float(lookup[(sequence, method)]["ate_rmse_m_std"]) for sequence in SEQUENCES]
        bars = axes[0].bar(x + (index - 0.5) * legacy_width, means, legacy_width,
                           yerr=stds, capsize=2, label=label, color=color)
        axes[0].bar_label(bars, labels=[f"{value:.3f}" for value in means],
                          padding=2, fontsize=7, rotation=90)
    axes[0].set_yscale("log")
    axes[0].set_ylim(0.12, 8.0)
    axes[0].set_title("Frontend vs legacy (log scale)")
    axes[0].set_ylabel("Translation ATE RMSE (m), lower is better")
    axes[0].legend(ncol=2, loc="upper left")

    modern_methods = ("frontend_only", "local_ba_only", "loop_pgo", "full_backend")
    modern_labels = ("Frontend", "Local BA", "Loop + PGO", "Full backend")
    modern_colors = (COLORS[0], COLORS[2], COLORS[3], COLORS[4])
    modern_width = 0.19
    for index, (method, label, color) in enumerate(zip(modern_methods, modern_labels, modern_colors, strict=True)):
        means = [float(lookup[(sequence, method)]["ate_rmse_m_mean"]) for sequence in SEQUENCES]
        stds = [float(lookup[(sequence, method)]["ate_rmse_m_std"]) for sequence in SEQUENCES]
        axes[1].bar(x + (index - 1.5) * modern_width, means, modern_width,
                    yerr=stds, capsize=2, label=label, color=color)
    axes[1].set_ylim(0.0, 0.68)
    axes[1].set_title("Controlled new-backend ablation (linear scale)")
    axes[1].set_ylabel("Translation ATE RMSE (m), lower is better")
    axes[1].legend(ncol=2, loc="upper left")
    for ax in axes:
        ax.set_xticks(x, SEQUENCES, rotation=10)
    fig.suptitle("Backend accuracy with an identical raw LIO frontend in every method", fontsize=12)
    path = args.output_dir / "m3dgr_backend_controlled_precision.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.5), constrained_layout=True)
    for index, (method, label, color) in enumerate(zip(METHODS, LABELS, COLORS, strict=True)):
        wall = [float(lookup[(sequence, method)]["wall_ms_per_output_frame_mean"]) for sequence in SEQUENCES]
        wall_std = [float(lookup[(sequence, method)]["wall_ms_per_output_frame_std"]) for sequence in SEQUENCES]
        rss = [float(lookup[(sequence, method)]["peak_rss_mb_mean"]) for sequence in SEQUENCES]
        axes[0].bar(x + (index - 2) * width, wall, width, yerr=wall_std, capsize=2, color=color, label=label)
        axes[1].bar(x + (index - 2) * width, rss, width, color=color, label=label)
    axes[0].set_title("Unpaced end-to-end processing cost")
    axes[0].set_ylabel("Wall time / output frame (ms)")
    axes[1].set_title("Peak resident memory")
    axes[1].set_ylabel("Peak RSS (MB)")
    for ax in axes:
        ax.set_xticks(x, SEQUENCES, rotation=10)
        ax.set_ylim(bottom=0)
    axes[0].legend(ncol=3, loc="upper left")
    path = args.output_dir / "m3dgr_backend_controlled_resources.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    effect_lookup = {(row["sequence"], row["comparison"]): row for row in effects}
    stages = (
        ("local_ba_minus_frontend", "Local BA", "#72B7B2"),
        ("loop_pgo_minus_local_ba", "Loop + PGO", "#F2CF5B"),
        ("hba_minus_loop_pgo", "HBA", "#7A5195"),
    )
    fig, ax = plt.subplots(figsize=(10.2, 4.7), constrained_layout=True)
    stage_width = 0.24
    for index, (comparison, label, color) in enumerate(stages):
        values = [float(effect_lookup[(sequence, comparison)]["ate_rmse_delta_percent"]) for sequence in SEQUENCES]
        bars = ax.bar(x + (index - 1) * stage_width, values, stage_width, color=color, label=label)
        for bar, value in zip(bars, values, strict=True):
            value_label = f"{value:+.3f}%" if abs(value) < 0.05 else f"{value:+.1f}%"
            ax.text(bar.get_x() + bar.get_width() / 2.0, value + (0.7 if value >= 0.0 else -0.7),
                    value_label, ha="center", va="bottom" if value >= 0.0 else "top", fontsize=7)
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xticks(x, SEQUENCES)
    ax.set_ylabel("Incremental ATE change (%)\nnegative = improvement")
    ax.set_title("Where backend accuracy changes occur (controlled incremental ablation)")
    ax.legend(ncol=3, loc="upper left")
    path = args.output_dir / "m3dgr_backend_controlled_ablation.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    full_activity = [row for row in activity if row["method"] == "full_backend"]
    means = []
    for sequence in SEQUENCES:
        rows = [row for row in full_activity if row["sequence"] == sequence]
        means.append(np.mean([float(row["loops_graph_inliers"]) for row in rows]))
    full_effect = [float(effect_lookup[(sequence, "full_minus_frontend")]["ate_rmse_delta_percent"]) for sequence in SEQUENCES]
    fig, ax = plt.subplots(figsize=(8.6, 4.6), constrained_layout=True)
    points = ax.scatter(means, full_effect, s=85, c=COLORS[-1], edgecolors="white", linewidths=0.8)
    del points
    for loop_count, effect, sequence in zip(means, full_effect, SEQUENCES, strict=True):
        ax.annotate(sequence, (loop_count, effect), xytext=(5, 5), textcoords="offset points")
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xlabel("Applied pose-graph loop inliers (mean count)")
    ax.set_ylabel("Full backend ATE change vs frontend (%)\nnegative = improvement")
    ax.set_title("Loop execution is proven, but loop count alone does not guarantee accuracy gain")
    path = args.output_dir / "m3dgr_backend_controlled_loop_effect.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    manifest = {
        "schema_version": 1,
        "generator": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
        "inputs": {
            path.name: {"path": str(path.resolve()), "sha256": sha256(path)}
            for path in (args.analysis_dir / "summary_metrics.csv", args.analysis_dir / "paired_effects.csv",
                         args.analysis_dir / "backend_activity.csv")
        },
        "outputs": {path.name: {"path": str(path.resolve()), "sha256": sha256(path)} for path in outputs},
    }
    manifest_path = args.analysis_dir / "figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
