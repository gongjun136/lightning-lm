#!/usr/bin/env python3
"""Build deterministic formal-report figures from validated matrix outputs."""

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


COLORS = {
    "lightning_lm": "#0072B2",
    "lightning_frontend": "#56B4E9",
    "fastlio": "#E69F00",
    "fastlivo2_lio": "#009E73",
    "voxel_slam_frontend": "#CC79A7",
    "lightning_legacy": "#D55E00",
    "lightning_new_backend": "#0072B2",
    "voxel_slam_full": "#CC79A7",
    "startup": "#0072B2",
    "forced_loss": "#D55E00",
}
DISPLAY = {
    "lightning_lm": "Lightning-LM",
    "lightning_frontend": "Lightning frontend",
    "fastlio": "FAST-LIO",
    "fastlivo2_lio": "FAST-LIVO2 LIO",
    "voxel_slam_frontend": "Voxel-SLAM LIO",
    "lightning_legacy": "Lightning legacy",
    "lightning_new_backend": "Lightning new",
    "voxel_slam_full": "Voxel-SLAM full",
    "startup": "Startup",
    "forced_loss": "Forced loss",
}


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    formal = repo / "runs" / "formal_report_20260717"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, default=formal / "analysis")
    parser.add_argument(
        "--output-dir", type=Path,
        default=repo / "doc" / "assets" / "formal_report_20260717",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in {"", "None", None} else np.nan


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


def grouped_summary(
    rows: list[dict[str, str]],
    group_key: str,
    groups: tuple[str, ...],
    method_key: str,
    methods: tuple[str, ...],
    output: Path,
    title: str,
) -> None:
    metrics = (
        ("ate_rmse_m_mean", "ate_rmse_m_std", "ATE RMSE (m)"),
        ("mean_cpu_cores_mean", "mean_cpu_cores_std", "Mean CPU use (cores)"),
        ("peak_rss_mb_mean", "peak_rss_mb_std", "Peak RSS (MB)"),
    )
    lookup = {(row[group_key], row[method_key]): row for row in rows}
    # Keep enough horizontal room for the long sequence labels on every panel.
    # A narrower canvas lets the first tick label of one panel intrude into the
    # preceding panel after ``bbox_inches='tight'`` is applied at save time.
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.4))
    x = np.arange(len(groups), dtype=float)
    width = 0.8 / len(methods)
    for axis, (mean_key, std_key, ylabel) in zip(axes, metrics, strict=True):
        for index, method in enumerate(methods):
            selected = [lookup[(group, method)] for group in groups]
            means = np.asarray([number(row, mean_key) for row in selected])
            stds = np.asarray([number(row, std_key) for row in selected])
            offset = (index - (len(methods) - 1) / 2) * width
            axis.bar(
                x + offset, means, width=width * 0.92, yerr=stds,
                color=COLORS[method], label=DISPLAY[method], capsize=2,
            )
        axis.set_xticks(x, groups, rotation=20, ha="right")
        axis.tick_params(axis="x", labelsize=8)
        axis.set_ylabel(ylabel)
        axis.set_ylim(bottom=0)
    axes[0].legend(frameon=False, ncol=2)
    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.text(0.5, 0.005, "Bars: mean of 3 runs; error bars: sample standard deviation.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95), w_pad=3.0)
    save(fig, output)


def trajectory_error_grid(
    rows: list[dict[str, str]],
    group_key: str,
    groups: tuple[str, ...],
    method_key: str,
    methods: tuple[str, ...],
    output: Path,
    title: str,
) -> None:
    fig, axes = plt.subplots(len(groups), 2, figsize=(12.5, 3.25 * len(groups)), squeeze=False)
    for row_index, group in enumerate(groups):
        group_rows = [row for row in rows if row[group_key] == group and int(row["repeat"]) == 1]
        truth_rows = [row for row in group_rows if row[method_key] == methods[0]]
        axes[row_index, 0].plot(
            [number(row, "truth_x_m") for row in truth_rows],
            [number(row, "truth_y_m") for row in truth_rows],
            color="#202020", linewidth=1.6, label="Reference",
        )
        for method in methods:
            selected = [row for row in group_rows if row[method_key] == method]
            times = np.asarray([number(row, "timestamp_s") for row in selected])
            axes[row_index, 0].plot(
                [number(row, "estimate_x_m") for row in selected],
                [number(row, "estimate_y_m") for row in selected],
                color=COLORS[method], linewidth=1.0, label=DISPLAY[method],
            )
            axes[row_index, 1].plot(
                times - times[0], [number(row, "translation_error_m") for row in selected],
                color=COLORS[method], linewidth=0.9, label=DISPLAY[method],
            )
        axes[row_index, 0].set_title(f"{group}: aligned XY trajectory")
        axes[row_index, 0].set_xlabel("X (m)")
        axes[row_index, 0].set_ylabel("Y (m)")
        axes[row_index, 0].axis("equal")
        axes[row_index, 1].set_title(f"{group}: translation error over time")
        axes[row_index, 1].set_xlabel("Elapsed reference time (s)")
        axes[row_index, 1].set_ylabel("Translation error (m)")
        axes[row_index, 1].set_ylim(bottom=0)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.997)
    fig.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.957),
        ncol=min(5, len(labels)), frameon=False,
    )
    fig.text(0.5, 0.002, "Predeclared representative run: repeat 1; alignment uses fixed-scale SE(3).", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.02, 1, 0.91))
    save(fig, output)


def sany_mapping_summary(rows: list[dict[str, str]], output: Path) -> None:
    order = (
        ("single114", "frontend", "Single FE"),
        ("four_lidar", "frontend", "Four FE"),
        ("four_no_noise", "frontend", "Four no-noise"),
        ("drop114", "frontend", "Drop 114"),
        ("drop127", "frontend", "Drop 127"),
        ("drop187", "frontend", "Drop 187"),
        ("drop195", "frontend", "Drop 195"),
        ("single114", "slam", "Single SLAM"),
        ("four_lidar", "slam", "Four SLAM"),
    )
    lookup = {(row["variant"], row["pipeline"]): row for row in rows}
    colors = ["#56B4E9" if pipeline == "frontend" else "#0072B2" for _, pipeline, _ in order]
    metrics = (
        ("ate_rmse_m_mean", "ate_rmse_m_std", "ATE RMSE (m)"),
        ("mean_cpu_cores_mean", "mean_cpu_cores_std", "Mean CPU use (cores)"),
        ("peak_rss_mb_mean", "peak_rss_mb_std", "Peak RSS (MB)"),
        ("partial_ratio_mean", "partial_ratio_std", "Partial assembled-frame ratio"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(14, 8.2))
    x = np.arange(len(order))
    labels = [item[2] for item in order]
    selected = [lookup[(variant, pipeline)] for variant, pipeline, _ in order]
    for metric_index, (axis, (mean_key, std_key, ylabel)) in enumerate(
        zip(axes.flat, metrics, strict=True)
    ):
        means = [number(row, mean_key) for row in selected]
        stds = [number(row, std_key) for row in selected]
        axis.bar(x, means, yerr=stds, color=colors, capsize=2)
        axis.set_xticks(x, labels, rotation=28, ha="right")
        axis.set_ylabel(ylabel)
        if metric_index == 0:
            # The no-noise ablation fails by four orders of magnitude relative
            # to the other cells. A logarithmic axis keeps the failure visible
            # without visually collapsing all successful configurations to 0.
            positive = [value for value in means if np.isfinite(value) and value > 0]
            axis.set_yscale("log")
            axis.set_ylabel("ATE RMSE (m, log scale)")
            axis.set_ylim(min(positive) / 2.0, max(positive) * 2.0)
        else:
            axis.set_ylim(bottom=0)
    fig.suptitle("SANY 20260701: multi-LiDAR accuracy, efficiency and robustness", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.005, "Bars: mean of 3 runs; error bars: sample standard deviation. FE = frontend only.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    save(fig, output)


def sany_mapping_trajectory(rows: list[dict[str, str]], output: Path) -> None:
    cells = (
        ("single114", "frontend", "Single LiDAR frontend", "#56B4E9"),
        ("four_lidar", "frontend", "Four LiDAR frontend", "#009E73"),
        ("single114", "slam", "Single LiDAR SLAM", "#D55E00"),
        ("four_lidar", "slam", "Four LiDAR SLAM", "#0072B2"),
    )
    selected_all = [row for row in rows if int(row["repeat"]) == 1]
    truth = [row for row in selected_all if row["variant"] == "single114" and row["pipeline"] == "frontend"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    axes[0].plot(
        [number(row, "truth_x_m") for row in truth], [number(row, "truth_y_m") for row in truth],
        color="#202020", linewidth=1.6, label="Voxel-SLAM 114 proxy",
    )
    for variant, pipeline, label, color in cells:
        selected = [row for row in selected_all if row["variant"] == variant and row["pipeline"] == pipeline]
        times = np.asarray([number(row, "timestamp_s") for row in selected])
        axes[0].plot(
            [number(row, "estimate_x_m") for row in selected],
            [number(row, "estimate_y_m") for row in selected], color=color, linewidth=1.0, label=label,
        )
        axes[1].plot(
            times - times[0], [number(row, "translation_error_m") for row in selected],
            color=color, linewidth=0.9, label=label,
        )
    axes[0].set(title="Aligned XY trajectory", xlabel="X (m)", ylabel="Y (m)")
    axes[0].axis("equal")
    axes[1].set(title="Translation error over time", xlabel="Elapsed proxy-reference time (s)", ylabel="Translation error (m)")
    axes[1].set_ylim(bottom=0)
    axes[0].legend(frameon=False, fontsize=7)
    fig.suptitle("SANY 20260701: single/four-LiDAR mapping against the Voxel-SLAM proxy", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.005, "Predeclared representative run: repeat 1; proxy reference is not absolute ground truth.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    save(fig, output)


def sany_relocalization_summary(rows: list[dict[str, str]], output: Path) -> None:
    order = tuple((dataset, mode) for dataset in ("data1", "data2") for mode in ("startup", "forced_loss"))
    lookup = {(row["dataset"], row["mode"]): row for row in rows}
    selected = [lookup[key] for key in order]
    labels = [f"{dataset}\n{DISPLAY[mode]}" for dataset, mode in order]
    colors = [COLORS[mode] for _, mode in order]
    metrics = (
        ("ate_rmse_m_mean", "ate_rmse_m_std", "ATE RMSE (m)"),
        ("startup_latency_s_mean", "startup_latency_s_std", "Initial acceptance latency (s)"),
        ("recovery_latency_s_mean", "recovery_latency_s_std", "Forced-loss recovery latency (s)"),
        ("peak_rss_mb_mean", "peak_rss_mb_std", "Peak RSS (MB)"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8))
    x = np.arange(len(order))
    for axis, (mean_key, std_key, ylabel) in zip(axes.flat, metrics, strict=True):
        means = np.asarray([number(row, mean_key) for row in selected])
        stds = np.asarray([number(row, std_key) for row in selected])
        mask = np.isfinite(means)
        axis.bar(x[mask], means[mask], yerr=stds[mask], color=np.asarray(colors)[mask], capsize=2)
        axis.set_xticks(x, labels)
        axis.set_ylabel(ylabel)
        axis.set_ylim(bottom=0)
    fig.suptitle("SANY 20260716: startup and forced-loss relocalization", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.005, "Bars: mean of 3 runs; error bars: sample standard deviation.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    save(fig, output)


def sany_relocalization_trajectory(rows: list[dict[str, str]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2))
    for index, dataset in enumerate(("data1", "data2")):
        dataset_rows = [row for row in rows if row["dataset"] == dataset and int(row["repeat"]) == 1]
        truth = [row for row in dataset_rows if row["mode"] == "startup"]
        axes[index, 0].plot(
            [number(row, "truth_x_m") for row in truth], [number(row, "truth_y_m") for row in truth],
            color="#202020", linewidth=1.6, label="Voxel-SLAM 114 proxy",
        )
        for mode in ("startup", "forced_loss"):
            selected = [row for row in dataset_rows if row["mode"] == mode]
            times = np.asarray([number(row, "timestamp_s") for row in selected])
            axes[index, 0].plot(
                [number(row, "estimate_x_m") for row in selected],
                [number(row, "estimate_y_m") for row in selected],
                color=COLORS[mode], linewidth=1.0, label=DISPLAY[mode],
            )
            axes[index, 1].plot(
                times - times[0], [number(row, "translation_error_m") for row in selected],
                color=COLORS[mode], linewidth=0.9, label=DISPLAY[mode],
            )
        axes[index, 0].set(title=f"{dataset}: aligned XY trajectory", xlabel="X (m)", ylabel="Y (m)")
        axes[index, 0].axis("equal")
        axes[index, 1].set(
            title=f"{dataset}: translation error over time",
            xlabel="Elapsed proxy-reference time (s)", ylabel="Translation error (m)",
        )
        axes[index, 1].set_ylim(bottom=0)
    axes[0, 0].legend(frameon=False)
    fig.suptitle("SANY 20260716: localization trajectories and error evolution", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.005, "Predeclared representative run: repeat 1; proxy reference is not absolute ground truth.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    save(fig, output)


def main() -> int:
    args = parse_args()
    inputs = {
        f"{name}/{filename}": args.analysis_root / name / filename
        for name in (
            "m3dgr_frontend", "m3dgr_backend", "m3dgr_localization",
            "sany_20260701_mapping", "sany_20260716_relocalization",
        )
        for filename in ("summary_metrics.csv", "aligned_samples.csv", "validation.json")
    }
    missing = [str(path) for path in inputs.values() if not path.is_file()]
    if missing:
        raise SystemExit("missing validated analysis artifacts:\n" + "\n".join(missing))
    failed = []
    for name in (
        "m3dgr_frontend", "m3dgr_backend", "m3dgr_localization",
        "sany_20260701_mapping", "sany_20260716_relocalization",
    ):
        validation = json.loads((args.analysis_root / name / "validation.json").read_text(encoding="utf-8"))
        if validation["status"] not in {"passed", "passed_with_warnings"}:
            failed.append(f"{name}: {validation['status']}")
    if failed:
        raise SystemExit("refusing to chart failed analyses:\n" + "\n".join(failed))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    configure_style()
    sequences = ("Grass02", "Outdoor04", "Z-Rough-Road01", "Dark01")

    frontend_summary = read_csv(args.analysis_root / "m3dgr_frontend" / "summary_metrics.csv")
    frontend_samples = read_csv(args.analysis_root / "m3dgr_frontend" / "aligned_samples.csv")
    frontend_methods = ("lightning_lm", "fastlio", "fastlivo2_lio", "voxel_slam_frontend")
    grouped_summary(
        frontend_summary, "sequence", sequences, "method", frontend_methods,
        args.output_dir / "m3dgr_frontend_summary.png",
        "M3DGR frontend: accuracy and resource trade-offs across four challenge sequences",
    )
    trajectory_error_grid(
        frontend_samples, "sequence", sequences, "method", frontend_methods,
        args.output_dir / "m3dgr_frontend_trajectory_error.png",
        "M3DGR frontend: trajectory and translation error",
    )
    del frontend_summary, frontend_samples

    backend_summary = read_csv(args.analysis_root / "m3dgr_backend" / "summary_metrics.csv")
    backend_samples = read_csv(args.analysis_root / "m3dgr_backend" / "aligned_samples.csv")
    backend_methods = (
        "lightning_frontend", "lightning_legacy", "lightning_new_backend", "voxel_slam_full"
    )
    grouped_summary(
        backend_summary, "sequence", sequences, "method", backend_methods,
        args.output_dir / "m3dgr_backend_summary.png",
        "M3DGR backend effect: frontend, legacy, new and Voxel-SLAM full comparison",
    )
    trajectory_error_grid(
        backend_samples, "sequence", sequences, "method", backend_methods,
        args.output_dir / "m3dgr_backend_trajectory_error.png",
        "M3DGR backend effect: trajectory and translation error",
    )
    del backend_summary, backend_samples

    localization_summary = read_csv(args.analysis_root / "m3dgr_localization" / "summary_metrics.csv")
    localization_samples = read_csv(args.analysis_root / "m3dgr_localization" / "aligned_samples.csv")
    map_variants = ("lightning_new_backend", "voxel_slam_full")
    grouped_summary(
        localization_summary, "sequence", sequences, "map_variant", map_variants,
        args.output_dir / "m3dgr_localization_summary.png",
        "M3DGR localization: dual map-source comparison",
    )
    trajectory_error_grid(
        localization_samples, "sequence", sequences, "map_variant", map_variants,
        args.output_dir / "m3dgr_localization_trajectory_error.png",
        "M3DGR localization: trajectory and translation error",
    )
    del localization_summary, localization_samples

    mapping_summary = read_csv(args.analysis_root / "sany_20260701_mapping" / "summary_metrics.csv")
    mapping_samples = read_csv(args.analysis_root / "sany_20260701_mapping" / "aligned_samples.csv")
    sany_mapping_summary(mapping_summary, args.output_dir / "sany_mapping_summary.png")
    sany_mapping_trajectory(mapping_samples, args.output_dir / "sany_mapping_trajectory_error.png")
    del mapping_summary, mapping_samples

    relocalization_summary = read_csv(args.analysis_root / "sany_20260716_relocalization" / "summary_metrics.csv")
    relocalization_samples = read_csv(args.analysis_root / "sany_20260716_relocalization" / "aligned_samples.csv")
    sany_relocalization_summary(relocalization_summary, args.output_dir / "sany_relocalization_summary.png")
    sany_relocalization_trajectory(
        relocalization_samples, args.output_dir / "sany_relocalization_trajectory_error.png"
    )

    # Keep this manifest scoped to the ten figures generated above.  The four
    # frontend-compute figures share the directory but have their own manifest.
    output_names = (
        "m3dgr_frontend_summary.png",
        "m3dgr_frontend_trajectory_error.png",
        "m3dgr_backend_summary.png",
        "m3dgr_backend_trajectory_error.png",
        "m3dgr_localization_summary.png",
        "m3dgr_localization_trajectory_error.png",
        "sany_mapping_summary.png",
        "sany_mapping_trajectory_error.png",
        "sany_relocalization_summary.png",
        "sany_relocalization_trajectory_error.png",
    )
    outputs = [args.output_dir / name for name in output_names]
    manifest = {
        "schema_version": 1,
        "chart_contract": {
            "summary_statistics": "mean of 3 independent runs; sample standard deviation; worst values remain in CSV",
            "trajectory_selection": "repeat 1 selected before result inspection",
            "alignment": "fixed-scale SE(3) on each analysis-defined common support",
            "resource_panels": "mean CPU cores and peak process-group RSS; ATE is comparison context",
            "bar_baseline": "zero",
        },
        "inputs": {str(path.resolve()): sha256(path) for path in sorted(inputs.values())},
        "outputs": {path.name: sha256(path) for path in outputs},
    }
    (args.output_dir / "figure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output_dir": str(args.output_dir), "figures": len(outputs)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
