#!/usr/bin/env python3
"""Build reproducible tables and figures for the SANY multi-LiDAR resource report."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager


CONFIG_LABELS = {
    "C0_single114_noise": "C0 单雷达",
    "C1_four_noise": "C1 四雷达主方案",
    "C2_four_no_noise": "C2 四雷达无噪声模型",
    "C3_drop114_lidar": "C3 缺前雷达",
    "C4_drop127": "C4 缺右雷达",
    "C5_drop187": "C5 缺左雷达",
    "C6_drop195": "C6 缺后雷达",
}

STAGE_LABELS = {
    "Preprocess (Standard)": "点云预处理",
    "Undistort Pcl": "点云去畸变",
    "ObsModel (Lidar Match)": "激光观测匹配",
    "Incremental Mapping": "增量建图",
    "IVox Add Points": "IVox 加点",
}

TIME_PATTERNS = {
    "user_s": re.compile(r"User time \(seconds\):\s*([0-9.]+)"),
    "system_s": re.compile(r"System time \(seconds\):\s*([0-9.]+)"),
    "cpu_percent": re.compile(r"Percent of CPU this job got:\s*([0-9.]+)%"),
    "peak_rss_kib": re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)"),
}

TIMER_PATTERN = re.compile(
    r"> \[\s*(?P<stage>[^]]+?)\s*\] average time usage: "
    r"(?P<average>[0-9.eE+-]+) ms, med: (?P<median>[0-9.eE+-]+) "
    r"95%: (?P<p95>[0-9.eE+-]+), called times\s*:\s*(?P<count>\d+)"
)


def parse_elapsed(value: str) -> float:
    parts = [float(part) for part in value.split(":")]
    if len(parts) == 2:
        return parts[0] * 60.0 + parts[1]
    if len(parts) == 3:
        return parts[0] * 3600.0 + parts[1] * 60.0 + parts[2]
    raise ValueError(f"Unsupported elapsed time: {value}")


def parse_gnu_time(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    row: dict[str, float] = {}
    for key, pattern in TIME_PATTERNS.items():
        match = pattern.search(text)
        if not match:
            raise ValueError(f"Missing {key} in {path}")
        row[key] = float(match.group(1))
    elapsed_match = re.search(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9:.]+)", text
    )
    if not elapsed_match:
        raise ValueError(f"Missing elapsed time in {path}")
    row["wall_time_s"] = parse_elapsed(elapsed_match.group(1))
    row["average_cpu_cores"] = row["cpu_percent"] / 100.0
    row["cpu_pct_of_8_core_allocation"] = row["cpu_percent"] / 8.0
    row["peak_rss_mib"] = row["peak_rss_kib"] / 1024.0
    return row


def parse_timer(path: Path) -> list[dict[str, float | int | str]]:
    rows = []
    text = path.read_text(encoding="utf-8", errors="replace")
    for match in TIMER_PATTERN.finditer(text):
        stage = match.group("stage").strip()
        if stage not in STAGE_LABELS:
            continue
        rows.append(
            {
                "stage": stage,
                "average_ms": float(match.group("average")),
                "median_ms": float(match.group("median")),
                "p95_ms": float(match.group("p95")),
                "retained_samples": int(match.group("count")),
            }
        )
    if len(rows) != len(STAGE_LABELS):
        raise ValueError(f"Expected {len(STAGE_LABELS)} timer stages in {path}, got {len(rows)}")
    return rows


def collect_formal_runs(formal_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    resource_rows = []
    stage_rows = []
    for time_path in sorted(formal_root.glob("runs/*/repeat_*/attempt_*/gnu_time.txt")):
        run_dir = time_path.parent
        config = run_dir.parents[1].name
        repeat = run_dir.parents[0].name
        resource = parse_gnu_time(time_path)
        resource_rows.append(
            {
                "config": config,
                "config_label": CONFIG_LABELS[config],
                "repeat": repeat,
                **resource,
            }
        )
        timer_path = run_dir / "logs" / "run_frontend_offline.stderr.log"
        for timing in parse_timer(timer_path):
            stage_rows.append(
                {
                    "config": config,
                    "config_label": CONFIG_LABELS[config],
                    "repeat": repeat,
                    "stage_label": STAGE_LABELS[str(timing["stage"])],
                    **timing,
                }
            )
    resources = pd.DataFrame(resource_rows)
    stages = pd.DataFrame(stage_rows)
    if len(resources) != 21:
        raise ValueError(f"Expected 21 formal runs, got {len(resources)}")
    if len(stages) != 21 * len(STAGE_LABELS):
        raise ValueError(f"Expected {21 * len(STAGE_LABELS)} stage rows, got {len(stages)}")
    return resources, stages


def configure_style() -> None:
    candidates = [
        Path("/mnt/c/Windows/Fonts/msyh.ttc"),
        Path("C:/Windows/Fonts/msyh.ttc"),
    ]
    for candidate in candidates:
        if candidate.exists():
            font_manager.fontManager.addfont(str(candidate))
            family = font_manager.FontProperties(fname=str(candidate)).get_name()
            plt.rcParams["font.family"] = family
            break
    plt.rcParams.update(
        {
            "axes.unicode_minus": False,
            "figure.facecolor": "white",
            "axes.facecolor": "#FAFBFC",
            "axes.edgecolor": "#AAB2BD",
            "axes.labelcolor": "#263238",
            "text.color": "#263238",
            "xtick.color": "#4B5563",
            "ytick.color": "#4B5563",
            "grid.color": "#D9DEE5",
            "grid.linewidth": 0.8,
        }
    )


def summarize_resources(resources: pd.DataFrame) -> pd.DataFrame:
    metrics = ["wall_time_s", "peak_rss_mib", "average_cpu_cores", "cpu_pct_of_8_core_allocation"]
    summary = resources.groupby(["config", "config_label"], sort=False)[metrics].agg(
        ["median", "min", "max"]
    )
    summary.columns = [f"{metric}_{stat}" for metric, stat in summary.columns]
    return summary.reset_index()


def plot_resource_overview(summary: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(17, 7.5), constrained_layout=True)
    labels = summary["config_label"].tolist()
    y = np.arange(len(labels))
    specs = [
        ("wall_time_s", "墙钟时间", "秒"),
        ("peak_rss_mib", "峰值常驻内存（RSS）", "MiB"),
        ("average_cpu_cores", "平均 CPU 占用", "逻辑核"),
    ]
    colors = ["#AAB7C4" if config == "C0_single114_noise" else "#3568A8" for config in summary["config"]]
    colors[summary.index[summary["config"] == "C1_four_noise"][0]] = "#D58A2B"
    for axis, (metric, title, unit) in zip(axes, specs):
        median = summary[f"{metric}_median"].to_numpy()
        lower = median - summary[f"{metric}_min"].to_numpy()
        upper = summary[f"{metric}_max"].to_numpy() - median
        axis.barh(y, median, color=colors, edgecolor="#263238", linewidth=0.6)
        axis.errorbar(median, y, xerr=np.vstack([lower, upper]), fmt="none", ecolor="#263238", capsize=3)
        axis.set_title(title, fontsize=14, weight="bold")
        axis.set_xlabel(unit)
        axis.set_yticks(y, labels if axis is axes[0] else [""] * len(labels))
        axis.invert_yaxis()
        axis.grid(axis="x")
        axis.set_axisbelow(True)
        axis.set_xlim(0, max(summary[f"{metric}_max"]) * 1.18)
        for yi, value in zip(y, median):
            label = f"{value:.2f}" if metric != "peak_rss_mib" else f"{value:,.0f}"
            axis.text(value, yi, f"  {label}", va="center", fontsize=9)
    fig.suptitle(
        "SANY 七种配置的离线资源开销\n"
        "每种配置 3 次完整运行；柱为中位数，误差线为最小值—最大值；CPU 固定在逻辑核 0–7",
        fontsize=18,
        weight="bold",
    )
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_resource_timeseries(samples: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True, constrained_layout=True)
    x = samples["elapsed_s"]
    series = [
        ("cpu_cores", "进程树 CPU 占用", "逻辑核", "#3568A8"),
        ("rss_mb", "进程树常驻内存（RSS）", "MiB", "#D58A2B"),
    ]
    for axis, (column, title, unit, color) in zip(axes, series):
        values = samples[column]
        axis.plot(x, values, color=color, linewidth=1.2)
        axis.fill_between(x, values, color=color, alpha=0.14)
        mean = values.mean()
        p95 = values.quantile(0.95)
        peak = values.max()
        axis.axhline(mean, color="#263238", linestyle="--", linewidth=1, label=f"均值 {mean:.2f}")
        axis.axhline(p95, color="#8B5A2B", linestyle=":", linewidth=1.2, label=f"P95 {p95:.2f}")
        axis.set_title(title, fontsize=14, weight="bold")
        axis.set_ylabel(unit)
        axis.set_ylim(bottom=0)
        axis.grid(axis="y")
        axis.legend(loc="upper left", frameon=False, ncol=2)
        axis.text(x.max(), peak, f" 峰值 {peak:.2f}", va="bottom", ha="right", fontsize=9)
    axes[-1].set_xlabel("运行经过时间（秒）")
    fig.suptitle(
        "C1 四雷达主方案完整运行的资源时序\n"
        "采样间隔约 0.2 s，共 656 点；统计对象为统一离线脚本启动的完整进程树",
        fontsize=18,
        weight="bold",
    )
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_stage_timing(stages: pd.DataFrame, output: Path) -> None:
    selected = stages[stages["config"].isin(["C0_single114_noise", "C1_four_noise", "C2_four_no_noise"])]
    medians = selected.groupby(["config", "stage"], sort=False)["average_ms"].median().reset_index()
    stage_order = list(STAGE_LABELS)
    config_order = ["C0_single114_noise", "C1_four_noise", "C2_four_no_noise"]
    colors = ["#AAB7C4", "#D58A2B", "#3568A8"]
    x = np.arange(len(stage_order))
    width = 0.24
    fig, axis = plt.subplots(figsize=(14, 7), constrained_layout=True)
    for index, (config, color) in enumerate(zip(config_order, colors)):
        values = []
        for stage in stage_order:
            value = medians[(medians["config"] == config) & (medians["stage"] == stage)]["average_ms"].iloc[0]
            values.append(value)
        positions = x + (index - 1) * width
        bars = axis.bar(
            positions,
            values,
            width,
            label=CONFIG_LABELS[config],
            color=color,
            edgecolor="#263238",
            linewidth=0.6,
        )
        axis.bar_label(bars, fmt="%.2f", padding=2, fontsize=8, rotation=90)
    axis.set_title(
        "核心算法阶段的单次调用平均耗时\n"
        "每根柱为 3 次完整运行中日志所报平均耗时的中位数；阶段可能嵌套，禁止相加推导墙钟时间",
        fontsize=16,
        weight="bold",
    )
    axis.set_ylabel("毫秒 / 调用")
    axis.set_xticks(x, [STAGE_LABELS[stage] for stage in stage_order])
    axis.set_ylim(0, medians["average_ms"].max() * 1.28)
    axis.grid(axis="y")
    axis.set_axisbelow(True)
    axis.legend(frameon=False, ncol=3, loc="upper left")
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_source_notes(output_dir: Path, formal_root: Path, manual_run: Path) -> None:
    text = f"""# 资源报告数据与图表说明

## 数据源

- 正式矩阵：`{formal_root}`
- 四雷达完整运行资源采样：`{manual_run / 'resource_samples.csv'}`
- 四雷达完整运行资源摘要：`{manual_run / 'resource_summary.json'}`

## 指标口径

- 墙钟时间：GNU time 的 `Elapsed (wall clock) time`，包含算法、ROS 2 bag 读取、结果写盘和脚本开销。
- 平均 CPU 逻辑核：GNU time 的 CPU 百分比除以 100；占 8 核配额比例再除以 8。
- 峰值 RSS：GNU time 的 `Maximum resident set size`，由 KiB 换算为 MiB。
- 时序 CPU/RSS：统一离线脚本对子进程树约每 0.2 秒采样一次。
- 阶段耗时：算法 Timer 的最近样本统计；不同阶段可能嵌套，不能相加为端到端耗时。

## 图表地图

| 报告段落 | 分析问题 | 图表 | 字段 | 支持的结论 |
|---|---|---|---|---|
| 总体资源 | 七种配置的资源代价有多大 | 横向分组比较 | 墙钟时间、峰值 RSS、平均 CPU 核 | 四雷达主要增加内存与总耗时，CPU 平均并未占满 8 核 |
| 运行时序 | C1 的资源峰值在何时出现 | 双面板时序折线 | elapsed_s、cpu_cores、rss_mb | 内存随地图增长累积，CPU 存在阶段性波动 |
| 阶段耗时 | 核心计算热点在哪里 | 分组柱状图 | average_ms、config、stage | 激光观测匹配是已埋点阶段中最重的单次调用 |

## 限制

- 正式矩阵的三次重复为确定性技术重复，不等同于三段独立数据或三台机器。
- 正式矩阵采自历史冻结安装树；当前统一分支增加了耗时导出，但未重新执行全部 21 个完整单元。
- WSL 2 + Windows 挂载盘的 I/O 与调度不代表目标车载原生 Linux 平台。
- GNU time 的 CPU 百分比是整个任务平均值，无法定位线程级并行效率。
"""
    (output_dir / "source_notes.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--formal-root", type=Path, required=True)
    parser.add_argument("--manual-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    configure_style()

    resources, stages = collect_formal_runs(args.formal_root)
    summary = summarize_resources(resources)
    samples = pd.read_csv(args.manual_run / "resource_samples.csv")

    resources.to_csv(args.output_dir / "formal_resource_runs.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    summary.to_csv(args.output_dir / "formal_resource_summary.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    stages.to_csv(args.output_dir / "formal_stage_timing.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    plot_resource_overview(summary, args.output_dir / "formal_resource_overview.png")
    plot_resource_timeseries(samples, args.output_dir / "c1_resource_timeseries.png")
    plot_stage_timing(stages, args.output_dir / "core_stage_timing.png")
    write_source_notes(args.output_dir, args.formal_root, args.manual_run)


if __name__ == "__main__":
    main()
