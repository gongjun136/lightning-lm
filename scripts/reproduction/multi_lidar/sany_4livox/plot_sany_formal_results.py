#!/usr/bin/env python3
"""Create publication-style SANY formal-result figures.

Conclusion encoded by the figure: four-LiDAR fusion with the heteroscedastic
measurement model improves proxy trajectory consistency and preserves broad
cross-source coverage; disabling the model causes joint trajectory/map failure.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ORDER = [
    "C0_single114_noise", "C1_four_noise", "C2_four_no_noise",
    "C3_drop114_lidar", "C4_drop127", "C5_drop187", "C6_drop195",
]
LABELS = ["C0\n单114", "C1\n四雷达", "C2\n无噪声模型", "C3\n缺114", "C4\n缺127", "C5\n缺187", "C6\n缺195"]
COLORS = ["#6B7280", "#0072B2", "#D55E00", "#56B4E9", "#009E73", "#CC79A7", "#E69F00"]


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory-json", type=Path, required=True)
    parser.add_argument("--map-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    trajectory = load(args.trajectory_json)["variants"]
    maps = load(args.map_json)["variants"]
    x = np.arange(len(ORDER))
    ape = np.asarray([trajectory[name]["ape_translation_rmse_m_median"] for name in ORDER])
    residual50 = np.asarray([
        np.nan if maps[name]["residual_median_m_median"] is None else maps[name]["residual_median_m_median"]
        for name in ORDER
    ])
    residual95 = np.asarray([
        np.nan if maps[name]["residual_p95_m_median"] is None else maps[name]["residual_p95_m_median"]
        for name in ORDER
    ])
    coverage2 = np.asarray([maps[name]["coverage_at_least_2_median"] for name in ORDER])
    coverage3 = np.asarray([maps[name]["coverage_at_least_3_median"] for name in ORDER])
    coverage4 = np.asarray([maps[name]["coverage_all_4_median"] for name in ORDER])

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.7), constrained_layout=True)

    ax = axes[0]
    ax.bar(x, ape, color=COLORS, width=0.72)
    ax.set_yscale("log")
    ax.set_ylabel("平移 APE RMSE（m，对数轴）")
    ax.set_title("a  轨迹代理一致性")
    ax.set_xticks(x, LABELS)
    ax.grid(axis="y", which="both", color="#D1D5DB", linewidth=0.6, alpha=0.8)
    for index, value in enumerate(ape):
        label = f"{value:.3f}" if value < 1 else f"{value:.1f}"
        ax.text(index, value * (1.22 if value < 1 else 0.62), label, ha="center", va="bottom", fontsize=7)

    ax = axes[1]
    width = 0.34
    ax.bar(x - width / 2, residual50, width=width, color="#0072B2", label="中位数")
    ax.bar(x + width / 2, residual95, width=width, color="#E69F00", label="P95")
    ax.set_ylabel("互惠跨源平面残差（m）")
    ax.set_title("b  地图跨源几何一致性")
    ax.set_xticks(x, LABELS)
    ax.legend(frameon=False, ncol=2, loc="upper left")
    ax.grid(axis="y", color="#D1D5DB", linewidth=0.6, alpha=0.8)
    ax.text(0, 0.01, "N/A", ha="center", va="bottom", fontsize=8, color="#4B5563")

    ax = axes[2]
    ax.plot(x, coverage2 * 100, "o-", color="#0072B2", label="≥2 源")
    ax.plot(x, coverage3 * 100, "s-", color="#009E73", label="≥3 源")
    ax.plot(x, coverage4 * 100, "^-", color="#D55E00", label="4 源")
    ax.set_ylim(-3, 103)
    ax.set_ylabel("共同覆盖率（%）")
    ax.set_title("c  多源空间覆盖")
    ax.set_xticks(x, LABELS)
    ax.legend(frameon=False, ncol=3, loc="upper center", fontsize=8)
    ax.grid(axis="y", color="#D1D5DB", linewidth=0.6, alpha=0.8)

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", labelsize=7.5)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.output_dir / "SANY_多雷达正式实验_轨迹与地图质量"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(stem)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
