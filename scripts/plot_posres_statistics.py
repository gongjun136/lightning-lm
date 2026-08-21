#!/usr/bin/env python3
"""Visualize timing and per-second receive rate produced by pos_res_recorder."""

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot /PosRes receive intervals, end-to-end latency, and message rate."
    )
    parser.add_argument("--timing", required=True, help="posres_timing.csv")
    parser.add_argument("--rate", required=True, help="posres_rate.csv")
    parser.add_argument(
        "--output", default="posres_statistics.png", help="output image (default: %(default)s)"
    )
    parser.add_argument("--dpi", type=int, default=160, help="output image DPI")
    parser.add_argument("--show", action="store_true", help="also open an interactive window")
    return parser.parse_args()


def read_numeric_csv(path: Path) -> Dict[str, List[float]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        result = {name: [] for name in reader.fieldnames}
        for row_number, row in enumerate(reader, start=2):
            try:
                for name in reader.fieldnames:
                    result[name].append(float(row[name]))
            except (TypeError, ValueError) as error:
                raise ValueError(f"invalid numeric value at {path}:{row_number}") from error
    return result


def require_columns(data: Dict[str, List[float]], path: Path, names: List[str]) -> None:
    missing = [name for name in names if name not in data]
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")


def finite_pairs(x_values: List[float], y_values: List[float]):
    return [(x, y) for x, y in zip(x_values, y_values) if math.isfinite(x) and math.isfinite(y)]


def percentile(values: List[float], probability: float) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    position = (len(finite) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return finite[lower]
    return finite[lower] * (upper - position) + finite[upper] * (position - lower)


def print_summary(timing: Dict[str, List[float]], rate: Dict[str, List[float]]) -> None:
    intervals = timing["interarrival_receive_ms"]
    latencies = timing["latency_ms"]
    counts = rate["message_count"]
    rates = rate["rate_hz"]
    print(f"messages={len(timing['sequence'])} windows={len(counts)}")
    print(
        "receive_interval_ms: "
        f"median={percentile(intervals, 0.5):.3f} "
        f"p95={percentile(intervals, 0.95):.3f} "
        f"max={percentile(intervals, 1.0):.3f}"
    )
    print(
        "latency_ms: "
        f"median={percentile(latencies, 0.5):.3f} "
        f"p95={percentile(latencies, 0.95):.3f} "
        f"max={percentile(latencies, 1.0):.3f}"
    )
    if counts:
        print(
            f"messages_per_window: min={min(counts):.0f} max={max(counts):.0f}; "
            f"rate_hz: median={percentile(rates, 0.5):.3f}"
        )


def main() -> int:
    args = parse_args()
    timing_path = Path(args.timing)
    rate_path = Path(args.rate)
    timing = read_numeric_csv(timing_path)
    rate = read_numeric_csv(rate_path)
    require_columns(
        timing,
        timing_path,
        ["sequence", "receive_elapsed_s", "latency_ms", "interarrival_receive_ms"],
    )
    require_columns(
        rate,
        rate_path,
        [
            "window_start_elapsed_s",
            "window_end_elapsed_s",
            "window_duration_s",
            "message_count",
            "rate_hz",
        ],
    )
    if not timing["sequence"]:
        raise ValueError(f"timing CSV has no samples: {timing_path}")

    elapsed = timing["receive_elapsed_s"]
    interval_points = finite_pairs(elapsed, timing["interarrival_receive_ms"])
    latency_points = finite_pairs(elapsed, timing["latency_ms"])
    rate_x = [
        (begin + end) * 0.5
        for begin, end in zip(rate["window_start_elapsed_s"], rate["window_end_elapsed_s"])
    ]

    figure, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True, constrained_layout=True)
    if interval_points:
        axes[0].plot(*zip(*interval_points), linewidth=0.8, color="#2563eb")
    axes[0].set_ylabel("Interval (ms)")
    axes[0].set_title("/PosRes receive interval")
    axes[0].grid(True, alpha=0.3)

    if latency_points:
        axes[1].plot(*zip(*latency_points), linewidth=0.8, color="#1d4ed8")
    else:
        axes[1].text(
            0.5,
            0.5,
            "No valid header timestamps",
            ha="center",
            va="center",
            transform=axes[1].transAxes,
        )
    axes[1].axhline(0.0, color="black", linewidth=0.6, alpha=0.5)
    axes[1].set_ylabel("Latency (ms)")
    axes[1].set_title("Receive time - message header timestamp")
    axes[1].grid(True, alpha=0.3)

    if rate_x:
        bar_width = [duration * 0.8 for duration in rate["window_duration_s"]]
        axes[2].bar(
            rate_x,
            rate["message_count"],
            width=bar_width,
            color="#93c5fd",
            edgecolor="#2563eb",
            linewidth=0.6,
            label="count",
        )
        rate_axis = axes[2].twinx()
        rate_axis.plot(rate_x, rate["rate_hz"], color="#1e3a8a", marker=".", label="rate")
        rate_axis.set_ylabel("Rate (Hz)", color="#1e3a8a")
        rate_axis.tick_params(axis="y", labelcolor="#1e3a8a")
    axes[2].set_xlabel("Elapsed receive time (s)")
    axes[2].set_ylabel("Messages / window")
    axes[2].set_title("Messages received per one-second window")
    axes[2].grid(True, axis="y", alpha=0.3)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=args.dpi)
    print_summary(timing, rate)
    print(f"plot={output_path}")
    if args.show:
        plt.show()
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
