#!/usr/bin/env python3
"""Plot /PosRes X/Y/speed time series or a 2D XY trajectory from MCAP bags."""

from __future__ import annotations

import argparse
import math
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import matplotlib.pyplot as plt
import numpy as np

try:
    from geosun_msgs.msg import PosRes
    from rclpy.serialization import deserialize_message
except ImportError as error:
    raise SystemExit(
        "ROS 2 Python messages are unavailable. Source /opt/ros/humble/setup.bash and "
        "the lightning-lm install/setup.bash before running this script."
    ) from error


MCAP_MAGIC = b"\x89MCAP0\r\n"
MCAP_CHANNEL = 0x04
MCAP_MESSAGE = 0x05
MCAP_CHUNK = 0x06
MESSAGE_HEADER_SIZE = struct.calcsize("<HIQQ")


@dataclass(frozen=True)
class Sample:
    record_time_ns: int
    message_time_ns: int
    x_m: float
    y_m: float
    speed_mps: float


@dataclass
class Series:
    path: Path
    label: str
    samples: list[Sample]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="MCAP file, rosbag directory, or directory recursively containing rosbag MCAP files",
    )
    parser.add_argument("--topic", default="/PosRes", help="topic to read (default: %(default)s)")
    parser.add_argument(
        "--time-source",
        choices=("message", "record"),
        default="message",
        help="time axis source (default: message header timestamp)",
    )
    parser.add_argument(
        "--plot",
        choices=("time", "xy"),
        default="time",
        help="plot X/Y/speed over time or the 2D XY trajectory (default: %(default)s)",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("posres_xy_speed.png"), help="output image"
    )
    parser.add_argument("--dpi", type=int, default=180, help="output image DPI")
    parser.add_argument("--show", action="store_true", help="also open an interactive window")
    return parser.parse_args()


def discover_mcap_files(inputs: Sequence[Path]) -> list[Path]:
    discovered: list[Path] = []
    for input_path in inputs:
        if input_path.is_file():
            if input_path.suffix.lower() != ".mcap":
                raise ValueError(f"input file is not MCAP: {input_path}")
            candidates = [input_path]
        elif input_path.is_dir():
            direct = sorted(input_path.glob("*.mcap"))
            candidates = direct if direct else sorted(input_path.rglob("*.mcap"))
            if not candidates:
                raise ValueError(f"no MCAP files found under: {input_path}")
        else:
            raise ValueError(f"input does not exist: {input_path}")
        discovered.extend(candidates)

    unique: list[Path] = []
    seen: set[Path] = set()
    for path in discovered:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(resolved)
    return unique


def read_u32_prefixed_string(data: bytes, offset: int) -> tuple[str, int]:
    if offset + 4 > len(data):
        raise ValueError("truncated MCAP string length")
    length = struct.unpack_from("<I", data, offset)[0]
    begin = offset + 4
    end = begin + length
    if end > len(data):
        raise ValueError("truncated MCAP string")
    return data[begin:end].decode("utf-8"), end


def iter_mcap_records(path: Path) -> Iterator[tuple[int, bytes]]:
    payload = path.read_bytes()
    if len(payload) < 2 * len(MCAP_MAGIC) or payload[:8] != MCAP_MAGIC or payload[-8:] != MCAP_MAGIC:
        raise ValueError(f"invalid or incomplete MCAP file: {path}")

    offset = len(MCAP_MAGIC)
    limit = len(payload) - len(MCAP_MAGIC)
    while offset < limit:
        if offset + 9 > limit:
            raise ValueError(f"truncated MCAP record header at byte {offset}: {path}")
        opcode = payload[offset]
        length = struct.unpack_from("<Q", payload, offset + 1)[0]
        body_begin = offset + 9
        body_end = body_begin + length
        if body_end > limit:
            raise ValueError(f"truncated MCAP record body at byte {offset}: {path}")
        yield opcode, payload[body_begin:body_end]
        offset = body_end


def read_posres_mcap(path: Path, topic: str) -> list[Sample]:
    """Read an unchunked MCAP directly, avoiding a rosbag2 MCAP plugin dependency."""
    channel_topics: dict[int, str] = {}
    samples: list[Sample] = []
    for opcode, body in iter_mcap_records(path):
        if opcode == MCAP_CHANNEL:
            if len(body) < 4:
                raise ValueError(f"invalid MCAP channel record: {path}")
            channel_id, _schema_id = struct.unpack_from("<HH", body, 0)
            channel_topic, _ = read_u32_prefixed_string(body, 4)
            channel_topics[channel_id] = channel_topic
        elif opcode == MCAP_CHUNK:
            raise ValueError(
                f"chunked MCAP is not supported by the direct reader: {path}. "
                "Install ros-humble-rosbag2-storage-mcap or record with MCAP chunking disabled."
            )
        elif opcode == MCAP_MESSAGE:
            if len(body) < MESSAGE_HEADER_SIZE:
                raise ValueError(f"invalid MCAP message record: {path}")
            channel_id, _sequence, log_time_ns, _publish_time_ns = struct.unpack_from(
                "<HIQQ", body, 0
            )
            if channel_topics.get(channel_id) != topic:
                continue
            message = deserialize_message(body[MESSAGE_HEADER_SIZE:], PosRes)
            message_time_ns = (
                int(message.header.stamp.sec) * 1_000_000_000 + int(message.header.stamp.nanosec)
            )
            samples.append(
                Sample(
                    record_time_ns=int(log_time_ns),
                    message_time_ns=message_time_ns,
                    x_m=float(message.f8enh[0]),
                    y_m=float(message.f8enh[1]),
                    speed_mps=float(message.f8vehiclespeed),
                )
            )
    if not samples:
        raise ValueError(f"topic {topic!r} has no messages in {path}")
    return samples


def make_unique_labels(paths: Sequence[Path]) -> list[str]:
    base_labels = [path.parent.name for path in paths]
    labels: list[str] = []
    for index, base in enumerate(base_labels):
        if base_labels.count(base) == 1:
            labels.append(base)
        else:
            labels.append(f"{base}/{paths[index].name}")
    return labels


def selected_time_ns(sample: Sample, source: str) -> int:
    if source == "record" or sample.message_time_ns <= 0:
        return sample.record_time_ns
    return sample.message_time_ns


def finite_samples(samples: Iterable[Sample], source: str) -> list[Sample]:
    result = [
        sample
        for sample in samples
        if selected_time_ns(sample, source) > 0
        and all(math.isfinite(value) for value in (sample.x_m, sample.y_m, sample.speed_mps))
    ]
    return sorted(result, key=lambda sample: selected_time_ns(sample, source))


def sample_key(sample: Sample) -> tuple[int, float, float, float]:
    return (sample.message_time_ns, sample.x_m, sample.y_m, sample.speed_mps)


def print_summary(series: Sequence[Series], time_source: str) -> None:
    previous_keys: set[tuple[int, float, float, float]] = set()
    for item in series:
        times = np.asarray([selected_time_ns(sample, time_source) for sample in item.samples])
        x = np.asarray([sample.x_m for sample in item.samples])
        y = np.asarray([sample.y_m for sample in item.samples])
        speed = np.asarray([sample.speed_mps for sample in item.samples])
        duration_s = (times[-1] - times[0]) * 1e-9 if len(times) > 1 else 0.0
        path_length_m = float(np.sum(np.hypot(np.diff(x), np.diff(y))))
        keys = {sample_key(sample) for sample in item.samples}
        duplicate_count = len(keys & previous_keys)
        duplicate_text = f" duplicate_with_earlier={duplicate_count}" if previous_keys else ""
        print(
            f"{item.label}: samples={len(item.samples)} duration_s={duration_s:.3f} "
            f"xy_path_m={path_length_m:.3f} speed_mps=[{speed.min():.3f}, {speed.max():.3f}]"
            f"{duplicate_text}"
        )
        previous_keys.update(keys)


def plot_series(
    series: Sequence[Series], time_source: str, output: Path, dpi: int, show: bool
) -> None:
    palette = ["#2563eb", "#ea580c", "#64748b", "#ca8a04", "#9333ea"]
    line_styles = ["-", "--", "-.", ":"]
    origin_ns = min(selected_time_ns(item.samples[0], time_source) for item in series)

    figure, axes = plt.subplots(
        3, 1, figsize=(13, 9), sharex=True, constrained_layout=True
    )
    for index, item in enumerate(series):
        color = palette[index % len(palette)]
        line_style = line_styles[index % len(line_styles)]
        x = np.asarray([sample.x_m for sample in item.samples])
        y = np.asarray([sample.y_m for sample in item.samples])
        speed = np.asarray([sample.speed_mps for sample in item.samples])
        elapsed_s = np.asarray(
            [(selected_time_ns(sample, time_source) - origin_ns) * 1e-9 for sample in item.samples]
        )
        duration_s = elapsed_s[-1] - elapsed_s[0] if len(elapsed_s) > 1 else 0.0
        legend_label = f"{item.label} (n={len(item.samples)}, {duration_s:.1f} s)"

        axes[0].plot(
            elapsed_s,
            x,
            color=color,
            linestyle=line_style,
            linewidth=1.35,
            alpha=0.9,
            label=legend_label,
        )
        axes[1].plot(
            elapsed_s,
            y,
            color=color,
            linestyle=line_style,
            linewidth=1.35,
            alpha=0.9,
            label=item.label,
        )
        axes[2].plot(
            elapsed_s,
            speed,
            color=color,
            linestyle=line_style,
            linewidth=1.35,
            alpha=0.9,
            label=item.label,
        )

    source_label = "message header" if time_source == "message" else "MCAP record"
    axes[0].set_title("X over time")
    axes[0].set_ylabel("X / East (m)")
    axes[0].legend(loc="best", fontsize=8)
    axes[1].set_title("Y over time")
    axes[1].set_ylabel("Y / North (m)")
    axes[2].set_title("Published signed vehicle speed over time")
    axes[2].set_ylabel("Speed (m/s)")
    axes[2].set_xlabel(f"Elapsed {source_label} time from earliest sample (s)")
    axes[2].axhline(0.0, color="#334155", linewidth=0.8, alpha=0.8)
    for axis in axes:
        axis.grid(True, color="#d1d5db", linewidth=0.7, alpha=0.7)

    figure.suptitle("/PosRes X, Y, and velocity over time", fontsize=15, color="#111827")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi, facecolor="white")
    if show:
        plt.show()
    plt.close(figure)


def plot_xy(series: Sequence[Series], time_source: str, output: Path, dpi: int, show: bool) -> None:
    palette = ["#2563eb", "#ea580c", "#64748b", "#ca8a04", "#9333ea"]
    line_styles = ["-", "--", "-.", ":"]
    figure, axis = plt.subplots(figsize=(10, 8), constrained_layout=True)

    for index, item in enumerate(series):
        color = palette[index % len(palette)]
        line_style = line_styles[index % len(line_styles)]
        x = np.asarray([sample.x_m for sample in item.samples])
        y = np.asarray([sample.y_m for sample in item.samples])
        times = np.asarray([selected_time_ns(sample, time_source) for sample in item.samples])
        duration_s = (times[-1] - times[0]) * 1e-9 if len(times) > 1 else 0.0
        label = f"{item.label} (n={len(item.samples)}, {duration_s:.1f} s)"

        axis.plot(
            x,
            y,
            color=color,
            linestyle=line_style,
            linewidth=1.5,
            alpha=0.9,
            label=label,
        )
        axis.scatter(
            x[0],
            y[0],
            s=58,
            marker="o",
            facecolor="white",
            edgecolor=color,
            linewidth=1.6,
            zorder=3,
        )
        axis.scatter(
            x[-1],
            y[-1],
            s=52,
            marker="s",
            facecolor=color,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )

    axis.set_title("/PosRes XY trajectory")
    axis.set_xlabel("X / East (m)")
    axis.set_ylabel("Y / North (m)")
    axis.set_aspect("equal", adjustable="datalim")
    axis.grid(True, color="#d1d5db", linewidth=0.7, alpha=0.7)
    axis.legend(loc="best", fontsize=8)
    axis.text(
        0.01,
        0.01,
        "open circle: start   filled square: end",
        transform=axis.transAxes,
        fontsize=8,
        color="#475569",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi, facecolor="white")
    if show:
        plt.show()
    plt.close(figure)


def main() -> int:
    args = parse_args()
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive")
    paths = discover_mcap_files(args.inputs)
    labels = make_unique_labels(paths)
    series: list[Series] = []
    for path, label in zip(paths, labels):
        samples = finite_samples(read_posres_mcap(path, args.topic), args.time_source)
        if not samples:
            raise ValueError(f"topic {args.topic!r} has no finite XY/speed samples in {path}")
        series.append(Series(path=path, label=label, samples=samples))

    if args.plot == "xy":
        plot_xy(series, args.time_source, args.output, args.dpi, args.show)
    else:
        plot_series(series, args.time_source, args.output, args.dpi, args.show)
    print_summary(series, args.time_source)
    print(f"plot={args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
