#!/usr/bin/env python3
"""Summarize IMU excitation in a ROS 2 bag without using ground truth."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

import numpy as np
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bag", type=Path, required=True)
    parser.add_argument("--topic", default="/livox/mid360/imu")
    parser.add_argument("--block-samples", type=int, default=20)
    return parser.parse_args()


def percentiles(values: np.ndarray) -> dict[str, float]:
    return {
        "median": float(np.median(values)),
        "p90": float(np.percentile(values, 90.0)),
        "p95": float(np.percentile(values, 95.0)),
        "p99": float(np.percentile(values, 99.0)),
        "max": float(np.max(values)),
    }


def main() -> None:
    args = parse_args()
    timestamps: list[float] = []
    acceleration: list[tuple[float, float, float]] = []
    angular_velocity: list[tuple[float, float, float]] = []
    database_paths = sorted(args.bag.glob("*.db3"))
    if not database_paths:
        raise SystemExit(f"no ROS 2 SQLite database found in: {args.bag}")
    message_type = None
    for database_path in database_paths:
        connection = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True)
        try:
            topic_row = connection.execute(
                "SELECT id, type FROM topics WHERE name = ?", (args.topic,)
            ).fetchone()
            if topic_row is None:
                continue
            topic_id, topic_type = topic_row
            if message_type is None:
                message_type = get_message(topic_type)
            for timestamp_ns, serialized in connection.execute(
                "SELECT timestamp, data FROM messages WHERE topic_id = ? ORDER BY timestamp",
                (topic_id,),
            ):
                message = deserialize_message(serialized, message_type)
                timestamps.append(timestamp_ns * 1e-9)
                acceleration.append(
                    (message.linear_acceleration.x, message.linear_acceleration.y, message.linear_acceleration.z)
                )
                angular_velocity.append(
                    (message.angular_velocity.x, message.angular_velocity.y, message.angular_velocity.z)
                )
        finally:
            connection.close()
    if message_type is None:
        raise SystemExit(f"topic not found: {args.topic}")

    acc = np.asarray(acceleration, dtype=np.float64)
    gyro = np.asarray(angular_velocity, dtype=np.float64)
    time = np.asarray(timestamps, dtype=np.float64)
    if len(time) < max(3, args.block_samples):
        raise SystemExit("not enough IMU samples")

    acc_norm = np.linalg.norm(acc, axis=1)
    gyro_norm = np.linalg.norm(gyro, axis=1)
    acc_delta = np.linalg.norm(np.diff(acc, axis=0), axis=1)
    gyro_delta = np.linalg.norm(np.diff(gyro, axis=0), axis=1)
    block_count = len(acc) // args.block_samples
    blocks = acc[: block_count * args.block_samples].reshape(block_count, args.block_samples, 3)
    block_rms = np.sqrt(np.mean(np.sum((blocks - blocks.mean(axis=1, keepdims=True)) ** 2, axis=2), axis=1))

    result = {
        "bag": str(args.bag.resolve()),
        "topic": args.topic,
        "samples": int(len(time)),
        "duration_s": float(time[-1] - time[0]),
        "sample_dt_s": percentiles(np.diff(time)),
        "acceleration": {
            "mean_xyz": [float(value) for value in acc.mean(axis=0)],
            "std_xyz": [float(value) for value in acc.std(axis=0)],
            "norm": percentiles(acc_norm),
            "successive_delta_norm": percentiles(acc_delta),
            f"block_{args.block_samples}_sample_rms": percentiles(block_rms),
        },
        "angular_velocity": {
            "mean_xyz": [float(value) for value in gyro.mean(axis=0)],
            "std_xyz": [float(value) for value in gyro.std(axis=0)],
            "norm": percentiles(gyro_norm),
            "successive_delta_norm": percentiles(gyro_delta),
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
