#!/usr/bin/env python3
"""Summarize raw IMU and lidar changes over time in a ROS2 SQLite bag."""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
from pathlib import Path

import numpy as np
from rclpy.serialization import deserialize_message
from sensor_msgs.msg import Imu, PointCloud2


def topic_id(connection: sqlite3.Connection, name: str) -> int:
    row = connection.execute("SELECT id FROM topics WHERE name = ?", (name,)).fetchone()
    if row is None:
        raise ValueError(f"topic not found: {name}")
    return int(row[0])


def xyz_array(message: PointCloud2) -> np.ndarray:
    offsets = {field.name: field.offset for field in message.fields}
    if not all(name in offsets for name in ("x", "y", "z")):
        raise ValueError("PointCloud2 lacks x/y/z fields")
    count = int(message.width) * int(message.height)
    dtype = np.dtype(
        {
            "names": ["x", "y", "z"],
            "formats": ["<f4", "<f4", "<f4"],
            "offsets": [offsets["x"], offsets["y"], offsets["z"]],
            "itemsize": int(message.point_step),
        }
    )
    points = np.frombuffer(message.data, dtype=dtype, count=count)
    xyz = np.column_stack((points["x"], points["y"], points["z"]))
    return xyz[np.all(np.isfinite(xyz), axis=1)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bag", type=Path, required=True)
    parser.add_argument("--lidar-topic", required=True)
    parser.add_argument("--imu-topic", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cloud-sample-period", type=float, default=1.0)
    args = parser.parse_args()
    if not math.isfinite(args.cloud_sample_period) or args.cloud_sample_period <= 0.0:
        raise ValueError("cloud-sample-period must be positive")

    databases = sorted(args.bag.glob("*.db3"))
    if len(databases) != 1:
        raise ValueError(f"expected one db3 payload in {args.bag}, found {len(databases)}")
    connection = sqlite3.connect(f"file:{databases[0]}?mode=ro", uri=True)
    try:
        imu_id = topic_id(connection, args.imu_topic)
        lidar_id = topic_id(connection, args.lidar_topic)
        imu_samples: list[list[float]] = []
        for timestamp, payload in connection.execute(
            "SELECT timestamp, data FROM messages WHERE topic_id = ? ORDER BY timestamp", (imu_id,)
        ):
            message = deserialize_message(payload, Imu)
            imu_samples.append(
                [
                    float(timestamp) * 1e-9,
                    message.angular_velocity.x,
                    message.angular_velocity.y,
                    message.angular_velocity.z,
                    message.linear_acceleration.x,
                    message.linear_acceleration.y,
                    message.linear_acceleration.z,
                ]
            )

        cloud_samples: list[dict[str, object]] = []
        last_sample_time = -math.inf
        for timestamp, payload in connection.execute(
            "SELECT timestamp, data FROM messages WHERE topic_id = ? ORDER BY timestamp", (lidar_id,)
        ):
            time_s = float(timestamp) * 1e-9
            if time_s - last_sample_time < args.cloud_sample_period:
                continue
            message = deserialize_message(payload, PointCloud2)
            points = xyz_array(message)
            if not len(points):
                continue
            ranges = np.linalg.norm(points, axis=1)
            cloud_samples.append(
                {
                    "time_s": time_s,
                    "points": int(len(points)),
                    "centroid_m": np.mean(points, axis=0).tolist(),
                    "range_quantiles_m": np.percentile(ranges, [10.0, 50.0, 90.0]).tolist(),
                }
            )
            last_sample_time = time_s
    finally:
        connection.close()

    imu = np.asarray(imu_samples, dtype=np.float64)
    if len(imu) == 0 or len(cloud_samples) == 0:
        raise ValueError("IMU or lidar topic has no usable samples")
    origin = min(float(imu[0, 0]), float(cloud_samples[0]["time_s"]))
    gyro_norm = np.linalg.norm(imu[:, 1:4], axis=1)
    acceleration_norm = np.linalg.norm(imu[:, 4:7], axis=1)
    bins: list[dict[str, float]] = []
    second = math.floor(float(imu[0, 0]))
    last_second = math.floor(float(imu[-1, 0]))
    while second <= last_second:
        mask = (imu[:, 0] >= second) & (imu[:, 0] < second + 1.0)
        if np.any(mask):
            bins.append(
                {
                    "time_from_start_s": second + 0.5 - origin,
                    "gyro_mean_rad_s": float(np.mean(gyro_norm[mask])),
                    "gyro_p95_rad_s": float(np.percentile(gyro_norm[mask], 95.0)),
                    "acceleration_norm_mean_m_s2": float(np.mean(acceleration_norm[mask])),
                    "acceleration_norm_std_m_s2": float(np.std(acceleration_norm[mask])),
                }
            )
        second += 1

    for sample in cloud_samples:
        sample["time_from_start_s"] = float(sample["time_s"]) - origin
    centroids = np.asarray([sample["centroid_m"] for sample in cloud_samples], dtype=np.float64)
    descriptors = np.asarray(
        [sample["centroid_m"] + sample["range_quantiles_m"] for sample in cloud_samples], dtype=np.float64
    )
    descriptor_step = np.linalg.norm(np.diff(descriptors, axis=0), axis=1)
    report = {
        "bag": str(args.bag.resolve()),
        "imu_topic": args.imu_topic,
        "lidar_topic": args.lidar_topic,
        "imu_samples": int(len(imu)),
        "cloud_samples": len(cloud_samples),
        "duration_s": float(max(imu[-1, 0], cloud_samples[-1]["time_s"]) - origin),
        "gyro_norm_rad_s": {
            "median": float(np.median(gyro_norm)),
            "p95": float(np.percentile(gyro_norm, 95.0)),
            "max": float(np.max(gyro_norm)),
        },
        "cloud_centroid_axis_range_m": np.ptp(centroids, axis=0).tolist(),
        "cloud_descriptor_step": {
            "median": float(np.median(descriptor_step)),
            "p95": float(np.percentile(descriptor_step, 95.0)),
            "max": float(np.max(descriptor_step)),
        },
        "imu_one_second_bins": bins,
        "cloud_time_series": cloud_samples,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
