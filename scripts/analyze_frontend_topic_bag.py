#!/usr/bin/env python3
"""Audit the five online frontend topics from a ROS 2 bag recording."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import rosbag2_py
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message


POSE = "/slamPoseRaw_topic"
CLOUD = "/final_points_topic"
SAFETY = "/slamSafety_topic"
STATE = "/slamState_topic"
SYSTEM = "/SystemState"
CLOCK = "/clock"
EXPECTED = (POSE, CLOUD, SAFETY, STATE, SYSTEM)
EXPECTED_TYPES = {
    POSE: "geometry_msgs/msg/PoseStamped",
    CLOUD: "sensor_msgs/msg/PointCloud2",
    SAFETY: "std_msgs/msg/Float64",
    STATE: "std_msgs/msg/Float64",
    SYSTEM: "std_msgs/msg/Int32",
    CLOCK: "rosgraph_msgs/msg/Clock",
}
EXPECTED_CLOUD_SCHEMA = {
    "x": (0, 7, 1),
    "y": (4, 7, 1),
    "z": (8, 7, 1),
    "intensity": (12, 7, 1),
    "time": (16, 8, 1),
    "lidar_id": (24, 2, 1),
}


def stamp_ns(stamp: Any) -> int:
    return int(stamp.sec) * 1_000_000_000 + int(stamp.nanosec)


def distribution(values: np.ndarray) -> dict[str, float | int | None]:
    if not len(values):
        return {
            "count": 0, "mean": None, "min": None, "p01": None, "p05": None,
            "median": None, "p95": None, "p99": None, "max": None,
        }
    return {
        "count": int(len(values)),
        "mean": float(np.mean(values)),
        "min": float(np.min(values)),
        "p01": float(np.percentile(values, 1)),
        "p05": float(np.percentile(values, 5)),
        "median": float(np.median(values)),
        "p95": float(np.percentile(values, 95)),
        "p99": float(np.percentile(values, 99)),
        "max": float(np.max(values)),
    }


def timing_summary(record_ns: list[int]) -> dict[str, Any]:
    times = np.asarray(record_ns, dtype=np.int64)
    intervals = np.diff(times).astype(np.float64) * 1e-9
    duration = float((times[-1] - times[0]) * 1e-9) if len(times) >= 2 else 0.0
    interval_summary = distribution(intervals)
    median_interval = interval_summary["median"]
    return {
        "message_count": int(len(times)),
        "duration_s": duration,
        "average_rate_hz": float((len(times) - 1) / duration) if duration > 0.0 else None,
        "nominal_receive_rate_hz": (
            float(1.0 / median_interval) if median_interval and median_interval > 0.0 else None
        ),
        "record_interval_s": interval_summary,
    }


def header_summary(headers_ns: list[int], causal_clock_ages_s: list[float]) -> dict[str, Any]:
    headers = np.asarray(headers_ns, dtype=np.int64)
    differences = np.diff(headers)
    unique_headers = np.unique(headers)
    header_duration = float((unique_headers[-1] - unique_headers[0]) * 1e-9) if len(unique_headers) >= 2 else 0.0
    summary = {
        "unique_header_count": int(len(unique_headers)),
        "duplicate_header_count": int(len(headers) - len(unique_headers)),
        "duplicate_header_ratio": float((len(headers) - len(unique_headers)) / len(headers)) if len(headers) else None,
        "nonmonotonic_or_duplicate_count": int(np.count_nonzero(differences <= 0)),
        "unique_header_rate_hz": float((len(unique_headers) - 1) / header_duration) if header_duration > 0.0 else None,
        "header_interval_s": distribution(np.diff(unique_headers).astype(np.float64) * 1e-9),
    }
    ages = np.asarray(causal_clock_ages_s, dtype=np.float64)
    summary["ros_clock_match_ratio"] = float(len(ages) / len(headers)) if len(headers) else 0.0
    summary["ros_clock_data_age_s"] = distribution(ages)
    summary["ros_clock_absolute_data_age_s"] = distribution(np.abs(ages))
    return summary


def cloud_schema_errors(message: Any) -> list[str]:
    errors: list[str] = []
    names = [field.name for field in message.fields]
    if len(names) != len(set(names)):
        errors.append("duplicate field names")
    fields = {field.name: field for field in message.fields}
    if set(fields) != set(EXPECTED_CLOUD_SCHEMA):
        errors.append(f"field names={sorted(fields)}")
    for name, (offset, datatype, count) in EXPECTED_CLOUD_SCHEMA.items():
        field = fields.get(name)
        if field is None:
            continue
        if (field.offset, field.datatype, field.count) != (offset, datatype, count):
            errors.append(
                f"{name} schema={(field.offset, field.datatype, field.count)} expected={(offset, datatype, count)}"
            )
    if message.height != 1:
        errors.append(f"height={message.height}")
    if message.point_step != 25:
        errors.append(f"point_step={message.point_step}")
    if message.row_step != message.point_step * message.width:
        errors.append(f"row_step={message.row_step}")
    if len(message.data) != message.row_step * message.height:
        errors.append(f"data_size={len(message.data)} expected={message.row_step * message.height}")
    return errors


def point_array(message: Any) -> np.ndarray:
    fields = {field.name: field for field in message.fields}
    endian = ">" if message.is_bigendian else "<"
    dtype = np.dtype(
        {
            "names": ["x", "y", "z", "intensity", "time", "lidar_id"],
            "formats": [endian + "f4", endian + "f4", endian + "f4", endian + "f4", endian + "f8", "u1"],
            "offsets": [fields[name].offset for name in ("x", "y", "z", "intensity", "time", "lidar_id")],
            "itemsize": message.point_step,
        }
    )
    return np.frombuffer(message.data, dtype=dtype, count=message.width)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bag", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--enforce", dest="enforce", action="store_true", default=True)
    parser.add_argument("--diagnostic-only", dest="enforce", action="store_false")
    parser.add_argument("--min-duration-s", type=float, default=60.0)
    parser.add_argument("--min-message-count", type=int, default=500)
    parser.add_argument("--expected-lidar-ids", default="0,1,2,3")
    parser.add_argument("--min-nonempty-cloud-ratio", type=float, default=0.99)
    parser.add_argument("--min-exact-lidar-id-frame-ratio", type=float, default=0.99)
    parser.add_argument("--min-points-per-cloud", type=int, default=1000)
    parser.add_argument("--min-points-per-lidar-per-cloud", type=int, default=100)
    parser.add_argument("--min-frame-time-span-s", type=float, default=0.01)
    parser.add_argument("--min-content-valid-frame-ratio", type=float, default=0.99)
    parser.add_argument("--min-clock-match-ratio", type=float, default=0.98)
    parser.add_argument("--max-clock-age-p95-s", type=float, default=0.20)
    parser.add_argument("--max-clock-age-s", type=float, default=0.30)
    parser.add_argument("--min-normal-state-ratio", type=float, default=0.80)
    parser.add_argument("--max-lifecycle-reorder-s", type=float, default=0.15)
    args = parser.parse_args()
    try:
        expected_lidar_ids = tuple(sorted({int(value) for value in args.expected_lidar_ids.split(",")}))
    except ValueError as error:
        raise SystemExit(f"invalid --expected-lidar-ids: {error}") from error
    if not expected_lidar_ids or any(value < 0 or value > 3 for value in expected_lidar_ids):
        raise SystemExit("--expected-lidar-ids must be a non-empty subset of 0,1,2,3")
    if args.min_duration_s <= 0.0 or args.min_message_count < 2:
        raise SystemExit("minimum duration/count must be positive and count must be at least two")
    ratios = (
        args.min_nonempty_cloud_ratio,
        args.min_exact_lidar_id_frame_ratio,
        args.min_content_valid_frame_ratio,
        args.min_clock_match_ratio,
        args.min_normal_state_ratio,
    )
    if any(value < 0.0 or value > 1.0 for value in ratios):
        raise SystemExit("all ratio thresholds must be within [0,1]")
    if args.min_points_per_cloud < 1 or args.min_points_per_lidar_per_cloud < 1 or args.min_frame_time_span_s <= 0.0:
        raise SystemExit("cloud content thresholds must be positive")
    if args.max_lifecycle_reorder_s < 0.0:
        raise SystemExit("--max-lifecycle-reorder-s must be non-negative")

    reader = rosbag2_py.SequentialReader()
    reader.open(
        rosbag2_py.StorageOptions(uri=str(args.bag), storage_id="sqlite3"),
        rosbag2_py.ConverterOptions("", ""),
    )
    topic_types = {item.name: item.type for item in reader.get_all_topics_and_types()}
    missing = [topic for topic in EXPECTED if topic not in topic_types]
    if missing:
        raise SystemExit(f"recording lacks required topics: {missing}")
    wrong_types = {
        topic: {"actual": topic_types.get(topic), "expected": expected_type}
        for topic, expected_type in EXPECTED_TYPES.items()
        if topic in topic_types and topic_types[topic] != expected_type
    }
    if wrong_types:
        raise SystemExit(f"recording has incorrect topic types: {wrong_types}")
    available_topics = (*EXPECTED, *((CLOCK,) if CLOCK in topic_types else ()))
    message_classes = {topic: get_message(topic_types[topic]) for topic in available_topics}
    record_times: dict[str, list[int]] = {topic: [] for topic in (*EXPECTED, CLOCK)}
    headers: dict[str, list[int]] = {POSE: [], CLOUD: []}
    state_values: dict[str, list[float | int]] = {SAFETY: [], STATE: [], SYSTEM: []}
    latest_clock_ns: int | None = None
    last_clock_ns: int | None = None
    clock_values_ns: list[int] = []
    clock_nonmonotonic_count = 0
    causal_clock_ages: dict[str, list[float]] = {POSE: [], CLOUD: []}
    header_record_times: dict[str, dict[int, list[int]]] = {POSE: {}, CLOUD: {}}
    pose_frames: set[str] = set()
    cloud_frames: set[str] = set()
    quaternion_max_norm_error = 0.0
    cloud_fields: list[dict[str, Any]] | None = None
    cloud_point_count = 0
    cloud_bad_xyz = 0
    cloud_bad_intensity = 0
    cloud_bad_time = 0
    cloud_time_min = math.inf
    cloud_time_max = -math.inf
    cloud_bad_lidar_id = 0
    cloud_frame_count = 0
    cloud_nonempty_frame_count = 0
    cloud_frames_exact_expected_ids = 0
    cloud_schema_error_count = 0
    cloud_schema_error_examples: list[str] = []
    lidar_id_counts = np.zeros(4, dtype=np.int64)
    cloud_frame_point_counts: list[int] = []
    cloud_frame_time_spans: list[float] = []
    cloud_frame_lidar_counts: list[list[int]] = [[], [], [], []]
    cloud_content_valid_frame_count = 0
    pose_nonfinite_count = 0
    first_pose_translation_norm: float | None = None
    first_pose_rotation_angle_deg: float | None = None

    while reader.has_next():
        topic, serialized, record_time = reader.read_next()
        if topic not in message_classes:
            continue
        message = deserialize_message(serialized, message_classes[topic])
        record_times[topic].append(int(record_time))
        if topic == CLOCK:
            clock_ns = stamp_ns(message.clock)
            if last_clock_ns is not None and clock_ns <= last_clock_ns:
                clock_nonmonotonic_count += 1
            last_clock_ns = clock_ns
            latest_clock_ns = clock_ns
            clock_values_ns.append(clock_ns)
        elif topic == POSE:
            header_ns = stamp_ns(message.header.stamp)
            headers[POSE].append(header_ns)
            header_record_times[POSE].setdefault(header_ns, []).append(int(record_time))
            if latest_clock_ns is not None:
                causal_clock_ages[POSE].append((latest_clock_ns - header_ns) * 1e-9)
            pose_frames.add(message.header.frame_id)
            quaternion = message.pose.orientation
            values = (
                message.pose.position.x, message.pose.position.y, message.pose.position.z,
                quaternion.x, quaternion.y, quaternion.z, quaternion.w,
            )
            if not all(math.isfinite(value) for value in values):
                pose_nonfinite_count += 1
            else:
                norm = math.sqrt(quaternion.x**2 + quaternion.y**2 + quaternion.z**2 + quaternion.w**2)
                quaternion_max_norm_error = max(quaternion_max_norm_error, abs(norm - 1.0))
                if first_pose_translation_norm is None:
                    first_pose_translation_norm = math.sqrt(sum(value * value for value in values[:3]))
                    first_pose_rotation_angle_deg = math.degrees(
                        2.0 * math.acos(min(1.0, abs(quaternion.w) / norm))
                    ) if norm > 0.0 else math.inf
        elif topic == CLOUD:
            header_ns = stamp_ns(message.header.stamp)
            headers[CLOUD].append(header_ns)
            header_record_times[CLOUD].setdefault(header_ns, []).append(int(record_time))
            if latest_clock_ns is not None:
                causal_clock_ages[CLOUD].append((latest_clock_ns - header_ns) * 1e-9)
            cloud_frames.add(message.header.frame_id)
            cloud_frame_count += 1
            if cloud_fields is None:
                cloud_fields = [
                    {"name": field.name, "offset": field.offset, "datatype": field.datatype, "count": field.count}
                    for field in message.fields
                ]
            schema_errors = cloud_schema_errors(message)
            if schema_errors:
                cloud_schema_error_count += 1
                if len(cloud_schema_error_examples) < 5:
                    cloud_schema_error_examples.append("; ".join(schema_errors))
                continue
            points = point_array(message)
            cloud_nonempty_frame_count += int(len(points) > 0)
            cloud_frame_point_counts.append(len(points))
            cloud_point_count += len(points)
            cloud_bad_xyz += int(np.count_nonzero(~np.isfinite(points["x"]) | ~np.isfinite(points["y"]) | ~np.isfinite(points["z"])))
            cloud_bad_intensity += int(np.count_nonzero(~np.isfinite(points["intensity"])))
            point_times = points["time"]
            finite_times = point_times[np.isfinite(point_times)]
            cloud_bad_time += int(len(point_times) - len(finite_times))
            if len(finite_times):
                cloud_time_min = min(cloud_time_min, float(np.min(finite_times)))
                cloud_time_max = max(cloud_time_max, float(np.max(finite_times)))
                frame_time_span = float(np.max(finite_times) - np.min(finite_times))
                cloud_frame_time_spans.append(frame_time_span)
                cloud_bad_time += int(np.count_nonzero((finite_times < -0.5) | (finite_times > 0.01)))
            else:
                frame_time_span = 0.0
                cloud_frame_time_spans.append(frame_time_span)
            ids = points["lidar_id"].astype(np.int64)
            cloud_bad_lidar_id += int(np.count_nonzero((ids < 0) | (ids > 3)))
            valid_ids = ids[(ids >= 0) & (ids <= 3)]
            counts = np.bincount(valid_ids, minlength=4)[:4]
            lidar_id_counts += counts
            for lidar_id in range(4):
                cloud_frame_lidar_counts[lidar_id].append(int(counts[lidar_id]))
            present_ids = tuple(np.flatnonzero(counts > 0).tolist())
            cloud_frames_exact_expected_ids += int(present_ids == expected_lidar_ids)
            content_valid = (
                len(points) >= args.min_points_per_cloud
                and all(counts[lidar_id] >= args.min_points_per_lidar_per_cloud for lidar_id in expected_lidar_ids)
                and frame_time_span >= args.min_frame_time_span_s
            )
            cloud_content_valid_frame_count += int(content_valid)
        else:
            state_values[topic].append(message.data)

    topic_summary = {topic: timing_summary(record_times[topic]) for topic in EXPECTED}
    topic_summary[POSE]["header"] = header_summary(headers[POSE], causal_clock_ages[POSE])
    topic_summary[CLOUD]["header"] = header_summary(headers[CLOUD], causal_clock_ages[CLOUD])
    pose_unique = np.unique(np.asarray(headers[POSE], dtype=np.int64))
    cloud_unique = np.unique(np.asarray(headers[CLOUD], dtype=np.int64))
    paired = np.intersect1d(pose_unique, cloud_unique)
    pair_denominator = max(len(pose_unique), len(cloud_unique), 1)
    pair_arrival_deltas = np.asarray(
        [abs(header_record_times[POSE][int(stamp)][0] - header_record_times[CLOUD][int(stamp)][0]) * 1e-9 for stamp in paired],
        dtype=np.float64,
    )
    topic_starts = [record_times[topic][0] for topic in EXPECTED if record_times[topic]]
    topic_ends = [record_times[topic][-1] for topic in EXPECTED if record_times[topic]]
    common_record_duration_s = (
        max(0.0, (min(topic_ends) - max(topic_starts)) * 1e-9)
        if len(topic_starts) == len(EXPECTED) else 0.0
    )

    safety_values = np.asarray(state_values[SAFETY], dtype=np.float64)
    state_array = np.asarray(state_values[STATE], dtype=np.float64)
    system_array = np.asarray(state_values[SYSTEM], dtype=np.int64)
    first_state_one_record_ns = next(
        (record_time for record_time, value in zip(record_times[STATE], state_array) if value == 1.0), None
    )
    first_system_one_record_ns = next(
        (record_time for record_time, value in zip(record_times[SYSTEM], system_array) if value == 1), None
    )
    lifecycle_first_one_delta_s = (
        (first_state_one_record_ns - first_system_one_record_ns) * 1e-9
        if first_state_one_record_ns is not None and first_system_one_record_ns is not None
        else None
    )
    if record_times[POSE] and record_times[CLOUD] and len(state_array):
        active_start_ns = max(record_times[POSE][0], record_times[CLOUD][0])
        active_end_ns = min(record_times[POSE][-1], record_times[CLOUD][-1])
        state_record_ns = np.asarray(record_times[STATE], dtype=np.int64)
        active_state = state_array[(state_record_ns >= active_start_ns) & (state_record_ns <= active_end_ns)]
    else:
        active_start_ns = 0
        active_end_ns = 0
        active_state = np.asarray([], dtype=np.float64)
    checks = {
        "pose_frame_ids": sorted(pose_frames),
        "cloud_frame_ids": sorted(cloud_frames),
        "quaternion_max_norm_error": quaternion_max_norm_error,
        "pose_nonfinite_count": pose_nonfinite_count,
        "first_rear_axle_pose_translation_norm_m": first_pose_translation_norm,
        "first_rear_axle_pose_rotation_angle_deg": first_pose_rotation_angle_deg,
        "pose_cloud_exact_unique_stamp_pair_count": int(len(paired)),
        "pose_cloud_exact_unique_stamp_pair_ratio": float(len(paired) / pair_denominator),
        "pose_cloud_record_arrival_delta_s": distribution(pair_arrival_deltas),
        "common_record_duration_s": common_record_duration_s,
        "clock_nonmonotonic_count": clock_nonmonotonic_count,
        "clock_value_duration_s": (
            (clock_values_ns[-1] - clock_values_ns[0]) * 1e-9 if len(clock_values_ns) >= 2 else 0.0
        ),
        "clock_first_value_ns": clock_values_ns[0] if clock_values_ns else None,
        "clock_last_value_ns": clock_values_ns[-1] if clock_values_ns else None,
        "safety_values_are_binary": bool(np.all(np.isin(safety_values, (0.0, 1.0)))),
        "safety_strict_alternation_ratio": (
            float(np.mean(np.diff(safety_values) != 0.0)) if len(safety_values) >= 2 else None
        ),
        "slam_state_values_are_binary": bool(np.all(np.isin(state_array, (0.0, 1.0)))),
        "slam_state_normal_ratio": float(np.mean(state_array == 1.0)) if len(state_array) else 0.0,
        "slam_state_final_value": float(state_array[-1]) if len(state_array) else None,
        "slam_state_active_window_start_record_ns": active_start_ns or None,
        "slam_state_active_window_end_record_ns": active_end_ns or None,
        "slam_state_active_sample_count": int(len(active_state)),
        "slam_state_active_normal_ratio": float(np.mean(active_state == 1.0)) if len(active_state) else 0.0,
        "slam_state_first_one_record_ns": first_state_one_record_ns,
        "slam_state_minus_system_state_first_one_receive_s": lifecycle_first_one_delta_s,
        "system_state_values_are_binary": bool(np.all(np.isin(system_array, (0, 1)))),
        "system_state_transition_count": int(np.count_nonzero(np.diff(system_array) != 0)) if len(system_array) >= 2 else 0,
        "system_state_zero_to_one_count": int(np.count_nonzero((system_array[:-1] == 0) & (system_array[1:] == 1))) if len(system_array) >= 2 else 0,
        "system_state_first_one_record_ns": first_system_one_record_ns,
        "cloud_fields": cloud_fields,
        "cloud_point_count": cloud_point_count,
        "cloud_bad_xyz_count": cloud_bad_xyz,
        "cloud_bad_intensity_count": cloud_bad_intensity,
        "cloud_bad_time_count": cloud_bad_time,
        "cloud_time_min_s": None if not math.isfinite(cloud_time_min) else cloud_time_min,
        "cloud_time_max_s": None if not math.isfinite(cloud_time_max) else cloud_time_max,
        "cloud_bad_lidar_id_count": cloud_bad_lidar_id,
        "cloud_frame_count": cloud_frame_count,
        "cloud_nonempty_frame_count": cloud_nonempty_frame_count,
        "cloud_nonempty_frame_ratio": cloud_nonempty_frame_count / cloud_frame_count if cloud_frame_count else 0.0,
        "expected_lidar_ids": list(expected_lidar_ids),
        "cloud_frames_exact_expected_ids": cloud_frames_exact_expected_ids,
        "cloud_frames_exact_expected_ids_ratio": cloud_frames_exact_expected_ids / cloud_frame_count if cloud_frame_count else 0.0,
        "cloud_frame_point_count": distribution(np.asarray(cloud_frame_point_counts, dtype=np.float64)),
        "cloud_frame_time_span_s": distribution(np.asarray(cloud_frame_time_spans, dtype=np.float64)),
        "cloud_frame_point_count_by_lidar": {
            str(lidar_id): distribution(np.asarray(cloud_frame_lidar_counts[lidar_id], dtype=np.float64))
            for lidar_id in range(4)
        },
        "cloud_content_valid_frame_count": cloud_content_valid_frame_count,
        "cloud_content_valid_frame_ratio": cloud_content_valid_frame_count / cloud_frame_count if cloud_frame_count else 0.0,
        "cloud_schema_error_count": cloud_schema_error_count,
        "cloud_schema_error_examples": cloud_schema_error_examples,
        "lidar_id_point_counts": lidar_id_counts.tolist(),
    }
    failures = []
    realtime_failures = []
    for topic in EXPECTED:
        timing = topic_summary[topic]
        interval = timing["record_interval_s"]
        if timing["message_count"] < args.min_message_count or timing["duration_s"] < args.min_duration_s:
            failures.append(f"{topic}: recording is shorter than the minimum sample contract")
        nominal_rate = timing["nominal_receive_rate_hz"]
        if nominal_rate is None or not 9.5 <= nominal_rate <= 10.5:
            failures.append(f"{topic}: nominal receive cadence outside [9.5,10.5] Hz")
        if interval["p99"] is None or interval["p99"] > 0.21:
            realtime_failures.append(f"{topic}: 99th-percentile receive interval exceeds 0.21 s")
        if timing["average_rate_hz"] is None or not 9.5 <= timing["average_rate_hz"] <= 10.5:
            realtime_failures.append(f"{topic}: end-to-end average receive rate outside [9.5,10.5] Hz")
        if interval["max"] is None or interval["max"] > 0.30:
            realtime_failures.append(f"{topic}: worst receive interval exceeds 0.30 s")
    for topic in (POSE, CLOUD):
        header = topic_summary[topic]["header"]
        if header["unique_header_rate_hz"] is None or not 9.5 <= header["unique_header_rate_hz"] <= 10.5:
            failures.append(f"{topic}: unique header rate outside [9.5,10.5] Hz")
        header_interval = header["header_interval_s"]
        if header_interval["p95"] is None or header_interval["p95"] > 0.13 or header_interval["p99"] > 0.20 or header_interval["max"] > 0.30:
            failures.append(f"{topic}: unique header interval contract failed")
        if header["duplicate_header_ratio"] is None or header["duplicate_header_ratio"] > 0.01:
            failures.append(f"{topic}: duplicate header ratio exceeds 1%")
        if header["nonmonotonic_or_duplicate_count"] != 0:
            failures.append(f"{topic}: header stamps are not strictly increasing")
        clock_age = header["ros_clock_absolute_data_age_s"]
        if header["ros_clock_match_ratio"] < args.min_clock_match_ratio:
            failures.append(f"{topic}: causal /clock match ratio is too low")
        if clock_age["p95"] is None or clock_age["p95"] > args.max_clock_age_p95_s:
            failures.append(f"{topic}: ROS-clock data age contract failed")
        if clock_age["max"] is None or clock_age["max"] > args.max_clock_age_s:
            realtime_failures.append(f"{topic}: worst ROS-clock data age exceeds {args.max_clock_age_s:.3f} s")
    if CLOCK not in topic_types:
        failures.append("recording lacks /clock")
    if checks["clock_nonmonotonic_count"]:
        failures.append("/clock is not strictly increasing")
    output_header_min = min((values[0] for values in headers.values() if values), default=None)
    output_header_max = max((values[-1] for values in headers.values() if values), default=None)
    if (
        checks["clock_value_duration_s"] < args.min_duration_s
        or checks["clock_first_value_ns"] is None
        or checks["clock_last_value_ns"] is None
        or output_header_min is None
        or output_header_max is None
        or checks["clock_first_value_ns"] > output_header_min + 200_000_000
        or checks["clock_last_value_ns"] < output_header_max - 200_000_000
    ):
        failures.append("/clock does not span the validated pose/cloud sensor-time window")
    if checks["common_record_duration_s"] < args.min_duration_s:
        failures.append("five topics do not share the minimum common recording window")
    if checks["pose_cloud_exact_unique_stamp_pair_ratio"] < 0.98:
        failures.append("pose/cloud exact unique stamp pair ratio below 98%")
    if checks["quaternion_max_norm_error"] > 1e-3:
        failures.append("pose quaternion norm error exceeds 1e-3")
    if checks["pose_nonfinite_count"]:
        failures.append("pose contains non-finite position or quaternion values")
    if (
        checks["first_rear_axle_pose_translation_norm_m"] is None
        or checks["first_rear_axle_pose_translation_norm_m"] > 0.001
        or checks["first_rear_axle_pose_rotation_angle_deg"] is None
        or checks["first_rear_axle_pose_rotation_angle_deg"] > 0.05
    ):
        failures.append("first valid rear-axle pose is not the required identity origin")
    pair_arrival = checks["pose_cloud_record_arrival_delta_s"]
    if pair_arrival["p95"] is None or pair_arrival["p95"] > 0.13:
        failures.append("pose/cloud paired-message arrival delay contract failed")
    if pair_arrival["p99"] is None or pair_arrival["p99"] > 0.20 or pair_arrival["max"] is None or pair_arrival["max"] > 0.30:
        realtime_failures.append("pose/cloud paired-message tail arrival delay exceeds the realtime diagnostic limit")
    if checks["cloud_schema_error_count"]:
        failures.append("cloud schema is inconsistent with the six-field public contract")
    if checks["cloud_bad_xyz_count"] or checks["cloud_bad_intensity_count"] or checks["cloud_bad_time_count"] or checks["cloud_bad_lidar_id_count"]:
        failures.append("cloud contains invalid XYZ, intensity, relative time, or lidar_id")
    if checks["cloud_nonempty_frame_ratio"] < args.min_nonempty_cloud_ratio:
        failures.append("non-empty cloud frame ratio is too low")
    if checks["cloud_frames_exact_expected_ids_ratio"] < args.min_exact_lidar_id_frame_ratio:
        failures.append("cloud frames do not contain exactly the configured lidar-id set")
    if checks["cloud_content_valid_frame_ratio"] < args.min_content_valid_frame_ratio:
        failures.append("cloud frames do not satisfy per-frame point-count or temporal-span content thresholds")
    if any(checks["lidar_id_point_counts"][lidar_id] <= 0 for lidar_id in expected_lidar_ids):
        failures.append("one or more expected lidars contributed no points")
    if checks["cloud_time_min_s"] is None or checks["cloud_time_max_s"] is None or checks["cloud_time_max_s"] - checks["cloud_time_min_s"] < 0.01:
        failures.append("cloud relative-time field has no meaningful temporal span")
    if checks["pose_frame_ids"] != ["map"] or checks["cloud_frame_ids"] != ["lidar_114"]:
        failures.append("frame_id contract failed")
    if not checks["safety_values_are_binary"] or checks["safety_strict_alternation_ratio"] != 1.0:
        failures.append("safety heartbeat is not a strict binary alternation")
    if not checks["slam_state_values_are_binary"] or not checks["system_state_values_are_binary"]:
        failures.append("state topics contain non-binary values")
    if checks["slam_state_active_normal_ratio"] < args.min_normal_state_ratio:
        failures.append("slamState does not sustain normal operation during the pose/cloud activity window")
    if checks["system_state_zero_to_one_count"] != 1 or checks["system_state_transition_count"] != 1:
        failures.append("SystemState does not contain exactly one 0-to-1 transition")
    if (
        checks["slam_state_first_one_record_ns"] is None
        or checks["system_state_first_one_record_ns"] is None
        or checks["slam_state_minus_system_state_first_one_receive_s"] < -args.max_lifecycle_reorder_s
    ):
        failures.append("slamState/SystemState initialization order exceeds the allowed cross-topic receive reordering")

    result = {
        "bag": str(args.bag),
        "topic_types": {topic: topic_types[topic] for topic in available_topics},
        "clock_topic_available": CLOCK in topic_types,
        "contract_parameters": {
            "time_basis": "sensor-time cadence is the hard 10 Hz contract; 1x wall receive throughput is reported separately because the approved first version may run below real time",
            "min_duration_s": args.min_duration_s,
            "min_message_count": args.min_message_count,
            "expected_lidar_ids": list(expected_lidar_ids),
            "min_nonempty_cloud_ratio": args.min_nonempty_cloud_ratio,
            "min_exact_lidar_id_frame_ratio": args.min_exact_lidar_id_frame_ratio,
            "min_points_per_cloud": args.min_points_per_cloud,
            "min_points_per_lidar_per_cloud": args.min_points_per_lidar_per_cloud,
            "min_frame_time_span_s": args.min_frame_time_span_s,
            "min_content_valid_frame_ratio": args.min_content_valid_frame_ratio,
            "min_clock_match_ratio": args.min_clock_match_ratio,
            "max_clock_age_p95_s": args.max_clock_age_p95_s,
            "max_clock_age_s": args.max_clock_age_s,
            "min_normal_state_ratio": args.min_normal_state_ratio,
            "max_lifecycle_reorder_s": args.max_lifecycle_reorder_s,
            "enforced": args.enforce,
        },
        "topic_timing": topic_summary,
        "message_contract": checks,
        "contract_failures": failures,
        "contract_passed": not failures,
        "realtime_contract_failures": realtime_failures,
        "realtime_contract_passed": not realtime_failures,
        "note": "pose/cloud are emitted once per unique LIO output; safety/state/system use an independent 100 ms publisher",
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("topic", "message_count", "duration_s", "average_rate_hz", "nominal_receive_rate_hz", "interval_median_s", "interval_p95_s", "interval_p99_s", "interval_max_s", "unique_header_rate_hz", "duplicate_header_ratio"),
        )
        writer.writeheader()
        for topic in EXPECTED:
            timing = topic_summary[topic]
            header = timing.get("header", {})
            writer.writerow(
                {
                    "topic": topic,
                    "message_count": timing["message_count"],
                    "duration_s": timing["duration_s"],
                    "average_rate_hz": timing["average_rate_hz"],
                    "nominal_receive_rate_hz": timing["nominal_receive_rate_hz"],
                    "interval_median_s": timing["record_interval_s"]["median"],
                    "interval_p95_s": timing["record_interval_s"]["p95"],
                    "interval_p99_s": timing["record_interval_s"]["p99"],
                    "interval_max_s": timing["record_interval_s"]["max"],
                    "unique_header_rate_hz": header.get("unique_header_rate_hz"),
                    "duplicate_header_ratio": header.get("duplicate_header_ratio"),
                }
            )
    print(json.dumps({"contract_passed": not failures, "failures": failures}, ensure_ascii=False))
    return 2 if args.enforce and failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
