#!/usr/bin/env python3
"""Audit fail-safe behavior and recovery from a deterministic SANY sensor outage."""

from __future__ import annotations

import argparse
import hashlib
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
EXPECTED = (POSE, CLOUD, SAFETY, STATE, SYSTEM, CLOCK)
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
MIN_WINDOW_HEADER_RATE_HZ = 8.0
NOMINAL_HEADER_RATE_HZ = 9.5
MIN_WINDOW_SPAN_RATIO = 0.8
MAX_CONTINUOUS_HEADER_GAP_S = 0.30
MAX_LIFECYCLE_REORDER_S = 0.15


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bag_file_identities(bag: Path) -> list[dict[str, Any]]:
    files = sorted(
        path for path in bag.iterdir()
        if path.is_file() and (path.name == "metadata.yaml" or path.suffix in {".db3", ".mcap"})
    )
    return [
        {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)}
        for path in files
    ]


def stamp_ns(stamp: Any) -> int:
    return int(stamp.sec) * 1_000_000_000 + int(stamp.nanosec)


def rate_hz(times: list[int]) -> float | None:
    if len(times) < 2 or times[-1] <= times[0]:
        return None
    return (len(times) - 1) / ((times[-1] - times[0]) * 1e-9)


def cloud_schema_errors(message: Any) -> list[str]:
    errors: list[str] = []
    names = [field.name for field in message.fields]
    fields = {field.name: field for field in message.fields}
    if len(names) != len(set(names)):
        errors.append("duplicate field names")
    if set(fields) != set(EXPECTED_CLOUD_SCHEMA):
        errors.append(f"field names={sorted(fields)}")
    for name, expected in EXPECTED_CLOUD_SCHEMA.items():
        field = fields.get(name)
        if field is not None and (field.offset, field.datatype, field.count) != expected:
            errors.append(f"{name} schema={(field.offset, field.datatype, field.count)} expected={expected}")
    if message.height != 1:
        errors.append(f"height={message.height}")
    if message.point_step != 25:
        errors.append(f"point_step={message.point_step}")
    if message.row_step != message.point_step * message.width:
        errors.append(f"row_step={message.row_step}")
    if len(message.data) != message.row_step * message.height:
        errors.append(f"data_size={len(message.data)} expected={message.row_step * message.height}")
    return errors


def audit_cloud(message: Any) -> tuple[set[int], list[str]]:
    errors = cloud_schema_errors(message)
    if errors:
        return set(), errors
    endian = ">" if message.is_bigendian else "<"
    dtype = np.dtype(
        {
            "names": ["x", "y", "z", "intensity", "time", "lidar_id"],
            "formats": [endian + "f4", endian + "f4", endian + "f4", endian + "f4", endian + "f8", "u1"],
            "offsets": [0, 4, 8, 12, 16, 24],
            "itemsize": message.point_step,
        }
    )
    points = np.frombuffer(message.data, dtype=dtype, count=message.width)
    if len(points) < 1000:
        errors.append(f"point_count={len(points)}")
    for name in ("x", "y", "z", "intensity", "time"):
        if np.any(~np.isfinite(points[name])):
            errors.append(f"non-finite {name}")
    finite_times = points["time"][np.isfinite(points["time"])]
    if len(finite_times):
        if np.any((finite_times < -0.5) | (finite_times > 0.01)):
            errors.append("relative time outside [-0.5,0.01] s")
        if float(np.max(finite_times) - np.min(finite_times)) < 0.01:
            errors.append("relative-time span below 0.01 s")
    ids = points["lidar_id"].astype(np.int64)
    if np.any((ids < 0) | (ids > 3)):
        errors.append("lidar_id outside [0,3]")
    valid_ids = ids[(ids >= 0) & (ids <= 3)]
    counts = np.bincount(valid_ids, minlength=4)[:4]
    present = set(np.flatnonzero(counts > 0).tolist())
    if any(counts[lidar_id] < 100 for lidar_id in present):
        errors.append("a present lidar contributes fewer than 100 points")
    return present, errors


def ratio_in_window(samples: list[tuple[int, float]], begin: int, end: int, value: float) -> tuple[float, int]:
    selected = [sample for clock, sample in samples if begin <= clock < end]
    return ((sum(sample == value for sample in selected) / len(selected)) if selected else 0.0, len(selected))


def header_window(headers: list[int], begin: int, end: int) -> dict[str, float | int | None | bool]:
    selected = [stamp for stamp in headers if begin <= stamp < end]
    window_duration_s = max(0.0, (end - begin) * 1e-9)
    observed_span_s = (selected[-1] - selected[0]) * 1e-9 if len(selected) >= 2 else 0.0
    intervals = np.diff(np.asarray(selected, dtype=np.int64)).astype(np.float64) * 1e-9
    rate = (len(selected) - 1) / observed_span_s if observed_span_s > 0.0 else None
    max_gap = float(np.max(intervals)) if len(intervals) else None
    min_count = int(math.floor(window_duration_s * MIN_WINDOW_HEADER_RATE_HZ))
    passed = (
        window_duration_s > 0.0
        and len(selected) >= max(2, min_count)
        and observed_span_s >= MIN_WINDOW_SPAN_RATIO * window_duration_s
        and rate is not None and rate >= MIN_WINDOW_HEADER_RATE_HZ
        and max_gap is not None and max_gap <= MAX_CONTINUOUS_HEADER_GAP_S
    )
    return {
        "count": len(selected),
        "minimum_count": min_count,
        "window_duration_s": window_duration_s,
        "observed_span_s": observed_span_s,
        "header_rate_hz": rate,
        "max_header_gap_s": max_gap,
        "passed": passed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bag", type=Path, help="recorded frontend output bag")
    parser.add_argument("--fault-contract", type=Path, required=True)
    parser.add_argument("--input-bag", type=Path, required=True, help="actual fault-injected playback bag")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diagnostic-only", action="store_true")
    args = parser.parse_args()
    contract = json.loads(args.fault_contract.read_text(encoding="utf-8"))
    args.input_bag = args.input_bag.resolve()

    binding_failures: list[str] = []
    if contract.get("complete_input_pass") is not True:
        binding_failures.append("fault contract is not from a complete input pass")
    if Path(contract.get("output_bag", "")).resolve() != args.input_bag:
        binding_failures.append("fault contract output_bag does not match the actual playback bag")
    actual_bag_files = bag_file_identities(args.input_bag)
    if contract.get("bag_files") != actual_bag_files:
        binding_failures.append("fault contract bag file hashes do not match the actual playback bag")

    reader = rosbag2_py.SequentialReader()
    reader.open(rosbag2_py.StorageOptions(uri=str(args.bag), storage_id="sqlite3"), rosbag2_py.ConverterOptions("", ""))
    topic_types = {item.name: item.type for item in reader.get_all_topics_and_types()}
    missing = [topic for topic in EXPECTED if topic not in topic_types]
    if missing:
        raise SystemExit(f"recording lacks topics: {missing}")
    wrong_types = {
        topic: {"actual": topic_types[topic], "expected": EXPECTED_TYPES[topic]}
        for topic in EXPECTED if topic_types[topic] != EXPECTED_TYPES[topic]
    }
    if wrong_types:
        raise SystemExit(f"recording has incorrect topic types: {wrong_types}")
    classes = {topic: get_message(topic_types[topic]) for topic in EXPECTED}
    record_times = {topic: [] for topic in EXPECTED}
    pose_headers: list[int] = []
    cloud_frames: list[tuple[int, set[int]]] = []
    pose_frames: set[str] = set()
    cloud_frame_ids: set[str] = set()
    safety_values: list[float] = []
    state_by_clock: list[tuple[int, float]] = []
    system_by_clock: list[tuple[int, float]] = []
    state_by_record: list[tuple[int, float]] = []
    system_by_record: list[tuple[int, float]] = []
    latest_clock: int | None = None
    clocks: list[int] = []
    pose_nonfinite = 0
    quaternion_max_norm_error = 0.0
    first_pose_translation_norm: float | None = None
    first_pose_rotation_angle_deg: float | None = None
    cloud_invalid_frame_count = 0
    cloud_error_examples: list[str] = []

    while reader.has_next():
        topic, serialized, record_time = reader.read_next()
        if topic not in classes:
            continue
        message = deserialize_message(serialized, classes[topic])
        record_times[topic].append(int(record_time))
        if topic == CLOCK:
            latest_clock = stamp_ns(message.clock)
            clocks.append(latest_clock)
        elif topic == POSE:
            pose_headers.append(stamp_ns(message.header.stamp))
            pose_frames.add(message.header.frame_id)
            quaternion = message.pose.orientation
            values = (
                message.pose.position.x, message.pose.position.y, message.pose.position.z,
                quaternion.x, quaternion.y, quaternion.z, quaternion.w,
            )
            if not all(math.isfinite(value) for value in values):
                pose_nonfinite += 1
            else:
                norm = math.sqrt(sum(value * value for value in values[3:]))
                quaternion_max_norm_error = max(quaternion_max_norm_error, abs(norm - 1.0))
                if first_pose_translation_norm is None:
                    first_pose_translation_norm = math.sqrt(sum(value * value for value in values[:3]))
                    first_pose_rotation_angle_deg = (
                        math.degrees(2.0 * math.acos(min(1.0, abs(quaternion.w) / norm)))
                        if norm > 0.0 else math.inf
                    )
        elif topic == CLOUD:
            cloud_frame_ids.add(message.header.frame_id)
            ids, errors = audit_cloud(message)
            if errors:
                cloud_invalid_frame_count += 1
                if len(cloud_error_examples) < 5:
                    cloud_error_examples.append("; ".join(errors))
            cloud_frames.append((stamp_ns(message.header.stamp), ids))
        elif topic == SAFETY:
            safety_values.append(float(message.data))
        elif topic == STATE:
            state_by_record.append((int(record_time), float(message.data)))
            if latest_clock is not None:
                state_by_clock.append((latest_clock, float(message.data)))
        elif topic == SYSTEM:
            system_by_record.append((int(record_time), float(message.data)))
            if latest_clock is not None:
                system_by_clock.append((latest_clock, float(message.data)))

    failures = list(binding_failures)
    realtime_failures: list[str] = []
    for topic in (SAFETY, STATE, SYSTEM):
        rate = rate_hz(record_times[topic])
        if len(record_times[topic]) < 500 or rate is None or rate < MIN_WINDOW_HEADER_RATE_HZ or rate > 10.5:
            failures.append(f"{topic}: lifecycle stream is below the approved minimum 8 Hz throughput")
        elif rate < NOMINAL_HEADER_RATE_HZ:
            realtime_failures.append(f"{topic}: average receive rate is below the nominal 9.5 Hz realtime target")
    if len(safety_values) < 2 or not all(value in (0.0, 1.0) for value in safety_values) or not all(
        left != right for left, right in zip(safety_values, safety_values[1:])
    ):
        failures.append("safety heartbeat is not strict binary alternation")
    clock_deltas = [right - left for left, right in zip(clocks, clocks[1:])]
    clock_backward_count = sum(delta < 0 for delta in clock_deltas)
    clock_duplicate_count = sum(delta == 0 for delta in clock_deltas)
    maximum_clock_backward_jump_s = (
        max((-delta for delta in clock_deltas if delta < 0), default=0) * 1e-9
    )
    if len(clocks) < 2 or clock_backward_count:
        failures.append("clock is absent or moves backward")
    if pose_nonfinite or any(right <= left for left, right in zip(pose_headers, pose_headers[1:])):
        failures.append("pose is non-finite or nonmonotonic")
    if pose_frames and pose_frames != {"map"}:
        failures.append("pose frame_id contract failed")
    if quaternion_max_norm_error > 1e-3:
        failures.append("pose quaternion norm error exceeds 1e-3")
    cloud_headers = [stamp for stamp, _ in cloud_frames]
    if any(right <= left for left, right in zip(cloud_headers, cloud_headers[1:])):
        failures.append("cloud headers are nonmonotonic")
    if cloud_frame_ids and cloud_frame_ids != {"lidar_114"}:
        failures.append("cloud frame_id contract failed")
    if cloud_invalid_frame_count:
        failures.append("cloud schema/content contract failed")
    if pose_headers != cloud_headers:
        failures.append("pose/cloud headers are not exactly paired")

    scenario = contract["scenario"]
    mode = contract["mode"]
    first_sensor = int(contract["first_target_header_ns"])
    last_sensor = int(contract["last_target_header_ns"])
    start = first_sensor if mode == "drop_all" else int(contract["fault_start_header_ns"])
    end = last_sensor + 1 if mode == "drop_all" else int(contract["fault_end_header_ns"])
    details: dict[str, Any] = {
        "scenario": scenario,
        "fault_start_ns": start,
        "fault_end_ns": end,
        "input_binding": {
            "complete_input_pass": contract.get("complete_input_pass"),
            "actual_bag_files": actual_bag_files,
            "passed": not binding_failures,
        },
        "output_contract": {
            "pose_frame_ids": sorted(pose_frames),
            "cloud_frame_ids": sorted(cloud_frame_ids),
            "pose_nonfinite_count": pose_nonfinite,
            "quaternion_max_norm_error": quaternion_max_norm_error,
            "first_pose_translation_norm_m": first_pose_translation_norm,
            "first_pose_rotation_angle_deg": first_pose_rotation_angle_deg,
            "cloud_invalid_frame_count": cloud_invalid_frame_count,
            "cloud_error_examples": cloud_error_examples,
        },
        "window_contract_parameters": {
            "minimum_header_rate_hz": MIN_WINDOW_HEADER_RATE_HZ,
            "nominal_header_rate_hz": NOMINAL_HEADER_RATE_HZ,
            "minimum_window_span_ratio": MIN_WINDOW_SPAN_RATIO,
            "maximum_header_gap_s": MAX_CONTINUOUS_HEADER_GAP_S,
        },
        "clock_backward_count": clock_backward_count,
        "clock_duplicate_count": clock_duplicate_count,
        "maximum_clock_backward_jump_s": maximum_clock_backward_jump_s,
    }

    system_values = [value for _, value in system_by_record]
    state_values = [value for _, value in state_by_record]
    if not all(value in (0.0, 1.0) for value in state_values + system_values):
        failures.append("state topics contain non-binary values")
    if mode == "drop_all":
        if any(value != 0.0 for value in system_values) or any(value != 0.0 for value in state_values):
            failures.append("missing primary IMU from startup did not remain fail-safe/uninitialized")
        if pose_headers or cloud_headers:
            failures.append("missing primary IMU from startup produced pose/cloud outputs")
    else:
        system_transitions = sum(left != right for left, right in zip(system_values, system_values[1:]))
        if not system_values or system_values[-1] != 1.0 or system_transitions != 1:
            failures.append("SystemState did not initialize exactly once and remain initialized")
        first_state_one = next((time for time, value in state_by_record if value == 1.0), None)
        first_system_one = next((time for time, value in system_by_record if value == 1.0), None)
        lifecycle_delta_s = (
            (first_state_one - first_system_one) * 1e-9
            if first_state_one is not None and first_system_one is not None else None
        )
        details["slam_state_minus_system_state_first_one_receive_s"] = lifecycle_delta_s
        if lifecycle_delta_s is None or lifecycle_delta_s < -MAX_LIFECYCLE_REORDER_S:
            failures.append("slamState/SystemState initialization order exceeds 0.15 s receive reordering")
        if (
            first_pose_translation_norm is None or first_pose_translation_norm > 0.001
            or first_pose_rotation_angle_deg is None or first_pose_rotation_angle_deg > 0.05
        ):
            failures.append("first valid rear-axle pose is not the required identity origin")

        pre_begin, pre_end = start - 10_000_000_000, start - 1_000_000_000
        core_begin, core_end = start + 1_000_000_000, end - 1_000_000_000
        post_begin, post_end = end + 2_000_000_000, min(last_sensor, end + 12_000_000_000)
        pre_one, pre_n = ratio_in_window(state_by_clock, pre_begin, pre_end, 1.0)
        core_one, core_n = ratio_in_window(state_by_clock, core_begin, core_end, 1.0)
        post_one, post_n = ratio_in_window(state_by_clock, post_begin, post_end, 1.0)
        details["slam_state_windows"] = {
            "pre": {"normal_ratio": pre_one, "count": pre_n},
            "core": {"normal_ratio": core_one, "count": core_n},
            "post": {"normal_ratio": post_one, "count": post_n},
        }
        if pre_n < 20 or pre_one < 0.9 or post_n < 20 or post_one < 0.8:
            failures.append("tracking was not normal before the outage or did not recover afterward")

        pose_windows = {
            "pre": header_window(pose_headers, pre_begin, pre_end),
            "post": header_window(pose_headers, post_begin, post_end),
        }
        cloud_windows = {
            "pre": header_window(cloud_headers, pre_begin, pre_end),
            "post": header_window(cloud_headers, post_begin, post_end),
        }
        details["pose_header_windows"] = pose_windows
        details["cloud_header_windows"] = cloud_windows
        if not all(window["passed"] for window in (*pose_windows.values(), *cloud_windows.values())):
            failures.append("pose/cloud did not sustain the required pre/post header rate, span, count, and continuity")

        target_topic = contract["target_topic"]
        is_imu = "/imu_" in target_topic
        if is_imu:
            if core_n < 20 or core_one > 0.1:
                failures.append("primary IMU outage did not drive slamState to zero")
            core_outputs = sum(core_begin <= stamp < core_end for stamp in pose_headers)
            details["core_pose_count"] = core_outputs
            if core_outputs > 2:
                failures.append("primary IMU outage continued publishing fresh poses")
            first_post = next((stamp for stamp in pose_headers if stamp >= end), None)
            details["first_post_fault_pose_ns"] = first_post
            if first_post is None or first_post > end + 3_000_000_000:
                failures.append("primary IMU recovery did not resume pose output within 3 s")
        else:
            target_id = int(target_topic.rsplit("_", 1)[1])
            id_by_ip = {114: 0, 127: 1, 187: 2, 195: 3}
            missing_id = id_by_ip[target_id]
            expected_core = {0, 1, 2, 3} - {missing_id}
            core_frames = [ids for stamp, ids in cloud_frames if core_begin <= stamp < core_end]
            post_frames = [ids for stamp, ids in cloud_frames if post_begin <= stamp < post_end]
            core_exact = sum(ids == expected_core for ids in core_frames) / len(core_frames) if core_frames else 0.0
            post_full = sum(ids == {0, 1, 2, 3} for ids in post_frames) / len(post_frames) if post_frames else 0.0
            details["cloud_id_windows"] = {
                "core_count": len(core_frames), "core_exact_expected_ratio": core_exact,
                "post_count": len(post_frames), "post_full_ratio": post_full,
            }
            if core_n < 20 or core_one < 0.8:
                failures.append("degraded lidar mode did not sustain normal tracking")
            core_pose_window = header_window(pose_headers, core_begin, core_end)
            core_cloud_window = header_window(cloud_headers, core_begin, core_end)
            details["pose_header_windows"]["core"] = core_pose_window
            details["cloud_header_windows"]["core"] = core_cloud_window
            if not core_pose_window["passed"] or not core_cloud_window["passed"] or core_exact < 0.9:
                failures.append("degraded output did not sustain rate/continuity or exactly the remaining lidar IDs")
            if not cloud_windows["post"]["passed"] or post_full < 0.9:
                failures.append("four-lidar cloud composition did not recover over the full post-fault window")
            continuity_headers = [stamp for stamp in cloud_headers if pre_begin <= stamp < post_end]
            max_gap_s = (
                float(np.max(np.diff(np.asarray(continuity_headers, dtype=np.int64)))) * 1e-9
                if len(continuity_headers) >= 2 else None
            )
            details["nonprimary_full_window_max_cloud_gap_s"] = max_gap_s
            if max_gap_s is None or max_gap_s > MAX_CONTINUOUS_HEADER_GAP_S:
                failures.append("lidar outage broke continuous cloud output")

    result = {
        "bag": str(args.bag),
        "input_bag": str(args.input_bag),
        "fault_contract": str(args.fault_contract),
        "contract_passed": not failures,
        "contract_failures": failures,
        "realtime_contract_passed": not realtime_failures,
        "realtime_contract_failures": realtime_failures,
        "message_counts": {topic: len(record_times[topic]) for topic in EXPECTED},
        "receive_rates_hz": {topic: rate_hz(record_times[topic]) for topic in EXPECTED},
        "details": details,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"contract_passed": not failures, "failures": failures}, ensure_ascii=False))
    return 0 if args.diagnostic_only or not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
