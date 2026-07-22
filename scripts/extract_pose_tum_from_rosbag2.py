#!/usr/bin/env python3
"""Extract a PoseStamped or Odometry topic from a ROS 2 bag as a TUM trajectory."""

import argparse
from pathlib import Path

import rosbag2_py
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bag", required=True, help="ROS 2 bag directory")
    parser.add_argument("--topic", default="/slamPoseRaw_topic")
    parser.add_argument("--output", required=True, help="output TUM file")
    return parser.parse_args()


def main():
    args = parse_args()
    storage_options = rosbag2_py.StorageOptions(uri=args.bag, storage_id="sqlite3")
    reader = rosbag2_py.SequentialReader()
    reader.open(storage_options, rosbag2_py.ConverterOptions("", ""))

    topic_types = {entry.name: entry.type for entry in reader.get_all_topics_and_types()}
    if args.topic not in topic_types:
        raise RuntimeError(f"topic not found in bag: {args.topic}")
    message_type = get_message(topic_types[args.topic])

    samples = []
    while reader.has_next():
        topic, serialized, _ = reader.read_next()
        if topic != args.topic:
            continue
        message = deserialize_message(serialized, message_type)
        stamp = message.header.stamp
        timestamp = stamp.sec + stamp.nanosec * 1e-9
        pose = message.pose.pose if hasattr(message.pose, "pose") else message.pose
        samples.append(
            (
                timestamp,
                pose.position.x,
                pose.position.y,
                pose.position.z,
                pose.orientation.x,
                pose.orientation.y,
                pose.orientation.z,
                pose.orientation.w,
            )
        )

    samples.sort(key=lambda sample: sample[0])
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        last_timestamp = None
        for sample in samples:
            if last_timestamp is not None and sample[0] <= last_timestamp:
                continue
            stream.write(f"{sample[0]:.9f} " + " ".join(f"{value:.12f}" for value in sample[1:]) + "\n")
            last_timestamp = sample[0]

    print(f"topic={args.topic} samples={len(samples)} output={output}")
    return 0 if samples else 2


if __name__ == "__main__":
    raise SystemExit(main())
