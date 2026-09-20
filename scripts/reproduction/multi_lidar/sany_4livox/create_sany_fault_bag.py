#!/usr/bin/env python3
"""Create one auditable ROS 2 bag with a deterministic topic outage."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import rosbag2_py
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stamp_ns(message: Any) -> int:
    return int(message.header.stamp.sec) * 1_000_000_000 + int(message.header.stamp.nanosec)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--target-topic", required=True)
    parser.add_argument("--mode", choices=("drop_all", "interrupt"), required=True)
    parser.add_argument("--start-offset-s", type=float, default=40.0)
    parser.add_argument("--duration-s", type=float, default=20.0)
    parser.add_argument("--max-messages", type=int, default=0, help="test-only input message limit")
    args = parser.parse_args()

    args.input = args.input.resolve()
    args.output = args.output.resolve()
    if not (args.input / "metadata.yaml").is_file():
        raise SystemExit(f"input is not a ROS 2 bag: {args.input}")
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite: {args.output}")
    if args.start_offset_s < 0.0 or args.duration_s <= 0.0:
        raise SystemExit("fault timing must be positive")

    reader = rosbag2_py.SequentialReader()
    reader.open(
        rosbag2_py.StorageOptions(uri=str(args.input), storage_id="sqlite3"),
        rosbag2_py.ConverterOptions("", ""),
    )
    topics = reader.get_all_topics_and_types()
    topic_types = {topic.name: topic.type for topic in topics}
    if args.target_topic not in topic_types:
        raise SystemExit(f"target topic is absent: {args.target_topic}")
    target_class = get_message(topic_types[args.target_topic])

    writer = rosbag2_py.SequentialWriter()
    writer.open(
        rosbag2_py.StorageOptions(uri=str(args.output), storage_id="sqlite3"),
        rosbag2_py.ConverterOptions("", ""),
    )
    for topic in topics:
        writer.create_topic(topic)

    input_counts = {topic.name: 0 for topic in topics}
    output_counts = {topic.name: 0 for topic in topics}
    dropped_count = 0
    first_target_header_ns: int | None = None
    last_target_header_ns: int | None = None
    fault_start_ns: int | None = None
    fault_end_ns: int | None = None
    total = 0
    while reader.has_next():
        topic, serialized, record_time = reader.read_next()
        input_counts[topic] = input_counts.get(topic, 0) + 1
        drop = False
        if topic == args.target_topic:
            header_ns = stamp_ns(deserialize_message(serialized, target_class))
            if first_target_header_ns is None:
                first_target_header_ns = header_ns
                fault_start_ns = header_ns + int(round(args.start_offset_s * 1e9))
                fault_end_ns = fault_start_ns + int(round(args.duration_s * 1e9))
            last_target_header_ns = header_ns
            drop = args.mode == "drop_all" or bool(fault_start_ns <= header_ns < fault_end_ns)
        if drop:
            dropped_count += 1
        else:
            writer.write(topic, serialized, int(record_time))
            output_counts[topic] = output_counts.get(topic, 0) + 1
        total += 1
        if args.max_messages and total >= args.max_messages:
            break

    del writer
    if first_target_header_ns is None:
        raise SystemExit("target topic contained no messages")
    if args.mode == "interrupt" and dropped_count == 0:
        raise SystemExit("fault interval removed no target messages")
    if args.mode == "drop_all" and output_counts[args.target_topic] != 0:
        raise SystemExit("drop_all left target messages in output")

    bag_files = sorted(
        path for path in args.output.iterdir()
        if path.is_file() and (path.name == "metadata.yaml" or path.suffix in {".db3", ".mcap"})
    )
    contract = {
        "created_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "scenario": args.scenario,
        "source_bag": str(args.input),
        "output_bag": str(args.output),
        "target_topic": args.target_topic,
        "mode": args.mode,
        "start_offset_s": args.start_offset_s,
        "duration_s": args.duration_s,
        "fault_start_header_ns": fault_start_ns,
        "fault_end_header_ns": fault_end_ns,
        "first_target_header_ns": first_target_header_ns,
        "last_target_header_ns": last_target_header_ns,
        "dropped_message_count": dropped_count,
        "input_counts": input_counts,
        "output_counts": output_counts,
        "complete_input_pass": args.max_messages == 0,
        "bag_files": [
            {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)} for path in bag_files
        ],
    }
    (args.output / "fault_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"scenario": args.scenario, "dropped": dropped_count, "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
