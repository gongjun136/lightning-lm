#!/usr/bin/env python3
"""Inspect a SQLite3 or MCAP ROS2 bag and derive the offline completion contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def primary_lidar_topic(config: dict) -> str:
    common = config.get("common") or {}
    multi = config.get("multi_lidar") or {}
    if bool(multi.get("enabled", False)):
        primary_id = int(multi.get("primary_lidar_id", 0))
        topic = (multi.get("topics") or {}).get(f"lidar_{primary_id}", "")
        if topic:
            return str(topic)
    return str(common.get("lidar_topic") or common.get("livox_lidar_topic") or "")


def primary_imu_topic(config: dict) -> str:
    common = config.get("common") or {}
    multi = config.get("multi_lidar") or {}
    if bool(multi.get("enabled", False)):
        primary_id = int(multi.get("primary_lidar_id", 0))
        topic = (multi.get("topics") or {}).get(f"imu_{primary_id}", "")
        if topic:
            return str(topic)
    return str(common.get("imu_topic") or "")


def topic_aliases(topic: str) -> list[str]:
    aliases = [topic]
    alternate = topic[1:] if topic.startswith("/") else f"/{topic}"
    if alternate and alternate not in aliases:
        aliases.append(alternate)
    return aliases


def inventory_contract(path: Path, sequence: str) -> tuple[float, float] | None:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    rows = payload if isinstance(payload, list) else payload.get("sequences", payload.get("items", []))
    row = next((item for item in rows if item.get("sequence") == sequence), None)
    if not row:
        raise SystemExit(f"sequence is absent from inventory: {sequence}")
    return float(row["lidar_last_end_ns"]) / 1e9, float(row["sensor_duration_s"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bag", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inventory", type=Path)
    parser.add_argument("--sequence", default="")
    args = parser.parse_args()

    metadata_path = args.bag / "metadata.yaml"
    metadata = yaml.safe_load(metadata_path.read_text(encoding="utf-8-sig"))
    info = metadata.get("rosbag2_bagfile_information") or {}
    storage_identifier = str(info.get("storage_identifier") or "")
    if storage_identifier not in {"sqlite3", "mcap"}:
        raise SystemExit(f"unsupported ROS2 bag storage: {storage_identifier!r}")
    relative_paths = list(info.get("relative_file_paths") or [])
    if not relative_paths:
        relative_paths = [item.get("path") for item in (info.get("files") or []) if item.get("path")]
    suffix = ".db3" if storage_identifier == "sqlite3" else ".mcap"
    payloads = [args.bag / path for path in relative_paths] or sorted(args.bag.glob(f"*{suffix}"))
    if not payloads or any(not path.is_file() or path.suffix != suffix for path in payloads):
        raise SystemExit(f"metadata does not identify valid {suffix} payloads")

    config = yaml.safe_load(args.config.read_text(encoding="utf-8-sig")) or {}
    configured_topic = primary_lidar_topic(config)
    configured_imu_topic = primary_imu_topic(config)
    if not configured_topic:
        raise SystemExit("configuration does not define a primary lidar topic")
    candidates = topic_aliases(configured_topic)
    imu_candidates = topic_aliases(configured_imu_topic) if configured_imu_topic else []

    global_min = global_max = topic_min = topic_max = imu_min = imu_max = None
    topic_count = 0
    imu_count = 0
    resolved_topic = configured_topic
    resolved_imu_topic = configured_imu_topic
    if storage_identifier == "sqlite3":
        for payload in payloads:
            connection = sqlite3.connect(f"file:{payload.as_posix()}?mode=ro", uri=True)
            try:
                row = connection.execute("SELECT MIN(timestamp), MAX(timestamp) FROM messages").fetchone()
                if row and row[0] is not None:
                    global_min = row[0] if global_min is None else min(global_min, row[0])
                    global_max = row[1] if global_max is None else max(global_max, row[1])
                for candidate in candidates:
                    row = connection.execute(
                        "SELECT MIN(messages.timestamp), MAX(messages.timestamp), COUNT(*) "
                        "FROM messages JOIN topics ON messages.topic_id=topics.id WHERE topics.name=?",
                        (candidate,),
                    ).fetchone()
                    if row and row[2]:
                        resolved_topic = candidate
                        topic_count += int(row[2])
                        topic_min = row[0] if topic_min is None else min(topic_min, row[0])
                        topic_max = row[1] if topic_max is None else max(topic_max, row[1])
                        break
                for candidate in imu_candidates:
                    row = connection.execute(
                        "SELECT MIN(messages.timestamp), MAX(messages.timestamp), COUNT(*) "
                        "FROM messages JOIN topics ON messages.topic_id=topics.id WHERE topics.name=?",
                        (candidate,),
                    ).fetchone()
                    if row and row[2]:
                        resolved_imu_topic = candidate
                        imu_count += int(row[2])
                        imu_min = row[0] if imu_min is None else min(imu_min, row[0])
                        imu_max = row[1] if imu_max is None else max(imu_max, row[1])
                        break
            finally:
                connection.close()
    else:
        global_min = int((info.get("starting_time") or {}).get("nanoseconds_since_epoch") or 0)
        duration_ns = int((info.get("duration") or {}).get("nanoseconds") or 0)
        if global_min <= 0 or duration_ns <= 0:
            raise SystemExit("MCAP metadata is missing a valid start time or duration")
        global_max = global_min + duration_ns
        topics = {
            str((item.get("topic_metadata") or {}).get("name") or ""): int(item.get("message_count") or 0)
            for item in (info.get("topics_with_message_count") or [])
        }
        for candidate in candidates:
            if topics.get(candidate, 0) > 0:
                resolved_topic = candidate
                topic_count = topics[candidate]
                topic_min, topic_max = global_min, global_max
                break
        for candidate in imu_candidates:
            if topics.get(candidate, 0) > 0:
                resolved_imu_topic = candidate
                imu_count = topics[candidate]
                imu_min, imu_max = global_min, global_max
                break
    if global_min is None or global_max is None:
        raise SystemExit("bag contains no messages")

    duration_ns = int((info.get("duration") or {}).get("nanoseconds") or (global_max - global_min))
    if topic_max is None:
        expected_end_ns = global_max
        source = "bag_end_fallback"
    elif imu_max is not None:
        # A lidar frame cannot be propagated past the last primary IMU sample.
        # Using only the final lidar timestamp makes complete multi-sensor runs
        # fail their contract when one recorder stops a little later than the
        # other, even though every processable frame was consumed.
        expected_end_ns = min(topic_max, imu_max)
        source = (
            "primary_lidar_and_imu_sqlite"
            if storage_identifier == "sqlite3"
            else "bag_metadata_mcap"
        )
    else:
        expected_end_ns = topic_max
        source = "primary_lidar_sqlite" if storage_identifier == "sqlite3" else "bag_metadata_mcap"
    sensor_duration_s = duration_ns / 1e9
    if args.inventory:
        expected_s, sensor_duration_s = inventory_contract(args.inventory, args.sequence)
        expected_end_ns = round(expected_s * 1e9)
        source = "inventory_json"

    payload_records = [
        {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)} for path in payloads
    ]
    aggregate = hashlib.sha256(
        "".join(f"{item['sha256']}  {Path(item['path']).name}\n" for item in payload_records).encode()
    ).hexdigest()
    result = {
        "schema_version": 1,
        "bag": str(args.bag.resolve()),
        "metadata_sha256": sha256(metadata_path),
        "storage_identifier": storage_identifier,
        "primary_lidar_topic_configured": configured_topic,
        "primary_lidar_topic": resolved_topic,
        "primary_lidar_message_count": topic_count,
        "primary_lidar_first_s": topic_min / 1e9 if topic_min is not None else None,
        "primary_lidar_last_s": topic_max / 1e9 if topic_max is not None else None,
        "primary_imu_topic_configured": configured_imu_topic,
        "primary_imu_topic": resolved_imu_topic,
        "primary_imu_message_count": imu_count,
        "primary_imu_first_s": imu_min / 1e9 if imu_min is not None else None,
        "primary_imu_last_s": imu_max / 1e9 if imu_max is not None else None,
        "bag_first_s": global_min / 1e9,
        "bag_last_s": global_max / 1e9,
        "expected_last_lidar_s": expected_end_ns / 1e9,
        "expected_end_source": source,
        "sensor_duration_s": sensor_duration_s,
        "payload_file_count": len(payload_records),
        "payload_sha256": aggregate,
        "payloads": payload_records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("primary_lidar_topic", "expected_end_source", "expected_last_lidar_s", "sensor_duration_s")}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
