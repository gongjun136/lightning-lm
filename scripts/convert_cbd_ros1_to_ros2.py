#!/usr/bin/env python3
"""Convert CBD_Building_03 ROS1 Livox bag to a ROS2 sqlite3 bag.

Only the topics needed by lightning-lm are converted:
  - /livox/lidar: livox_ros_driver/msg/CustomMsg -> livox_ros_driver2/msg/CustomMsg
  - /livox/imu: sensor_msgs/msg/Imu -> sensor_msgs/msg/Imu
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from rosbags.rosbag1 import Reader
from rosbags.rosbag2 import Writer
from rosbags.typesys import Stores, get_types_from_msg, get_typestore


LIDAR_TOPIC = "/livox/lidar"
IMU_TOPIC = "/livox/imu"
SRC_LIVOX_TYPE = "livox_ros_driver/msg/CustomMsg"
DST_LIVOX_TYPE = "livox_ros_driver2/msg/CustomMsg"
IMU_TYPE = "sensor_msgs/msg/Imu"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", required=True, type=Path, help="ROS1 .bag input path")
    parser.add_argument("--dst", required=True, type=Path, help="ROS2 bag output directory")
    parser.add_argument("--force", action="store_true", help="delete existing output directory first")
    parser.add_argument(
        "--max-messages",
        type=int,
        default=0,
        help="convert at most this many messages; useful for smoke tests",
    )
    parser.add_argument("--progress-every", type=int, default=5000, help="progress print interval")
    return parser.parse_args()


def register_livox_types(src_store, dst_store, reader: Reader, repo_root: Path) -> None:
    lidar_conn = next((conn for conn in reader.connections if conn.topic == LIDAR_TOPIC), None)
    if lidar_conn is None:
        raise RuntimeError(f"missing required topic: {LIDAR_TOPIC}")

    src_types = get_types_from_msg(lidar_conn.msgdef.data, lidar_conn.msgtype)
    src_store.register({name: spec for name, spec in src_types.items() if name.startswith("livox_ros_driver/")})

    msg_dir = repo_root / "thirdparty" / "livox_ros_driver" / "msg"
    for msg_name in ("CustomPoint", "CustomMsg"):
        msg_text = (msg_dir / f"{msg_name}.msg").read_text(encoding="utf-8")
        dst_store.register(get_types_from_msg(msg_text, f"livox_ros_driver2/msg/{msg_name}"))


def ros2_header(src_header, dst_store):
    time_cls = dst_store.types["builtin_interfaces/msg/Time"]
    header_cls = dst_store.types["std_msgs/msg/Header"]
    stamp = time_cls(src_header.stamp.sec, src_header.stamp.nanosec)
    return header_cls(stamp, src_header.frame_id)


def convert_livox(src_msg, dst_store):
    dst_msg_cls = dst_store.types[DST_LIVOX_TYPE]
    return dst_msg_cls(
        ros2_header(src_msg.header, dst_store),
        src_msg.timebase,
        src_msg.point_num,
        src_msg.lidar_id,
        src_msg.rsvd,
        src_msg.points,
    )


def convert_imu(src_msg, dst_store):
    src_msg.header = ros2_header(src_msg.header, dst_store)
    return src_msg


def prepare_output(path: Path, force: bool) -> None:
    if not path.exists():
        return
    if not force:
        raise RuntimeError(f"output path already exists: {path}; pass --force to replace it")
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def patch_humble_metadata(path: Path) -> None:
    """Make rosbags metadata readable by ROS2 Humble rosbag2 tools."""
    metadata = path / "metadata.yaml"
    text = metadata.read_text(encoding="utf-8")
    metadata.write_text(text.replace("offered_qos_profiles: []", 'offered_qos_profiles: ""'), encoding="utf-8")


def main() -> None:
    args = parse_args()
    if not args.src.is_file():
        raise RuntimeError(f"input bag does not exist: {args.src}")

    repo_root = Path(__file__).resolve().parents[1]
    prepare_output(args.dst, args.force)

    src_store = get_typestore(Stores.ROS1_NOETIC)
    dst_store = get_typestore(Stores.ROS2_HUMBLE)

    converted = {LIDAR_TOPIC: 0, IMU_TOPIC: 0}
    with Reader(args.src) as reader:
        register_livox_types(src_store, dst_store, reader, repo_root)
        input_connections = [conn for conn in reader.connections if conn.topic in converted]
        topics = {conn.topic for conn in input_connections}
        missing = sorted(set(converted) - topics)
        if missing:
            raise RuntimeError(f"missing required topics: {', '.join(missing)}")

        with Writer(args.dst, version=9) as writer:
            output_connections = {
                IMU_TOPIC: writer.add_connection(IMU_TOPIC, IMU_TYPE, typestore=dst_store),
                LIDAR_TOPIC: writer.add_connection(LIDAR_TOPIC, DST_LIVOX_TYPE, typestore=dst_store),
            }

            total = 0
            for conn, timestamp, raw in reader.messages(connections=input_connections):
                if conn.topic == IMU_TOPIC:
                    src_msg = src_store.deserialize_ros1(raw, IMU_TYPE)
                    dst_msg = convert_imu(src_msg, dst_store)
                    data = dst_store.serialize_cdr(dst_msg, IMU_TYPE)
                elif conn.topic == LIDAR_TOPIC:
                    src_msg = src_store.deserialize_ros1(raw, SRC_LIVOX_TYPE)
                    dst_msg = convert_livox(src_msg, dst_store)
                    data = dst_store.serialize_cdr(dst_msg, DST_LIVOX_TYPE)
                else:
                    continue

                writer.write(output_connections[conn.topic], timestamp, data)
                converted[conn.topic] += 1
                total += 1

                if args.progress_every > 0 and total % args.progress_every == 0:
                    print(f"converted {total} messages")
                if args.max_messages > 0 and total >= args.max_messages:
                    break

    patch_humble_metadata(args.dst)
    print(f"done: {args.dst}")
    print(f"{IMU_TOPIC}: {converted[IMU_TOPIC]}")
    print(f"{LIDAR_TOPIC}: {converted[LIDAR_TOPIC]}")


if __name__ == "__main__":
    main()
