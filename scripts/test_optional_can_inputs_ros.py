#!/usr/bin/env python3
"""Exercise the real startup probe with received LiDAR/IMU and absent CAN.

Source the ROS/workspace environment, then pass check_localization_inputs.
Uses localhost-only ROS domain 173; publishes no vehicle control messages.
"""
import os
os.environ["ROS_LOCALHOST_ONLY"] = "1"
os.environ["ROS_DOMAIN_ID"] = "173"

import subprocess
import sys
import time
from pathlib import Path
import rclpy
from rclpy.qos import qos_profile_sensor_data
from sensor_msgs.msg import Imu, PointCloud2


def main():
    executable = str(Path(sys.argv[1]).resolve())
    rclpy.init()
    node = rclpy.create_node("optional_can_probe_fixture")
    lidar = node.create_publisher(PointCloud2, "/test/optional_can/lidar", qos_profile_sensor_data)
    imu = node.create_publisher(Imu, "/test/optional_can/imu", qos_profile_sensor_data)
    try:
        for option, send_imu, send_lidar, expected in (
            ("--optional-wheel-speed", True, True, 0),
            ("--wheel-speed", True, True, 1),
            ("--optional-wheel-speed", False, True, 1),
            ("--optional-wheel-speed", True, False, 1),
        ):
            command = [executable, "--timeout", "3", "--lidar", "/test/optional_can/lidar",
                       "--imu", "/test/optional_can/imu", option, "/test/optional_can/can"]
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            try:
                deadline = time.monotonic() + 10
                while process.poll() is None and time.monotonic() < deadline:
                    stamp = node.get_clock().now().to_msg()
                    if send_lidar:
                        cloud = PointCloud2()
                        cloud.header.stamp = stamp
                        cloud.width = cloud.height = 1
                        cloud.data = [0]
                        lidar.publish(cloud)
                    if send_imu:
                        msg = Imu()
                        msg.header.stamp = stamp
                        imu.publish(msg)
                    rclpy.spin_once(node, timeout_sec=0.02)
                output, _ = process.communicate(timeout=2)
                assert process.returncode == expected, (command, process.returncode, output)
                assert "messages=" in output, output
                print(f"PASS {option}, IMU={send_imu}, LiDAR={send_lidar}, exit={process.returncode}")
            finally:
                if process.poll() is None:
                    process.kill()
                    process.communicate()
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == "__main__":
    main()
