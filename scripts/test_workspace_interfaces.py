#!/usr/bin/env python3
"""Regression checks for the two actual deployment interface failures."""
import contextlib
import io
from pathlib import Path
import tempfile
import unittest
import yaml
from check_workspace_repos import check_interfaces


class InterfaceTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.ws = Path(self.temp.name)
        self.repo = self.ws / "src/lightning-lm"
        self.loc = self.repo / "config/sany/localization.yaml"
        self.source = self.ws / "src/livox_ros_driver2/src/livox_ros_driver2.cpp"
        self.launch_dir = self.ws / "src/livox_ros_driver2/launch_ROS2"
        self.msg = self.ws / "src/common_msgs/diagnostic_monitor_interfaces/msg/NodeHeartbeat.msg"
        self.lidar_msg = self.ws / "src/common_msgs/livox_ros_driver/msg/LidarDiagnosticState.msg"
        for path in (self.loc, self.source, self.msg, self.lidar_msg):
            path.parent.mkdir(parents=True, exist_ok=True)
        self.launch_dir.mkdir(parents=True)
        self.loc.write_text(yaml.safe_dump({"diagnostics": {"radar_role_source_ids": [7, 12, 8, 9]}}))
        self.source.write_text('"livox-dtc-v3-bounded-sync"; dtc::PublishDiagnosticSnapshot();')
        for role, node_id in zip(("front", "rear", "left", "right"), (7, 12, 8, 9)):
            (self.launch_dir / f"lidar_{role}.launch.py").write_text(
                f"heartbeat_node_id = {node_id}\nheartbeat_topic = '/diagnostics/heartbeat/livox_{role}'\n")
        (self.launch_dir / "master_lidar.launch.py").write_text("# Four source drivers only\n")
        self.msg.write_text("uint32[] active_fault\n")
        self.lidar_msg.write_text("bool bounded_time_sync_warning\nint64 time_sync_offset_ns\nuint64 time_sync_evidence_age_us\n")
        self.identity_msg = self.lidar_msg.with_name("LidarDeviceInfo.msg")
        self.identity_msg.write_text("std_msgs/Header header\nstring role\nstring ip\nstring serial_number\nuint8 device_type\nbool identity_valid\n")

    def result(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return check_interfaces(self.repo, self.ws)

    def test_compatible(self):
        self.assertEqual(self.result(), 0)

    def test_missing_identity_interface(self):
        self.identity_msg.unlink()
        self.assertEqual(self.result(), 1)

    def test_incomplete_identity_interface(self):
        self.identity_msg.write_text("string serial_number\n")
        self.assertEqual(self.result(), 1)

    def test_old_bounded_proof_schema(self):
        self.lidar_msg.write_text("uint8 usability\n")
        self.assertEqual(self.result(), 1)

    def test_old_driver_semantics(self):
        self.source.write_text('"livox-dtc-v2-sync-source"; dtc::PublishDiagnosticSnapshot();')
        self.assertEqual(self.result(), 1)

    def test_driver_without_own_dtc_reporting(self):
        self.source.write_text('"livox-dtc-v3-bounded-sync";')
        self.assertEqual(self.result(), 1)

    def test_old_message(self):
        self.msg.write_text("uint8 state\n")
        self.assertEqual(self.result(), 1)

    def test_unregistered_aggregate_node(self):
        (self.launch_dir / "master_lidar.launch.py").write_text("lidar_dtc_aggregator.launch.py")
        self.assertEqual(self.result(), 1)

    def test_wrong_source(self):
        p = self.launch_dir / "lidar_front.launch.py"
        p.write_text(p.read_text().replace("= 7", "= 12"))
        self.assertEqual(self.result(), 1)

    def test_wrong_role_order(self):
        self.loc.write_text(yaml.safe_dump({"diagnostics": {"radar_role_source_ids": [7, 8, 9, 12]}}))
        self.assertEqual(self.result(), 1)

    def test_retired_aggregate_dependency(self):
        self.loc.write_text(yaml.safe_dump({"diagnostics": {
            "radar_role_source_ids": [7, 12, 8, 9], "radar_source_id": 13}}))
        self.assertEqual(self.result(), 1)


if __name__ == "__main__":
    unittest.main()
