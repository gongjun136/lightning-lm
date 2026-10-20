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
        self.config = self.ws / "src/livox_ros_driver2/config/dtc_aggregator.yaml"
        self.source = self.ws / "src/livox_ros_driver2/src/lidar_dtc_aggregator.cpp"
        self.msg = self.ws / "src/common_msgs/diagnostic_monitor_interfaces/msg/NodeHeartbeat.msg"
        for path in (self.loc, self.config, self.source, self.msg):
            path.parent.mkdir(parents=True, exist_ok=True)
        self.loc.write_text(yaml.safe_dump({"diagnostics": {
            "internal_radar_topic": "/internal/diagnostics/livox_dtc",
            "radar_source_id": 13, "radar_role_source_ids": [7, 12, 8, 9]}}))
        self.params = {"internal_output_topic": "/internal/diagnostics/livox_dtc",
                       "functional_safety.node_id": 13, "expected_source_ids": [7, 12, 8, 9]}
        self.write_config()
        self.source.write_text('declare_parameter<std::string>("internal_output_topic", "default");')
        self.msg.write_text("uint32[] active_fault\n")

    def write_config(self):
        self.config.write_text(yaml.safe_dump({"lidar_dtc_aggregator": {"ros__parameters": self.params}}))

    def result(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return check_interfaces(self.repo, self.ws)

    def test_compatible(self):
        self.assertEqual(self.result(), 0)

    def test_old_driver_topic(self):
        self.params["output_topic"] = "/diagnostics/source_dtc/livox"
        del self.params["internal_output_topic"]
        self.write_config()
        self.assertEqual(self.result(), 1)

    def test_old_message(self):
        self.msg.write_text("uint8 state\n")
        self.assertEqual(self.result(), 1)

    def test_old_producer_ignores_new_yaml(self):
        self.source.write_text('declare_parameter<std::string>("output_topic", "old");')
        self.assertEqual(self.result(), 1)

    def test_wrong_source(self):
        self.params["functional_safety.node_id"] = 12
        self.write_config()
        self.assertEqual(self.result(), 1)

    def test_wrong_role_order(self):
        self.params["expected_source_ids"] = [7, 8, 9, 12]
        self.write_config()
        self.assertEqual(self.result(), 1)


if __name__ == "__main__":
    unittest.main()
