#!/usr/bin/env python3
"""Synthetic contract/geometry tests; optional real rosbag sqlite3/MCAP smoke test.

python3 scripts/test_accuracy_export_bag.py
python3 scripts/test_accuracy_export_bag.py --ros-smoke NEW_DIR --storage mcap
python3 scripts/test_accuracy_export_bag.py --ros-live-smoke NEW_DIR --domain-id 221
python3 scripts/test_accuracy_export_bag.py --core-smoke NEW_DIR
The second command requires sourced ROS 2 Humble and cgi430_interfaces overlay.
It only writes synthetic bags/results under NEW_DIR; no graph or hardware access.
"""
import argparse
import copy
import csv
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np
import yaml

import accuracy_export_bag as ex


def config():
    c = yaml.safe_load((Path(__file__).resolve().parents[1] / "config/accuracy_capture.yaml").read_text(encoding="utf-8"))
    c.update(source_recording_id="synthetic-only-not-field-evidence", epoch="epoch_1",
             time_basis="unix_utc_receive_proxy", world_frame_id="synthetic_enu_0_0_0")
    for key in ("georeference", "height", "reference_point", "clock", "assembly", "quality", "estimate", "capture"):
        c[key].update(confirmed=True, evidence="Synthetic fixture only; not vehicle calibration")
    c["georeference"]["origin_llh"] = [0, 0, 0]
    c["height"].update(datum="ellipsoidal", geoid_separation_m=0)
    c["clock"].update(reference_header_offset_ns=0, reference_bag_offset_ns=0,
                       estimate_header_offset_ns=0, estimate_bag_offset_ns=0,
                       uncertainty_ns=1_000_000, latency_evidence="Synthetic clock")
    # Intentionally synthetic gate values, never installed into field config.
    c["quality"].update(max_age_sec=.1, max_skew_sec=.02, max_future_sec=.001,
                         max_differential_age_sec=1.3, max_horizontal_sigma_m=.08,
                         max_vertical_sigma_m=.17, max_attitude_sigma_deg=.4,
                         max_velocity_sigma_mps=.12, max_site_radius_m=1000., recovery_samples=1)
    c["estimate"]["world_frame_id"] = "map:SYNTHETIC_SAME_ASSET"
    c["capture"]["expected_can_bitrate"] = 500000
    c["capture"]["input_topics"] = ["/synthetic/lidar", "/synthetic/imu"]
    c["capture"]["required_output_topics"] = ["/slamPoseRaw_topic"]
    return c


def navigation_cycle(base_ns=1_700_000_000_000_000_000, overrides=None):
    vals = {
        "latitude": {"latitude_deg": 0.}, "longitude": {"longitude_deg": 0.},
        "altitude": {"altitude_m": 10.},
        "attitude": {"heading_deg": 90., "pitch_deg": 0., "roll_deg": 0.},
        "velocity": {"east_mps": 0., "north_mps": 0., "up_mps": 0.},
        "position_sigma": {"east_m": .01, "north_m": .01, "up_m": .02},
        "attitude_sigma": {"heading_deg": .03, "pitch_deg": .02, "roll_deg": .02},
        "status": {"system_state": 2, "satellite_status": 4, "differential_age_s": .1},
        "velocity_sigma": {"east_mps": .01, "north_mps": .01, "up_mps": .01},
    }
    for field, changes in (overrides or {}).items():
        vals[field].update(changes)
    records = []
    for index, field in enumerate([k for k in ex.FIELDS if k != "latitude"] + ["latitude"]):
        stamp = base_ns + index * 100_000
        frame = {"header": {"stamp": {"sec": stamp // ex.NS, "nanosec": stamp % ex.NS}, "frame_id": "cgi430_link"},
                 "valid": True, "can_id": 800 + index, "is_extended": False, "dlc": 8, "data": [0] * 8}
        message = dict(vals[field], frame=frame)
        records.append(ex.Record("/cgi430/" + ex.FIELDS[field][0], "cgi430_interfaces/msg/" + ex.FIELDS[field][1], stamp + 1_000_000, message, 80))
    return records


def one_row(records, c=None):
    assembly = ex.ReferenceAssembler(c or config())
    rows = [row for record in records if (row := assembly.ingest(record)) is not None]
    return rows[-1]


def quat_matrix(q):
    x, y, z, w = q
    return np.array([[1 - 2 * (y*y + z*z), 2 * (x*y - z*w), 2 * (x*z + y*w)],
                     [2 * (x*y + z*w), 1 - 2 * (x*x + z*z), 2 * (y*z - x*w)],
                     [2 * (x*z - y*w), 2 * (y*z + x*w), 1 - 2 * (x*x + y*y)]])


class AccuracyExportTests(unittest.TestCase):
    def test_wgs84_analytic_and_fixed_basis(self):
        p, q = ex.convert_pose([0, 0, 10], [90, 0, 0], [0, 0, 0])
        np.testing.assert_allclose(p, [0, 0, 10], atol=1e-9)
        np.testing.assert_allclose(quat_matrix(q), np.eye(3), atol=1e-12)
        angle, a = 1e-5, 6378137.
        p, _ = ex.convert_pose([0, math.degrees(angle), 0], [90, 0, 0], [0, 0, 0])
        np.testing.assert_allclose(p, [a * math.sin(angle), 0, a * (math.cos(angle)-1)], atol=1e-8)
        _, q = ex.convert_pose([0, 90, 0], [90, 0, 0], [0, 0, 0])
        np.testing.assert_allclose(quat_matrix(q) @ [1, 0, 0], [0, 0, -1], atol=1e-12)
        np.testing.assert_allclose(quat_matrix(q) @ [0, 0, 1], [1, 0, 0], atol=1e-12)

    def test_attitude_cardinals_pitch_roll_and_pi(self):
        for heading, forward in [(0, [0, 1, 0]), (90, [1, 0, 0]), (180, [0, -1, 0]), (270, [-1, 0, 0])]:
            r = ex.vehicle_attitude(heading, 0, 0)
            np.testing.assert_allclose(r @ [1, 0, 0], forward, atol=1e-12)
            np.testing.assert_allclose(quat_matrix(ex.quaternion_xyzw(r)), r, atol=1e-12)
        self.assertGreater((ex.vehicle_attitude(90, 20, 0) @ [1, 0, 0])[2], 0)
        self.assertGreater((ex.vehicle_attitude(90, 0, 20) @ [0, 1, 0])[2], 0)

    def test_good_sample_and_no_second_lever_arm(self):
        row = one_row(navigation_cycle())
        self.assertEqual(row["reference_valid"], 1)
        self.assertEqual(row["reason"], "good")
        np.testing.assert_allclose([row[k] for k in "xyz"], [0, 0, 10], atol=1e-9)
        self.assertEqual(row["timestamp_ns"], 1_700_000_000_000_800_000)
        self.assertEqual(row["field_skew_ns"], 800000)
        self.assertEqual(row["latitude_header_timestamp_ns"], row["timestamp_ns"])

    def test_only_exact_integrated_fixed_status_passes(self):
        for state in (0, 1, 3):
            self.assertIn("not_integrated_navigation", one_row(navigation_cycle(overrides={"status": {"system_state": state}}))["reason"])
        for sat in (0, 1, 2, 3, 5, 6, 7, 8, 9):
            self.assertIn("not_rtk_fixed_and_heading_valid", one_row(navigation_cycle(overrides={"status": {"satellite_status": sat}}))["reason"])

    def test_sigma_differential_age_and_missing_each_field(self):
        for field, changes, reason in [
            ("position_sigma", {"east_m": float("nan")}, "invalid_sigma"),
            ("position_sigma", {"up_m": .18}, "position_sigma"),
            ("attitude_sigma", {"heading_deg": .5}, "attitude_sigma"),
            ("velocity_sigma", {"up_mps": -.1}, "invalid_sigma"),
            ("status", {"differential_age_s": 1.4}, "differential_age"),
        ]:
            self.assertIn(reason, one_row(navigation_cycle(overrides={field: changes}))["reason"])
        for field in ex.FIELDS:
            if field != "latitude":
                records = [r for r in navigation_cycle() if r.topic != "/cgi430/" + ex.FIELDS[field][0]]
                self.assertIn("missing_field:" + field, one_row(records)["reason"])

    def test_age_skew_and_unfixed_37_second_clock_offset(self):
        records = navigation_cycle()
        records[-1].bag_ns += 37 * ex.NS
        self.assertIn("stale_field", one_row(records)["reason"])
        c = config()
        c["clock"]["reference_bag_offset_ns"] = -37 * ex.NS
        for record in records[:-1]:
            record.bag_ns += 37 * ex.NS
        self.assertEqual(one_row(records, c)["reference_valid"], 1)
        c = config()
        c["quality"]["max_skew_sec"] = .0005
        self.assertIn("field_time_skew", one_row(navigation_cycle(), c)["reason"])

    def test_primary_must_advance_and_recovery_reset(self):
        c = config()
        c["quality"]["recovery_samples"] = 2
        assembly = ex.ReferenceAssembler(c)
        rows = []
        for index in range(4):
            override = {"status": {"satellite_status": 5}} if index == 2 else None
            for rec in navigation_cycle(1_700_000_000_000_000_000 + index * 10_000_000, override):
                row = assembly.ingest(rec)
                if row is not None:
                    rows.append(row)
        self.assertEqual([r["reference_valid"] for r in rows], [0, 1, 0, 0])
        self.assertEqual(rows[0]["reason"], "recovering")
        repeated = navigation_cycle(1_700_000_000_040_000_000)[-1]
        self.assertIn("navigation_field_not_advanced", assembly.ingest(repeated)["reason"])

    def test_invalid_duplicate_and_clock_rollback_retained(self):
        records = navigation_cycle()
        assembly = ex.ReferenceAssembler(config())
        for record in records:
            assembly.ingest(record)
        duplicate = assembly.ingest(records[-1])
        self.assertEqual(duplicate["reference_valid"], 0)
        self.assertIn("duplicate_reference_timestamp", duplicate["reason"])
        rollback = navigation_cycle(1_600_000_000_000_000_000)[-1]
        self.assertIn("reference_clock_rollback", assembly.ingest(rollback)["reason"])

    def test_transient_status_loss_before_next_anchor_breaks_reference(self):
        first = navigation_cycle()
        assembly = ex.ReferenceAssembler(config())
        for record in first:
            assembly.ingest(record)
        bad = next(r for r in navigation_cycle(1_700_000_000_005_000_000, {"status": {"satellite_status": 5}}) if r.topic.endswith("/ins/status"))
        assembly.ingest(bad)
        for record in navigation_cycle(1_700_000_000_010_000_000):
            row = assembly.ingest(record)
        self.assertEqual(row["reference_valid"], 0)
        self.assertIn("intervening_rejection:nonfixed_status:status", row["reason"])

    def test_height_conversion_is_explicit(self):
        c = config()
        c["height"].update(datum="orthometric", geoid_separation_m=20.)
        self.assertAlmostEqual(one_row(navigation_cycle(), c)["z"], 30., places=8)

    def test_fail_closed_export_but_raw_capture_can_calibrate(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.yaml"
            raw = yaml.safe_load((Path(__file__).resolve().parents[1] / "config/accuracy_capture.yaml").read_text())
            raw.update(source_recording_id="first-calibration", epoch="first")
            raw["capture"].update(confirmed=True, evidence="Verified local topics", expected_can_bitrate=500000)
            path.write_text(yaml.safe_dump(raw), encoding="utf-8")
            ex.load_config(path, "capture", "reference")
            with self.assertRaises(ValueError):
                ex.load_config(path)
            c = config()
            path.write_text(yaml.safe_dump(c), encoding="utf-8")
            ex.load_config(path, require_estimate=True)
            c["clock"]["confirmed"] = False
            path.write_text(yaml.safe_dump(c), encoding="utf-8")
            with self.assertRaises(ValueError):
                ex.load_config(path)
            c["clock"]["assumption_acknowledged"] = True
            path.write_text(yaml.safe_dump(c), encoding="utf-8")
            ex.load_config(path)
            c["capture"]["input_topics"].append("/tf")
            path.write_text(yaml.safe_dump(c), encoding="utf-8")
            with self.assertRaises(ValueError):
                ex.load_config(path, "capture", "algorithm")

    def test_raw_timeline_rejections_and_recording_tail_preserved(self):
        c = config()
        records = navigation_cycle() + navigation_cycle(1_700_000_000_010_000_000, {"status": {"satellite_status": 5}})
        records.append(ex.Record("/unknown/heartbeat", "missing_interfaces/msg/Heartbeat", 1_700_000_020_000_000_000, None, 100))
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            summary = ex.export_stream(records, c, out, "reference")
            self.assertEqual(summary["rows"], 2)
            self.assertEqual(summary["valid_rows"], 1)
            self.assertEqual(summary["end_ns"], records[-1].bag_ns)
            with (out / "reference.csv").open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([r["reference_valid"] for r in rows], ["1", "0"])
            with (out / "reference.raw_frames.csv").open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 18)
            with (out / "reference.timeline.csv").open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 19)
            self.assertTrue((out / "reference.tum").read_text().startswith("1700000000.000800000 "))
            meta = ex.interval_metadata(c, "reference", "synthetic", None, summary["start_ns"], summary["end_ns"])
            self.assertEqual(meta["recorded_intervals"][0]["end_ns"], records[-1].bag_ns)
            self.assertIn("cgi_field_measurement_epoch_unverified", meta["conditions"])

    def test_playback_uses_source_interval_not_today(self):
        c = config()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "capture.json"
            path.write_text(json.dumps({"source_recording_id": c["source_recording_id"], "epoch": c["epoch"],
                                        "capture_start_ns": 10**18, "capture_end_ns": 10**18 + 1000,
                                        "exit_code": 0, "role": "playback",
                                        "playback_source_interval": {"start_ns": 100, "end_ns": 1000}}))
            meta = ex.interval_metadata(c, "estimate", "synthetic", path, 10**18, 10**18 + 1000)
            self.assertEqual(meta["recorded_intervals"], [{"start_ns": 100, "end_ns": 1000, "epoch": c["epoch"]}])

    def test_recording_boundary_never_expands_to_bad_clock(self):
        c = config()
        start, end = 1_700_000_000_000_000_000, 1_700_000_020_000_000_000
        for bounds in ((0, end), (start - 37 * ex.NS, end), (start, end + 37 * ex.NS)):
            with self.assertRaisesRegex(ValueError, "interval|boundary"):
                ex.interval_metadata(c, "reference", "synthetic", None, start, end, bounds)
        meta = ex.interval_metadata(c, "reference", "synthetic", None, start, end, (start - 1_000_000, end))
        self.assertEqual(meta["recorded_intervals"][0]["end_ns"] - meta["recorded_intervals"][0]["start_ns"], 20 * ex.NS + 1_000_000)

    def test_invalid_zero_anchor_preserves_audit_but_no_metadata(self):
        c = config()
        records = navigation_cycle()
        bad = copy.deepcopy(records[-1])
        bad.message["frame"]["header"]["stamp"] = {"sec": 0, "nanosec": 0}
        records.insert(0, bad)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "result"
            args = SimpleNamespace(output=str(out), reference_bag="ref", estimate_bag=None, reference_capture=None, estimate_capture=None)
            with patch.object(ex, "bag_records", lambda *_: iter(records)):
                with self.assertRaisesRegex(ValueError, "timestamp <=0"):
                    ex.export_bags(args, c)
            self.assertTrue((out / "reference.csv").exists())
            self.assertTrue((out / "reference.raw_frames.csv").exists())
            self.assertTrue((out / "reference.timeline.csv").exists())
            self.assertTrue((out / "export.failure.json").exists())
            self.assertFalse((out / "reference.metadata.json").exists())

    def test_missing_pose_publisher_retains_empty_tum_and_full_denominator(self):
        c = config()
        c["capture"]["required_output_topics"] = []
        status = [ex.Record("/localization/loc_status", "std_msgs/msg/String", 1_700_000_000_000_000_000 + i * ex.NS,
                            {"data": "initialization_failed"}, 32) for i in range(3)]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "result"
            args = SimpleNamespace(output=str(out), reference_bag="ref", estimate_bag="est", reference_capture=None, estimate_capture=None)
            def records(bag, required):
                if bag == "est":
                    self.assertEqual(required, [])
                    return iter(status)
                return iter(navigation_cycle() + navigation_cycle(1_700_000_000_010_000_000))
            with patch.object(ex, "bag_records", records), patch("builtins.print"):
                rc = ex.export_bags(args, c)
            self.assertEqual(rc, 0)
            self.assertEqual((out / "estimated.tum").read_text(), "")
            metadata = json.loads((out / "estimated.metadata.json").read_text())
            self.assertEqual(metadata["recorded_intervals"][0]["end_ns"], status[-1].bag_ns)
            summary = json.loads((out / "export.summary.json").read_text())["estimate"]
            self.assertEqual(summary["rows"], 0)
            self.assertFalse(summary["satisfied"])
            self.assertIn("/slamPoseRaw_topic", summary["missing_output_topics"])
            self.assertIn("ESTIMATE_STREAM_ABSENT", metadata["conditions"])
            self.assertEqual(len((out / "estimated.status.jsonl").read_text().splitlines()), 3)


def ros_smoke(directory, storage):
    import rosbag2_py
    from rclpy.serialization import serialize_message
    from rosidl_runtime_py.set_message import set_message_fields
    from rosidl_runtime_py.utilities import get_message
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=False)
    c = config()
    c["capture"]["storage_id"] = storage
    config_path = out / "config.yaml"
    config_path.write_text(yaml.safe_dump(c), encoding="utf-8")
    writer = rosbag2_py.SequentialWriter()
    writer.open(rosbag2_py.StorageOptions(uri=str(out / "reference_bag"), storage_id=storage), rosbag2_py.ConverterOptions("", ""))
    for key, (suffix, typename, _) in ex.FIELDS.items():
        writer.create_topic(rosbag2_py.TopicMetadata(name="/cgi430/" + suffix, type="cgi430_interfaces/msg/" + typename, serialization_format="cdr"))
    records = navigation_cycle() + navigation_cycle(1_700_000_000_010_000_000, {"status": {"satellite_status": 5}}) + navigation_cycle(1_700_000_000_020_000_000)
    for record in records:
        message = get_message(record.type_name)()
        set_message_fields(message, record.message)
        writer.write(record.topic, serialize_message(message), record.bag_ns)
    del writer
    writer = rosbag2_py.SequentialWriter()
    writer.open(rosbag2_py.StorageOptions(uri=str(out / "estimate_bag"), storage_id=storage), rosbag2_py.ConverterOptions("", ""))
    typename = "geometry_msgs/msg/PoseStamped"
    writer.create_topic(rosbag2_py.TopicMetadata(name=c["estimate"]["topic"], type=typename, serialization_format="cdr"))
    for n in range(3):
        timestamp = 1_700_000_000_000_800_000 + n * 10_000_000
        message = get_message(typename)()
        message.header.stamp.sec, message.header.stamp.nanosec = divmod(timestamp, ex.NS)
        message.header.frame_id = "map"
        message.pose.position.z = 10.
        message.pose.orientation.w = 1.
        writer.write(c["estimate"]["topic"], serialize_message(message), timestamp + 1_000_000)
    del writer
    rc = ex.main(["export", "--config", str(config_path), "--reference-bag", str(out / "reference_bag"),
                  "--estimate-bag", str(out / "estimate_bag"), "--output", str(out / "export")])
    if rc != 0:
        raise RuntimeError(f"real {storage} export failed: {rc}")
    summary = json.loads((out / "export/export.summary.json").read_text())
    assert summary["reference"]["rows"] == 3 and summary["reference"]["valid_rows"] == 2
    assert summary["estimate"]["rows"] == 3
    recording = json.loads((out / "export/estimated.recording.metadata.json").read_text())
    assert not {"world_frame_id", "body_frame_id", "trajectory_kind"}.intersection(recording)
    assert recording["evidence"]["applied_time_corrections"]["header_offset_ns"] == 0
    print("ROS storage smoke passed:", storage, out)


def ros_live_smoke(directory, domain):
    """Actual wrapper lifecycle in one explicitly isolated localhost ROS domain."""
    if not 1 <= domain <= 232:
        raise ValueError("choose an unused test domain 1..232")
    os.environ.update(ROS_DOMAIN_ID=str(domain), ROS_LOCALHOST_ONLY="1")
    import rclpy
    from rclpy.node import Node
    from rclpy.qos import QoSProfile, ReliabilityPolicy
    from rclpy.serialization import serialize_message
    from geometry_msgs.msg import PoseStamped
    import rosbag2_py

    out = Path(directory).resolve()
    out.mkdir(parents=True, exist_ok=False)
    c = config()
    c["capture"].update(input_topics=["/synthetic/input_pose"], output_topics=["/slamPoseRaw_topic"],
                        required_output_topics=["/slamPoseRaw_topic"], evidence_topics=[], storage_id="sqlite3")
    path = out / "config.yaml"
    path.write_text(yaml.safe_dump(c), encoding="utf-8")
    script = Path(__file__).resolve().parent / "accuracy_capture.sh"
    source_bag = out / "historical_input_bag"
    writer = rosbag2_py.SequentialWriter()
    writer.open(rosbag2_py.StorageOptions(uri=str(source_bag), storage_id="sqlite3"), rosbag2_py.ConverterOptions("", ""))
    writer.create_topic(rosbag2_py.TopicMetadata(name="/synthetic/input_pose", type="geometry_msgs/msg/PoseStamped", serialization_format="cdr"))
    source_start = 1_700_000_000_000_000_000
    for index in range(61):
        stamp = source_start + index * 50_000_000
        message = PoseStamped()
        message.header.stamp.sec, message.header.stamp.nanosec = divmod(stamp, ex.NS)
        message.header.frame_id = "map"
        message.pose.orientation.w = 1.
        writer.write("/synthetic/input_pose", serialize_message(message), stamp)
    del writer
    rclpy.init()
    node = Node("accuracy_capture_synthetic_smoke", enable_rosout=False, start_parameter_services=False)
    qos = QoSProfile(depth=100, reliability=ReliabilityPolicy.BEST_EFFORT)
    publisher = node.create_publisher(PoseStamped, "/slamPoseRaw_topic", qos)
    subscription = node.create_subscription(PoseStamped, "/synthetic/input_pose", publisher.publish, qos)
    child = None
    try:
        for mode in ("record", "playback"):
            target = out / mode
            command = ["bash", str(script), mode, "--config", str(path), "--output", str(target)]
            command += ["--role", "outputs"] if mode == "record" else ["--bag", str(source_bag), "--domain-id", str(domain), "--post-wait-seconds", "1"]
            with (out / (mode + ".wrapper.log")).open("w") as log:
                child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=os.environ.copy())
                deadline, begun, stopped, published_after_ready = time.monotonic() + 70, None, False, 0
                while child.poll() is None:
                    if time.monotonic() > deadline:
                        raise RuntimeError(f"{mode} wrapper timed out; see {log.name}")
                    if mode == "record":
                        stamp = time.time_ns()
                        message = PoseStamped()
                        message.header.stamp.sec, message.header.stamp.nanosec = divmod(stamp, ex.NS)
                        message.header.frame_id = "map"
                        message.pose.orientation.w = 1.
                        publisher.publish(message)
                        if begun is not None:
                            published_after_ready += 1
                        recorder_log = target / "logs/recorder.stderr.log"
                        # capture_start includes CLI/discovery latency. A fixed
                        # delay from that mark can expire before ROS subscribes.
                        # Require this recorder's observed readiness and a DDS
                        # match, then exercise a real three-second data window.
                        if begun is None and publisher.get_subscription_count() > 0 and recorder_log.exists():
                            if "Recording..." in recorder_log.read_text(encoding="utf-8"):
                                begun = time.monotonic()
                        if begun is not None and time.monotonic() - begun >= 3 and not stopped:
                            ex.write_json(target / "fixture.readiness.json", {
                                "matched_subscriptions": publisher.get_subscription_count(),
                                "published_after_ready": published_after_ready,
                                "seconds_after_ready": time.monotonic() - begun})
                            child.send_signal(signal.SIGINT)
                            stopped = True
                    rclpy.spin_once(node, timeout_sec=.02)
                if child.returncode != 0:
                    raise RuntimeError(f"{mode} wrapper exited {child.returncode}; see {log.name}")
                child = None
            metadata = json.loads((target / "capture.metadata.json").read_text())
            assert metadata["exit_code"] == 0 and (target / "bag/metadata.yaml").exists()
            messages = list(ex.bag_records(target / "bag", [c["estimate"]["topic"]]))
            assert len(messages) > 0, mode + " wrote no synthetic output poses"
            if mode == "playback":
                assert all(ex.stamp_ns(m.message) < 1_701_000_000_000_000_000 for m in messages)
                meta = ex.interval_metadata(c, "estimate", target / "bag", target / "capture.metadata.json", min(m.bag_ns for m in messages), max(m.bag_ns for m in messages))
                assert meta["recorded_intervals"][0]["start_ns"] == source_start
                assert meta["recorded_intervals"][0]["end_ns"] == source_start + 3 * ex.NS
                ex.write_json(target / "checked.estimated.metadata.json", meta)
    finally:
        if child is not None and child.poll() is None:
            child.send_signal(signal.SIGINT)
            try:
                child.wait(timeout=15)
            except subprocess.TimeoutExpired:
                child.terminate()
                child.wait(timeout=5)
        node.destroy_subscription(subscription)
        node.destroy_node()
        rclpy.shutdown()
    # Failure injection: a recorder can create metadata and still fail flushing.
    # PATH is changed only for this child; this fake CLI never creates a ROS node.
    fake_bin = out / "fake_recorder_bin"
    fake_bin.mkdir()
    fake_ros = fake_bin / "ros2"
    fake_ros.write_text("""#!/usr/bin/env python3
import pathlib, signal, sys, time
if sys.argv[1:3] == ['topic', 'list']:
    raise SystemExit(0)
if sys.argv[1:3] == ['bag', 'record']:
    target = pathlib.Path(sys.argv[sys.argv.index('--output') + 1])
    target.mkdir()
    (target / 'metadata.yaml').write_text('synthetic_final_flush_failure: true\\n')
    def fail_flush(signum, frame):
        raise SystemExit(23)
    signal.signal(signal.SIGINT, fail_flush)
    while True:
        time.sleep(.05)
raise SystemExit(2)
""", encoding="utf-8")
    fake_ros.chmod(0o755)
    fake_env = dict(os.environ, PATH=str(fake_bin) + os.pathsep + os.environ["PATH"])
    target = out / "flush_failure"
    command = ["bash", str(script), "record", "--config", str(path), "--role", "outputs", "--output", str(target), "--duration-seconds", "1"]
    proc = subprocess.run(command, env=fake_env, text=True, capture_output=True, timeout=50, check=False)
    (out / "flush_failure.wrapper.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
    metadata = json.loads((target / "capture.metadata.json").read_text())
    assert proc.returncode == 3 and metadata["exit_code"] == 3, proc.returncode
    assert metadata["child_exit_codes"]["recorder"] == 23
    print("Actual record/SIGINT/playback wrapper smoke passed in localhost domain", domain, out)


def core_smoke(directory):
    """Decoded mock records -> actual CSV/TUM/metadata -> core fit/evaluate CLIs."""
    out = Path(directory).resolve()
    out.mkdir(parents=True, exist_ok=False)
    core = Path(__file__).resolve().parent / "accuracy_eval/accuracy_eval.py"
    if not core.exists():
        raise RuntimeError("accuracy_eval core is not available")
    configs = {}
    for run, time_offset in (("B", 0), ("C", 100 * ex.NS)):
        folder = out / run
        folder.mkdir()
        c = config()
        c["source_recording_id"] = "SYNTHETIC_CAPTURE_CORE_" + run
        refs, estimates = [], []
        for index in range(60):
            t = index * .09
            lat, lon, altitude = 0.0001 * math.sin(t), 0.0001 * math.cos(t), 10 + math.sin(t * 2)
            heading, pitch, roll = (90 - index * 4) % 360, math.sin(t), math.cos(t)
            base = 1_700_000_000_000_000_000 + time_offset + index * 100_000_000
            cycle = navigation_cycle(base, {"latitude": {"latitude_deg": lat}, "longitude": {"longitude_deg": lon},
                                           "altitude": {"altitude_m": altitude}, "attitude": {"heading_deg": heading, "pitch_deg": pitch, "roll_deg": roll}})
            refs.extend(cycle)
            stamp = ex.stamp_ns(cycle[-1].message, "frame.header.stamp")
            p, q = ex.convert_pose([lat, lon, altitude], [heading, pitch, roll], c["georeference"]["origin_llh"])
            if run == "C":
                p += [.2, -.1, .05]
            pose = {"header": {"stamp": {"sec": stamp // ex.NS, "nanosec": stamp % ex.NS}, "frame_id": "map"},
                    "pose": {"position": dict(zip("xyz", map(float, p))), "orientation": dict(zip("xyzw", map(float, q)))}}
            estimates.append(ex.Record(c["estimate"]["topic"], "geometry_msgs/msg/PoseStamped", stamp + 1_000_000, pose, 80))
        for role, records, stem in (("reference", refs, "reference"), ("estimate", estimates, "estimated")):
            summary = ex.export_stream(records, c, folder, role)
            meta = ex.interval_metadata(c, role, "SYNTHETIC_DECODED_RECORDS", None, summary["start_ns"], summary["end_ns"],
                                        (summary["all_pose_start_ns"], summary["all_pose_end_ns"]))
            ex.write_json(folder / (stem + ".metadata.json"), meta)
        reference_meta = json.loads((folder / "reference.metadata.json").read_text())
        item = {"name": "localization", "tum": run + "/estimated.tum", "metadata": run + "/estimated.metadata.json"}
        if run == "C":
            item["alignment"] = "frozen.json"
        settings = {"schema_version": 1, "role": "calibration" if run == "B" else "evaluation",
                    "reference": {"csv": run + "/reference.csv", "metadata": run + "/reference.metadata.json"},
                    "estimates": [item], "analysis_intervals": reference_meta["recorded_intervals"],
                    "time": {"offset_status": "known", "estimate_to_reference_offset_ns": 0, "reference_epoch": c["epoch"],
                             "evidence": "Synthetic timestamps only; no field precision claim", "allow_conditional": False},
                    "pairing": {"max_reference_gap_sec": .12, "availability_hold_sec": .12},
                    "conditions": ["SYNTHETIC_EXPORT_CORE_SMOKE"], "thresholds": None}
        if run == "B":
            settings["fit"] = {"method": "position_attitude", "min_pairs": 30, "min_secondary_span_m": .5,
                               "attitude_weight_m": 2., "max_attitude_dispersion_deg": 1.}
        target = out / (run + ".yaml")
        target.write_text(yaml.safe_dump(settings, sort_keys=False), encoding="utf-8")
        configs[run] = target
    commands = [[sys.executable, str(core), "validate", "--config", str(configs["B"])],
                [sys.executable, str(core), "fit", "--config", str(configs["B"]), "--estimate", "localization", "--output", str(out / "frozen.json")],
                [sys.executable, str(core), "validate", "--config", str(configs["C"])],
                [sys.executable, str(core), "evaluate", "--config", str(configs["C"]), "--output", str(out / "report")]]
    for index, command in enumerate(commands):
        proc = subprocess.run(command, text=True, capture_output=True, timeout=30, check=False)
        (out / f"core_{index}.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
        if proc.returncode:
            raise RuntimeError(f"core CLI failed: {command}\n{proc.stdout}\n{proc.stderr}")
    report = json.loads((out / "report/report.json").read_text(encoding="utf-8"))
    actual = report["estimates"][0]["metrics"]["position_3d_m"]["rmse"]
    assert abs(actual - math.sqrt(.2**2 + .1**2 + .05**2)) < 1e-8, actual
    assert "reference_timestamp_is_socketcan_receive_proxy" in json.dumps(report)
    print("Export -> actual core validate/fit/evaluate CLI smoke passed:", out)


if __name__ == "__main__":
    if "--ros-live-smoke" in sys.argv:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--ros-live-smoke", required=True)
        parser.add_argument("--domain-id", type=int, default=221)
        args = parser.parse_args()
        ros_live_smoke(args.ros_live_smoke, args.domain_id)
    elif "--core-smoke" in sys.argv:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--core-smoke", required=True)
        args = parser.parse_args()
        core_smoke(args.core_smoke)
    elif "--ros-smoke" in sys.argv:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--ros-smoke", required=True)
        parser.add_argument("--storage", choices=["mcap", "sqlite3"], default="mcap")
        args = parser.parse_args()
        ros_smoke(args.ros_smoke, args.storage)
    else:
        unittest.main()
