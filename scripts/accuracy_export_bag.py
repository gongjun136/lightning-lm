#!/usr/bin/env python3
"""Read-only ROS 2 accuracy capture planning and audited CGI reference export.

ROS imports are lazy: --help, configuration checks and synthetic tests need only
Python, NumPy and PyYAML. Humble rosbag2_py reads sqlite3 and MCAP through installed
storage plugins. No ROS publisher, CAN transmitter or production node is created.

Coordinate/quality semantics mirror core/navigation/ins_navigation.{h,cc}:
fixed WGS84 ENU, R_fixed_local @ Rz(90-heading) @ Ry(-pitch) @ Rx(roll), system=2
and satellite=4. Position remains at the declared CGI output point: a confirmed
fixed point or an already compensated rear axle. No lever arm is applied here.
Quaternion axes always describe rear-body FLU. Offline arrival checks use
recorded bag time, not a fabricated steady clock or measurement epoch.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass

import numpy as np
import yaml

FIELDS = {
    "latitude": ("position/latitude", "Latitude", ("latitude_deg",)),
    "longitude": ("position/longitude", "Longitude", ("longitude_deg",)),
    "altitude": ("position/altitude", "Altitude", ("altitude_m",)),
    "attitude": ("attitude", "Attitude", ("heading_deg", "pitch_deg", "roll_deg")),
    "velocity": ("velocity", "EarthVelocity", ("east_mps", "north_mps", "up_mps")),
    "position_sigma": ("position/sigma", "PositionSigma", ("east_m", "north_m", "up_m")),
    "attitude_sigma": ("attitude/sigma", "AttitudeSigma", ("heading_deg", "pitch_deg", "roll_deg")),
    "status": ("ins/status", "InsStatus", ("system_state", "satellite_status", "differential_age_s")),
    "velocity_sigma": ("velocity/sigma", "VelocitySigma", ("east_mps", "north_mps", "up_mps")),
}
PRIMARY_FIELDS = tuple(FIELDS)[:5]
RAW_SUFFIXES = ["time", "imu/angular_rate_raw", "imu/acceleration_raw", "ins/status",
                "position/altitude", "position/sigma", "velocity", "velocity/sigma",
                "vehicle/acceleration", "attitude", "attitude/sigma",
                "vehicle/angular_rate", "position/longitude", "position/latitude"]
CORE_COLUMNS = ["timestamp_ns", "x", "y", "z", "qx", "qy", "qz", "qw",
                "reference_valid", "reason", "epoch", "time_basis"]
NS = 1_000_000_000


def fail(message):
    raise ValueError(message)


def finite(value, label, minimum=None, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        fail(f"{label} must be an explicitly configured finite number")
    if minimum is not None and (value < minimum or (positive and value == minimum)):
        fail(f"{label} must be {'>' if positive else '>='} {minimum}")
    return value


def nonempty(value, label):
    if not isinstance(value, str) or not value.strip() or value.strip().lower() in {"todo", "unknown", "null"}:
        fail(f"{label} must name verified experiment evidence")
    return value


def validate_capture(c, role=None):
    """Raw acquisition is allowed before clock/geodetic calibration exists."""
    cap = c.get("capture", {})
    if cap.get("confirmed") is not True:
        fail("capture.confirmed must be true after topic ownership/settings verification")
    nonempty(cap.get("evidence"), "capture.evidence")
    if cap.get("storage_id") not in {"mcap", "sqlite3"}:
        fail("capture.storage_id must be mcap or sqlite3")
    if not re.fullmatch(r"/[A-Za-z0-9_/]+", str(cap.get("cgi_prefix", ""))) or cap["cgi_prefix"].endswith("/"):
        fail("capture.cgi_prefix must be an absolute namespace without trailing slash")
    if not re.fullmatch(r"[A-Za-z0-9_]+", str(cap.get("can_interface", ""))):
        fail("invalid capture.can_interface")
    bitrate = cap.get("expected_can_bitrate")
    if role in {"reference", "combined"} and (isinstance(bitrate, bool) or not isinstance(bitrate, int) or bitrate <= 0):
        fail("capture.expected_can_bitrate must be verified, not guessed")
    for key in ("input_topics", "output_topics", "evidence_topics", "required_output_topics"):
        if not isinstance(cap.get(key), list) or (key in {"input_topics", "output_topics"} and not cap[key]):
            fail("capture." + key + " must be a topic list")
        for topic in cap[key]:
            if not isinstance(topic, str) or not re.fullmatch(r"/[A-Za-z0-9_/]+", topic):
                fail(f"invalid ROS topic {topic!r}")
        if len(cap[key]) != len(set(cap[key])):
            fail("duplicate topics in capture." + key)
    prohibited = set(cap["output_topics"] + cap["evidence_topics"] + ["/clock", "/tf", "/tf_static", "/PosRes", "/slamPoseRaw_topic"])
    if prohibited.intersection(cap["input_topics"]) or any(t.startswith(cap["cgi_prefix"] + "/") or t.startswith("/localization/") or t.startswith("/diagnostics/") for t in cap["input_topics"]):
        fail("input replay allowlist contains an output, CGI reference, TF, clock or diagnostics topic")
    if not set(cap["required_output_topics"]) <= set(cap["output_topics"]):
        fail("required output topics must be listed in capture.output_topics")


def reference_body_frame(ref):
    """Identify the output point without translating it or rotating body axes."""
    if ref.get("attitude_axes") != "FLU":
        fail("reference_point.attitude_axes must be FLU for the unified rear-body attitude")
    if ref.get("position") == "rear_axle_center":
        if ref.get("device_lever_arm_already_applied") is not True:
            fail("rear_axle_center requires device_lever_arm_already_applied:true; no lever arm is applied")
        return "rear_axle_flu"
    if ref.get("position") == "fixed_point":
        nonempty(ref.get("point_id"), "reference_point.point_id")
        if ref.get("device_lever_arm_already_applied") is not False:
            fail("fixed_point requires device_lever_arm_already_applied:false; no lever arm is applied")
        if ref.get("rigid_to_rear_axle_confirmed") is not True:
            fail("fixed_point requires rigid_to_rear_axle_confirmed:true; a static lever arm cannot cross articulation")
        return "cgi_output_flu"
    fail("reference_point.position must be fixed_point or rear_axle_center")


def load_config(path, purpose="export", role=None, require_estimate=False):
    with Path(path).open(encoding="utf-8") as stream:
        c = yaml.safe_load(stream)
    if not isinstance(c, dict) or c.get("schema_version") != 1:
        fail("configuration schema_version must be 1")
    for key in ("source_recording_id", "epoch"):
        nonempty(c.get(key), key)
    validate_capture(c, role)
    if purpose == "capture":
        return c
    for key in ("time_basis", "world_frame_id"):
        nonempty(c.get(key), key)
    for section in ("georeference", "height", "reference_point", "quality") + (("estimate",) if require_estimate else ()):
        block = c.get(section, {})
        if block.get("confirmed") is not True:
            fail(f"{section}.confirmed must be true after field verification")
        nonempty(block.get("evidence"), section + ".evidence")
    geo = c["georeference"]
    if geo.get("datum") != "WGS84" or not isinstance(geo.get("origin_llh"), list) or len(geo["origin_llh"]) != 3:
        fail("georeference requires WGS84 origin_llh [latitude, longitude, ellipsoid height]")
    lat, lon, height = [finite(v, "origin_llh") for v in geo["origin_llh"]]
    if abs(lat) >= 90 or abs(lon) > 180:
        fail("invalid origin; poles are unsupported, matching ins::GeoReference")
    h = c["height"]
    if h.get("datum") not in {"ellipsoidal", "orthometric"}:
        fail("height.datum must be ellipsoidal or orthometric")
    finite(h.get("geoid_separation_m"), "height.geoid_separation_m")
    if h["datum"] == "ellipsoidal" and h["geoid_separation_m"] != 0:
        fail("ellipsoidal input requires geoid_separation_m=0")
    reference_body_frame(c["reference_point"])
    clock = c["clock"]
    for section in ("clock", "assembly"):
        block = c[section]
        if block.get("confirmed") is not True and block.get("assumption_acknowledged") is not True:
            fail(f"{section} must be confirmed or explicitly assumption_acknowledged; reports will retain this condition")
        nonempty(block.get("evidence"), section + ".evidence")
    if clock.get("reference_header_semantics") != "socketcan_kernel_receive":
        fail("typed CGI header is socketcan_kernel_receive, not measurement epoch")
    offset_keys = ("reference_header_offset_ns", "reference_bag_offset_ns", "uncertainty_ns")
    if require_estimate:
        offset_keys += ("estimate_header_offset_ns", "estimate_bag_offset_ns")
    for key in offset_keys:
        if isinstance(clock.get(key), bool) or not isinstance(clock.get(key), int):
            fail(f"clock.{key} must be an explicit integer, including an audited zero")
    if clock["uncertainty_ns"] < 0:
        fail("clock.uncertainty_ns cannot be negative")
    nonempty(clock.get("latency_evidence"), "clock.latency_evidence")
    if c["assembly"].get("method") != "causal_record_order_at_latitude":
        fail("unsupported assembly method")
    for key in ("max_age_sec", "max_skew_sec", "max_differential_age_sec", "max_horizontal_sigma_m", "max_vertical_sigma_m", "max_attitude_sigma_deg", "max_velocity_sigma_mps", "max_site_radius_m"):
        finite(c["quality"].get(key), "quality." + key, 0, True)
    finite(c["quality"].get("max_future_sec"), "quality.max_future_sec", 0)
    n = c["quality"].get("recovery_samples")
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        fail("quality.recovery_samples must be an explicitly chosen integer >=1")
    if require_estimate:
        est = c["estimate"]
        if est.get("body_frame_id") != "rear_axle_flu":
            fail("estimate.body_frame_id must be rear_axle_flu; normalize internal IMU TUM separately")
        nonempty(est.get("world_frame_id"), "estimate.world_frame_id")
        if not est["world_frame_id"].startswith(("map:", "run:")):
            fail("estimate.world_frame_id must identify map:<asset hash> or run:<unique run id>")
        if est.get("trajectory_kind") not in {"frontend", "localization"}:
            fail("published pose export requires estimate.trajectory_kind frontend or localization")
        if est.get("topic") not in c["capture"]["output_topics"]:
            fail("estimate topic must be listed in capture.output_topics")
    return c


def ecef(llh):
    """WGS84 ellipsoid conversion; uses the same ellipsoid as GeographicLib."""
    lat, lon, h = llh
    lat, lon = math.radians(lat), math.radians(lon)
    a, f = 6378137.0, 1.0 / 298.257223563
    e2 = f * (2 - f)
    n = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    return np.array([(n + h) * math.cos(lat) * math.cos(lon),
                     (n + h) * math.cos(lat) * math.sin(lon),
                     (n * (1 - e2) + h) * math.sin(lat)])


def enu_from_ecef(lat_deg, lon_deg):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sl, cl, so, co = math.sin(lat), math.cos(lat), math.sin(lon), math.cos(lon)
    return np.array([[-so, co, 0], [-sl * co, -sl * so, cl], [cl * co, cl * so, sl]])


def vehicle_attitude(heading, pitch, roll):
    yaw, pitch, roll = map(math.radians, (90 - heading, -pitch, roll))
    cy, sy, cp, sp, cr, sr = math.cos(yaw), math.sin(yaw), math.cos(pitch), math.sin(pitch), math.cos(roll), math.sin(roll)
    rz = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]])
    ry = np.array([[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]])
    rx = np.array([[1, 0, 0], [0, cr, -sr], [0, sr, cr]])
    return rz @ ry @ rx


def quaternion_xyzw(matrix):
    """Stable matrix-to-quaternion conversion, including rotations near pi."""
    m = matrix
    trace = float(np.trace(m))
    if trace > 0:
        s = math.sqrt(trace + 1) * 2
        q = np.array([(m[2, 1] - m[1, 2]) / s, (m[0, 2] - m[2, 0]) / s,
                      (m[1, 0] - m[0, 1]) / s, s / 4])
    else:
        i = int(np.argmax(np.diag(m)))
        j, k = (i + 1) % 3, (i + 2) % 3
        s = math.sqrt(max(0, 1 + m[i, i] - m[j, j] - m[k, k])) * 2
        q = np.zeros(4)
        q[i], q[j], q[k], q[3] = s / 4, (m[j, i] + m[i, j]) / s, (m[k, i] + m[i, k]) / s, (m[k, j] - m[j, k]) / s
    q /= np.linalg.norm(q)
    return q if q[3] >= 0 else -q


def convert_pose(llh, attitude, origin):
    fixed = enu_from_ecef(*origin[:2])
    p = fixed @ (ecef(llh) - ecef(origin))
    fixed_from_local = fixed @ enu_from_ecef(*llh[:2]).T
    q = quaternion_xyzw(fixed_from_local @ vehicle_attitude(*attitude))
    return p, q


def get(obj, path, default=None):
    for key in path.split("."):
        if isinstance(obj, dict):
            obj = obj.get(key, default)
        else:
            obj = getattr(obj, key, default)
        if obj is None:
            return default
    return obj


def stamp_ns(message, prefix="header.stamp"):
    sec, ns = get(message, prefix + ".sec"), get(message, prefix + ".nanosec")
    if sec is None or ns is None:
        return None
    return int(sec) * NS + int(ns)


@dataclass
class Record:
    topic: str
    type_name: str
    bag_ns: int
    message: object
    size: int = 0


class ReferenceAssembler:
    """Conservative causal assembly; every latitude arrival produces an audit row."""
    def __init__(self, config):
        self.c = config
        self.latest = {}
        self.last_stamp = {}
        self.consumed = {}
        self.good_count = 0
        self.previous_anchor = None
        self.intervening_rejections = []
        self.topic_field = {config["capture"]["cgi_prefix"] + "/" + v[0]: k for k, v in FIELDS.items()}

    def ingest(self, rec):
        key = self.topic_field.get(rec.topic)
        if key is None:
            return None
        expected = "cgi430_interfaces/msg/" + FIELDS[key][1]
        if rec.type_name != expected:
            fail(f"{rec.topic} has type {rec.type_name}, expected {expected}")
        raw_ns = stamp_ns(rec.message, "frame.header.stamp")
        if raw_ns is None:
            fail("typed CGI message has no frame.header.stamp")
        corrected = raw_ns + self.c["clock"]["reference_header_offset_ns"]
        valid = get(rec.message, "frame.valid") is True
        issue = ""
        if raw_ns <= 0 or corrected <= 0:
            issue = "invalid_timestamp"
        elif key in self.last_stamp and raw_ns <= self.last_stamp[key]:
            issue = "nonmonotonic_timestamp"
        if not valid:
            issue = "invalid_device_frame"
        self.last_stamp[key] = max(raw_ns, self.last_stamp.get(key, raw_ns))
        self.latest[key] = {"message": rec.message, "raw_ns": raw_ns, "stamp_ns": corrected,
                            "bag_ns": rec.bag_ns, "bag_common_ns": rec.bag_ns + self.c["clock"]["reference_bag_offset_ns"],
                            "valid": valid, "issue": issue}
        event_issue = issue or self.instant_quality_issue(key, rec.message)
        if event_issue:
            self.good_count = 0
            # A transient loss must break interpolation even if a later status
            # packet recovers before the next latitude anchor arrives.
            self.intervening_rejections.append(event_issue + ":" + key)
        if key != "latitude":
            return None
        return self.make_row(corrected)

    def instant_quality_issue(self, key, message):
        q = self.c["quality"]
        values = [get(message, attr, float("nan")) for attr in FIELDS[key][2]]
        if not all(isinstance(v, (int, float)) and math.isfinite(v) for v in values):
            return "nonfinite_field"
        if key == "latitude" and abs(values[0]) >= 90:
            return "invalid_navigation_value"
        if key == "longitude" and abs(values[0]) > 180:
            return "invalid_navigation_value"
        if key == "attitude" and not (0 <= values[0] < 360 and abs(values[1]) <= 90 and abs(values[2]) <= 180):
            return "invalid_navigation_value"
        if key == "status":
            state, sat, age = values
            if state != 2 or sat != 4:
                return "nonfixed_status"
            if not 0 <= age <= q["max_differential_age_sec"]:
                return "differential_age"
        if key.endswith("sigma"):
            if min(values) < 0:
                return "invalid_sigma"
            if key == "position_sigma" and (max(values[:2]) > q["max_horizontal_sigma_m"] or values[2] > q["max_vertical_sigma_m"]):
                return key
            if key == "attitude_sigma" and max(values) > q["max_attitude_sigma_deg"]:
                return key
            if key == "velocity_sigma" and max(values) > q["max_velocity_sigma_mps"]:
                return key
        return ""

    def make_row(self, anchor):
        c, quality = self.c, self.c["quality"]
        row = dict.fromkeys(CORE_COLUMNS, "")
        row.update(timestamp_ns=anchor, reference_valid=0, reason="", epoch=c["epoch"], time_basis=c["time_basis"])
        reasons = ["intervening_rejection:" + issue for issue in dict.fromkeys(self.intervening_rejections)]
        self.intervening_rejections.clear()
        stamps = []
        at = self.latest["latitude"]["bag_common_ns"]
        row["anchor_bag_timestamp_ns"] = self.latest["latitude"]["bag_ns"]
        row["reference_timestamp_semantics"] = "corrected_socketcan_receive_proxy"
        for key in FIELDS:
            f = self.latest.get(key)
            if f is None:
                reasons.append("missing_field:" + key)
                continue
            row[key + "_header_timestamp_ns"] = f["raw_ns"]
            row[key + "_corrected_timestamp_ns"] = f["stamp_ns"]
            row[key + "_bag_timestamp_ns"] = f["bag_ns"]
            row[key + "_frame_valid"] = int(f["valid"])
            row[key + "_age_ns"] = at - f["stamp_ns"]
            stamps.append(f["stamp_ns"])
            if f["issue"]:
                reasons.append(f["issue"] + ":" + key)
            if at - f["stamp_ns"] > round(quality["max_age_sec"] * NS) or at - f["bag_common_ns"] > round(quality["max_age_sec"] * NS):
                reasons.append("stale_field:" + key)
            if f["stamp_ns"] - at > round(quality["max_future_sec"] * NS) or f["bag_common_ns"] > at:
                reasons.append("clock_mismatch:" + key)
            if key in PRIMARY_FIELDS and f["stamp_ns"] <= self.consumed.get(key, 0):
                reasons.append("navigation_field_not_advanced:" + key)
            for name in FIELDS[key][2]:
                row[key + "_" + name] = get(f["message"], name, "")
        for key in PRIMARY_FIELDS:
            if key in self.latest:
                self.consumed[key] = self.latest[key]["stamp_ns"]
        row["field_skew_ns"] = max(stamps) - min(stamps) if stamps else ""
        if row["field_skew_ns"] != "" and row["field_skew_ns"] > round(quality["max_skew_sec"] * NS):
            reasons.append("field_time_skew")
        if len(self.latest) == len(FIELDS):
            values = {k: [get(self.latest[k]["message"], attr, float("nan")) for attr in FIELDS[k][2]] for k in FIELDS}
            state, sat, age = values["status"]
            if state != 2:
                reasons.append("not_integrated_navigation")
            if sat != 4:
                reasons.append("not_rtk_fixed_and_heading_valid")
            if not math.isfinite(age) or age < 0 or age > quality["max_differential_age_sec"]:
                reasons.append("differential_age")
            lat, lon, alt = values["latitude"][0], values["longitude"][0], values["altitude"][0]
            heading, pitch, roll = values["attitude"]
            nav = [lat, lon, alt, heading, pitch, roll] + values["velocity"]
            valid_nav = all(math.isfinite(v) for v in nav) and abs(lat) < 90 and abs(lon) <= 180 and 0 <= heading < 360 and abs(pitch) <= 90 and abs(roll) <= 180
            if not valid_nav:
                reasons.append("invalid_navigation_value")
            for key in ("position_sigma", "attitude_sigma", "velocity_sigma"):
                if not all(math.isfinite(v) and v >= 0 for v in values[key]):
                    reasons.append("invalid_sigma:" + key)
            if max(values["position_sigma"][:2]) > quality["max_horizontal_sigma_m"] or values["position_sigma"][2] > quality["max_vertical_sigma_m"]:
                reasons.append("position_sigma")
            if max(values["attitude_sigma"]) > quality["max_attitude_sigma_deg"]:
                reasons.append("attitude_sigma")
            if max(values["velocity_sigma"]) > quality["max_velocity_sigma_mps"]:
                reasons.append("velocity_sigma")
            if valid_nav:
                p, q = convert_pose([lat, lon, alt + c["height"]["geoid_separation_m"]], [heading, pitch, roll], c["georeference"]["origin_llh"])
                row.update(zip(["x", "y", "z", "qx", "qy", "qz", "qw"], [*p, *q]))
                if np.linalg.norm(p) > quality["max_site_radius_m"]:
                    reasons.append("outside_site_radius")
        if self.previous_anchor is not None and anchor <= self.previous_anchor:
            reasons.append("reference_clock_rollback" if anchor < self.previous_anchor else "duplicate_reference_timestamp")
        self.previous_anchor = max(anchor, self.previous_anchor or anchor)
        if reasons:
            self.good_count = 0
        else:
            self.good_count += 1
            if self.good_count < quality["recovery_samples"]:
                reasons.append("recovering")
        row["reference_valid"] = int(not reasons)
        row["reason"] = ";".join(dict.fromkeys(reasons)) or "good"
        return row


def reference_columns():
    columns = CORE_COLUMNS + ["anchor_bag_timestamp_ns", "reference_timestamp_semantics", "field_skew_ns"]
    for key, (_, _, attrs) in FIELDS.items():
        columns += [key + suffix for suffix in ("_header_timestamp_ns", "_corrected_timestamp_ns", "_bag_timestamp_ns", "_frame_valid", "_age_ns")]
        columns += [key + "_" + attr for attr in attrs]
    return columns


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def audit_clock():
    result = {"realtime_ns": time.time_ns(), "monotonic_ns": time.monotonic_ns(), "hostname": platform.node()}
    for name, command in {"chrony_tracking": ["chronyc", "tracking"], "chrony_sources": ["chronyc", "sources", "-v"],
                          "timedatectl": ["timedatectl", "show"], "kernel": ["uname", "-a"]}.items():
        try:
            proc = subprocess.run(command, capture_output=True, text=True, timeout=5, check=False)
            result[name] = {"returncode": proc.returncode, "stdout": proc.stdout, "stderr": proc.stderr}
        except (OSError, subprocess.TimeoutExpired) as exc:
            result[name] = {"unavailable": str(exc)}
    return result


def capture_plan(c, role, directory, config_path):
    prefix = c["capture"]["cgi_prefix"]
    reference = [prefix + "/" + s for s in RAW_SUFFIXES + ["driver_status", "navsat_fix", "twist", "imu/data_raw"]]
    cap = c["capture"]
    groups = {"reference": reference, "algorithm": cap["input_topics"] + cap["output_topics"],
              "combined": reference + cap["input_topics"] + cap["output_topics"],
              "outputs": cap["output_topics"], "playback": cap["output_topics"]}
    topics = list(dict.fromkeys(groups[role] + cap["evidence_topics"]))
    required = []
    if role in {"reference", "combined"}:
        required += [prefix + "/" + s for s in RAW_SUFFIXES + ["driver_status"]]
    if role in {"algorithm", "combined"}:
        required += cap["input_topics"]
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=False)
    metadata = {"schema_version": 1, "source_recording_id": c["source_recording_id"], "role": role,
                "epoch": c["epoch"], "time_basis": c.get("time_basis"), "hostname": platform.node(),
                "config_sha256": hashlib.sha256(Path(config_path).read_bytes()).hexdigest(),
                "config": c, "topics": topics, "required_topics": required,
                "input_topics": cap["input_topics"], "storage_id": cap["storage_id"],
                "environment": {key: os.environ.get(key) for key in ("ROS_DISTRO", "ROS_DOMAIN_ID", "ROS_LOCALHOST_ONLY", "RMW_IMPLEMENTATION")},
                "clock_start": audit_clock(), "capture_start_ns": None, "capture_end_ns": None,
                "raw_can_log": "can_raw.log" if role in {"reference", "combined"} else None}
    write_json(out / "capture.metadata.json", metadata)
    (out / "config.yaml").write_bytes(Path(config_path).read_bytes())
    for name, values in (("topics.txt", topics), ("required_topics.txt", required), ("input_topics.txt", cap["input_topics"])):
        (out / name).write_text("".join(v + "\n" for v in values), encoding="utf-8")
    qos = {topic: {"history": "keep_last", "depth": 4096, "reliability": "best_effort", "durability": "volatile"} for topic in topics}
    if "/tf_static" in qos:
        qos["/tf_static"]["durability"] = "transient_local"
    (out / "qos.yaml").write_text(yaml.safe_dump(qos, sort_keys=False), encoding="utf-8")
    print(json.dumps({"directory": str(out), "role": role, "topics": len(topics), "required": len(required)}))


def bag_records(path, mandatory_topics):
    import rosbag2_py
    from rclpy.serialization import deserialize_message
    from rosidl_runtime_py.utilities import get_message
    reader = rosbag2_py.SequentialReader()
    # Empty storage_id lets rosbag2 use metadata.yaml, including MCAP split files.
    reader.open(rosbag2_py.StorageOptions(uri=str(path), storage_id=""), rosbag2_py.ConverterOptions("", ""))
    types = {entry.name: entry.type for entry in reader.get_all_topics_and_types()}
    missing = set(mandatory_topics) - set(types)
    if missing:
        fail("missing required bag topics: " + ", ".join(sorted(missing)))
    classes = {}
    for topic, typ in types.items():
        try:
            classes[topic] = get_message(typ)
        except (ImportError, AttributeError, ModuleNotFoundError, ValueError):
            if topic in mandatory_topics:
                raise
            classes[topic] = None
    while reader.has_next():
        topic, serialized, timestamp = reader.read_next()
        cls = classes[topic]
        # Unavailable optional interfaces remain visible in timeline.csv; bytes
        # remain in the original immutable bag, not silently discarded.
        message = deserialize_message(serialized, cls) if cls else None
        yield Record(topic, types[topic], int(timestamp), message, len(serialized))


def plain_message(message):
    from rosidl_runtime_py.convert import message_to_ordereddict
    return message_to_ordereddict(message)


def estimate_pose(rec, c):
    if rec.type_name not in {"geometry_msgs/msg/PoseStamped", "nav_msgs/msg/Odometry", "geometry_msgs/msg/PoseWithCovarianceStamped"}:
        fail("estimated trajectory requires PoseStamped, PoseWithCovarianceStamped or Odometry")
    timestamp = stamp_ns(rec.message)
    if timestamp is None or timestamp <= 0:
        fail("estimated pose has invalid header stamp")
    pose = get(rec.message, "pose")
    pose = get(pose, "pose", pose)
    xyz = [float(get(pose, "position." + k, float("nan"))) for k in "xyz"]
    quat = [float(get(pose, "orientation." + k, float("nan"))) for k in "xyzw"]
    norm = np.linalg.norm(quat)
    if not all(math.isfinite(v) for v in xyz + quat) or abs(norm - 1) > 1e-3:
        fail("estimated pose has invalid position/quaternion; original bag retained")
    row = dict(zip(["timestamp_ns", "x", "y", "z", "qx", "qy", "qz", "qw"],
                   [timestamp + c["clock"]["estimate_header_offset_ns"], *xyz, *(np.array(quat) / norm)]))
    row.update(epoch=c["epoch"], time_basis=c["time_basis"], header_timestamp_ns=timestamp,
               bag_timestamp_ns=rec.bag_ns, header_frame_id=get(rec.message, "header.frame_id", ""))
    return row


def interval_metadata(c, role, bag_path, capture_path, min_bag, max_bag, sample_bounds=None):
    evidence = {"config": c, "bag_path": str(Path(bag_path).resolve()), "measurement_epoch_available": False,
                "assembly": "causal_record_order_at_latitude", "timestamp_offset_inference": "none"}
    evidence["applied_time_corrections"] = {"header_offset_ns": c["clock"][role + "_header_offset_ns"],
                                           "bag_offset_ns": c["clock"][role + "_bag_offset_ns"],
                                           "output_time_basis": c["time_basis"],
                                           "formula": "canonical_timestamp_ns = source_timestamp_ns + explicit_offset_ns"}
    offset = c["clock"][role + "_bag_offset_ns"]
    if capture_path:
        capture = json.loads(Path(capture_path).read_text(encoding="utf-8"))
        if capture["source_recording_id"] != c["source_recording_id"] or capture["epoch"] != c["epoch"]:
            fail("capture session/epoch does not match export configuration")
        start, end = capture.get("capture_start_ns"), capture.get("capture_end_ns")
        if not isinstance(start, int) or not isinstance(end, int) or start >= end:
            fail("capture metadata has no completed recording interval")
        if capture.get("exit_code") != 0:
            fail("capture did not finish cleanly; review recorder logs and create a reviewed interval before export")
        evidence["capture_metadata"] = capture
        evidence["interval_source"] = "capture_process_start_end"
        if capture.get("role") == "playback":
            source = capture.get("playback_source_interval", {})
            start, end = source.get("start_ns"), source.get("end_ns")
            if not isinstance(start, int) or not isinstance(end, int) or start >= end:
                fail("playback output capture requires original input bag's recorded time interval")
            evidence["interval_source"] = "original_playback_input_bag_bounds"
    else:
        start, end = min_bag, max_bag
        evidence["interval_source"] = "all_bag_messages_bounds; silence outside bag bounds is unknown"
    if start is None or end is None or start >= end:
        fail("bag has no nonempty time interval")
    start, end = start + offset, end + offset
    evidence["original_recording_boundary_common_ns"] = {"start_ns": start, "end_ns": end}
    if sample_bounds and sample_bounds[0] is not None:
        # A captured header can precede recorder receipt by transport latency.
        # Preserve the entire original recording boundary AND every canonical
        # timestamp (including rejected reference samples); never trim to fixes.
        early_margin = c["clock"]["uncertainty_ns"] + round(c["quality"]["max_age_sec"] * NS)
        late_margin = c["clock"]["uncertainty_ns"] + round(c["quality"]["max_future_sec"] * NS)
        evidence["header_boundary_margin_ns"] = {"early": early_margin, "late": late_margin,
                                                 "basis": "explicit clock uncertainty plus declared field age/future gates"}
        if sample_bounds[0] <= 0 or sample_bounds[1] <= 0:
            fail("invalid canonical clock timestamp <=0; raw/CSV timeline retained, recording interval was not expanded")
        if sample_bounds[0] < start - early_margin or sample_bounds[1] > end + late_margin:
            fail("canonical header extent lies outside the audited recording boundary margin; check explicit clock offsets/epochs; raw/CSV retained, interval was not expanded")
        evidence["all_canonical_sample_extent_ns"] = {"start_ns": sample_bounds[0], "end_ns": sample_bounds[1]}
        start, end = min(start, sample_bounds[0]), max(end, sample_bounds[1])
    conditions = ["reference_timestamp_is_socketcan_receive_proxy", "cgi_field_measurement_epoch_unverified"] if role == "reference" else []
    if not c["clock"].get("confirmed"):
        conditions.append("clock_alignment_assumed")
    if role == "reference" and not c["assembly"].get("confirmed"):
        conditions.append("causal_field_assembly_unverified")
    if not capture_path:
        conditions.append("capture_process_boundaries_unavailable_message_extent_only")
    meta = {"schema_version": 1, "source_recording_id": c["source_recording_id"],
            "world_frame_id": c["world_frame_id"] if role == "reference" else c["estimate"]["world_frame_id"],
            "body_frame_id": reference_body_frame(c["reference_point"]) if role == "reference" else "rear_axle_flu",
            "time_basis": c["time_basis"], "clock_epoch": c["epoch"], "conditions": conditions,
            "trajectory_kind": "reference" if role == "reference" else c["estimate"]["trajectory_kind"],
            "recorded_intervals": [{"start_ns": start, "end_ns": end, "epoch": c["epoch"]}],
            "evidence": evidence}
    if role == "reference":
        meta["reference_point"] = dict(c["reference_point"])
    return meta


def export_stream(records, c, out, role):
    """Stream complete timeline/raw/status evidence while keeping no cloud history."""
    assembler = ReferenceAssembler(c) if role == "reference" else None
    columns = reference_columns() if assembler else ["timestamp_ns", "x", "y", "z", "qx", "qy", "qz", "qw", "epoch", "time_basis", "header_timestamp_ns", "bag_timestamp_ns", "header_frame_id"]
    stem = "reference" if assembler else "estimated"
    counts, rejects = Counter(), Counter()
    low, high, rows, valid, last_estimate = None, None, 0, 0, None
    pose_low, pose_high = None, None
    status_topics = set(c["capture"]["output_topics"]) - {c["estimate"]["topic"]}
    status_topics |= {c["capture"]["cgi_prefix"] + "/driver_status", c["capture"]["cgi_prefix"] + "/time"}
    status_topics |= {c["capture"]["cgi_prefix"] + "/" + FIELDS[k][0] for k in ("status", "position_sigma", "attitude_sigma", "velocity_sigma")}
    with (out / (stem + ".csv")).open("w", newline="", encoding="utf-8") as poses, \
         (out / (stem + ".timeline.csv")).open("w", newline="", encoding="utf-8") as timeline, \
         (out / (stem + ".raw_frames.csv")).open("w", newline="", encoding="utf-8") as raw, \
         (out / (stem + ".status.jsonl")).open("w", encoding="utf-8") as status, \
         (out / (stem + ".tum")).open("w", encoding="utf-8") as tum:
        writer = csv.DictWriter(poses, fieldnames=columns)
        writer.writeheader()
        timeline_writer = csv.writer(timeline)
        timeline_writer.writerow(["topic", "type", "bag_timestamp_ns", "header_timestamp_ns", "comm_timestamp_ns", "serialized_bytes", "decoded", "epoch", "time_basis"])
        raw_writer = csv.writer(raw)
        raw_writer.writerow(["topic", "bag_timestamp_ns", "header_timestamp_ns", "can_id", "is_extended", "dlc", "valid", "data_hex", "epoch"])
        for rec in records:
            counts[rec.topic] += 1
            low = rec.bag_ns if low is None else min(low, rec.bag_ns)
            high = rec.bag_ns if high is None else max(high, rec.bag_ns)
            msg = rec.message
            hs = stamp_ns(msg, "frame.header.stamp") if get(msg, "frame") is not None else stamp_ns(msg)
            comm = get(msg, "comm_header.stamp_us")
            timeline_writer.writerow([rec.topic, rec.type_name, rec.bag_ns, hs, None if comm is None else int(comm) * 1000, rec.size, int(msg is not None), c["epoch"], c["time_basis"]])
            frame = get(msg, "frame")
            if frame is not None and rec.type_name.startswith("cgi430_interfaces/msg/"):
                data = get(frame, "data", [])
                raw_writer.writerow([rec.topic, rec.bag_ns, hs, get(frame, "can_id"), get(frame, "is_extended"), get(frame, "dlc"), get(frame, "valid"), bytes(data).hex(), c["epoch"]])
            if rec.topic in status_topics and msg is not None:
                payload = msg if isinstance(msg, dict) else plain_message(msg)
                status.write(json.dumps({"topic": rec.topic, "bag_timestamp_ns": rec.bag_ns, "header_timestamp_ns": hs,
                                         "epoch": c["epoch"], "time_basis": c["time_basis"], "message": payload}, default=lambda v: v.tolist() if hasattr(v, "tolist") else str(v)) + "\n")
            row = assembler.ingest(rec) if assembler else (estimate_pose(rec, c) if rec.topic == c["estimate"]["topic"] else None)
            if row is None:
                continue
            writer.writerow(row)
            rows += 1
            pose_low = row["timestamp_ns"] if pose_low is None else min(pose_low, row["timestamp_ns"])
            pose_high = row["timestamp_ns"] if pose_high is None else max(pose_high, row["timestamp_ns"])
            is_valid = row["reference_valid"] == 1 if assembler else True
            if not assembler and last_estimate is not None and row["timestamp_ns"] <= last_estimate:
                fail("estimated pose stamps repeat/rollback; split clock epochs explicitly, never deduplicate silently")
            if not assembler:
                last_estimate = row["timestamp_ns"]
            if is_valid:
                valid += 1
                stamp = int(row["timestamp_ns"])
                tum.write(f"{stamp // NS}.{stamp % NS:09d} " + " ".join(f"{row[k]:.12g}" for k in ["x", "y", "z", "qx", "qy", "qz", "qw"]) + "\n")
            elif assembler:
                rejects.update(row["reason"].split(";"))
    required_outputs = set(c["capture"]["required_output_topics"]) | {c["estimate"]["topic"]}
    missing_outputs = sorted(required_outputs - set(counts)) if not assembler else []
    return {"start_ns": low, "end_ns": high, "rows": rows, "valid_rows": valid,
            "all_pose_start_ns": pose_low, "all_pose_end_ns": pose_high,
            "topics": dict(counts), "rejection_reasons": dict(rejects),
            "missing_output_topics": missing_outputs,
            "satisfied": bool(valid) and not missing_outputs,
            "satisfied_meaning": "usable reference or required output streams present; not an accuracy acceptance decision"}


def export_bags(args, c):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    write_json(out / "export.config.json", c)
    summaries = {}
    for role, bag, capture in (("reference", args.reference_bag, args.reference_capture), ("estimate", args.estimate_bag, args.estimate_capture)):
        if not bag:
            continue
        # Complete algorithm failure may never create its pose topic. Preserve
        # status/timeline and an empty TUM with the full recording denominator.
        required = [c["capture"]["cgi_prefix"] + "/" + item[0] for item in FIELDS.values()] if role == "reference" else []
        summary = export_stream(bag_records(bag, required), c, out, role)
        try:
            meta = interval_metadata(c, role, bag, capture, summary["start_ns"], summary["end_ns"],
                                     (summary["all_pose_start_ns"], summary["all_pose_end_ns"]))
        except ValueError as exc:
            write_json(out / "export.failure.json", {"role": role, "error": str(exc), "summary": summary,
                                                     "reference_and_raw_timeline_retained": True})
            raise
        meta["evidence"]["export_summary"] = summary
        if role == "estimate" and summary["rows"] == 0:
            meta["conditions"].append("ESTIMATE_STREAM_ABSENT")
        elif role == "estimate" and summary["missing_output_topics"]:
            meta["conditions"].append("REQUIRED_STATUS_STREAMS_ABSENT")
        stem = "reference" if role == "reference" else "estimated"
        write_json(out / (stem + ".metadata.json"), meta)
        # A mapping IMU TUM needs the capture interval/time evidence without
        # inheriting the published rear-axle pose's frame declaration.
        recording_meta = {key: meta[key] for key in ("schema_version", "source_recording_id", "time_basis", "clock_epoch",
                                                     "recorded_intervals", "conditions", "evidence")}
        write_json(out / (stem + ".recording.metadata.json"), recording_meta)
        summaries[role] = summary
    write_json(out / "export.summary.json", summaries)
    print(json.dumps(summaries, ensure_ascii=False, indent=2))
    if summaries["reference"]["rejection_reasons"].get("reference_clock_rollback", 0):
        fail("reference clock rolled back within epoch; raw audit exported, split epochs explicitly before evaluation")
    # Missing estimates are a real zero-availability result. Keep the command
    # chain usable so the evaluator reports NO_PAIRED_DATA with the full interval.
    # Return 0 denotes export completion, never algorithm quality or acceptance.
    return 0 if summaries["reference"]["valid_rows"] else 3


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check-config", help="fail closed until experiment contracts are confirmed")
    check.add_argument("--config", required=True)
    check.add_argument("--purpose", choices=["capture", "export"], default="export")
    check.add_argument("--role", choices=["reference", "algorithm", "combined", "outputs", "playback"])
    plan = sub.add_parser("plan", help="write capture plan, QoS and initial local session metadata")
    plan.add_argument("--config", required=True)
    plan.add_argument("--role", required=True, choices=["reference", "algorithm", "combined", "outputs", "playback"])
    plan.add_argument("--output", required=True)
    exp = sub.add_parser("export", help="export CGI reference, optional published estimate and full evidence timelines")
    exp.add_argument("--config", required=True)
    exp.add_argument("--reference-bag", required=True)
    exp.add_argument("--reference-capture", help="reference computer's capture.metadata.json")
    exp.add_argument("--estimate-bag")
    exp.add_argument("--estimate-capture", help="estimate computer's capture.metadata.json")
    exp.add_argument("--output", required=True, help="new directory; existing results are never overwritten")
    mark = sub.add_parser("mark", help="record local capture start/end and shutdown result")
    mark.add_argument("--metadata", required=True)
    mark.add_argument("--event", required=True, choices=["start", "end"])
    mark.add_argument("--exit-code", type=int, default=0)
    clocks = sub.add_parser("clock-log", help="write realtime/monotonic pairs; stop with SIGINT/TERM")
    clocks.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "clock-log":
            with Path(args.output).open("x", encoding="utf-8") as stream:
                while True:
                    stream.write(json.dumps({"realtime_ns": time.time_ns(), "monotonic_ns": time.monotonic_ns()}) + "\n")
                    stream.flush()
                    time.sleep(1)
        elif args.command == "mark":
            path = Path(args.metadata)
            metadata = json.loads(path.read_text(encoding="utf-8"))
            metadata["capture_" + args.event + "_ns"] = time.time_ns()
            metadata["capture_" + args.event + "_monotonic_ns"] = time.monotonic_ns()
            if args.event == "end":
                metadata["exit_code"] = args.exit_code
                metadata["clock_end"] = audit_clock()
                statuses = path.parent / "logs/child_exit_status.tsv"
                if statuses.exists():
                    metadata["child_exit_codes"] = {name: int(value) for name, value in
                                                    (line.split("\t") for line in statuses.read_text().splitlines())}
            write_json(path, metadata)
        else:
            purpose = "capture" if args.command == "plan" else getattr(args, "purpose", "export")
            c = load_config(args.config, purpose, getattr(args, "role", None), bool(getattr(args, "estimate_bag", None)))
            if args.command == "plan":
                capture_plan(c, args.role, args.output, args.config)
            elif args.command == "export":
                return export_bags(args, c)
            else:
                print(f"{purpose} contract checked; no clock correction was inferred")
        return 0
    except KeyboardInterrupt:
        return 0
    except (ValueError, OSError, ImportError, KeyError, RuntimeError) as exc:
        print("accuracy_export_bag: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
