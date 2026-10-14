#!/usr/bin/env python3
"""Offline, frozen-frame SE(3) reference disagreement assessment.

TUM poses and canonical reference CSV must already represent rear_axle_flu.
No lever arm, time offset, scale, or evaluation-set alignment is estimated here.
See templates for the explicit data, frame, clock, and recording contracts.
"""
from __future__ import annotations

import argparse
import bisect
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, ROUND_HALF_EVEN
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
from typing import Any

import numpy as np
import yaml

VERSION = "1.0.0"
NS = 1_000_000_000
BODY_FRAME = "rear_axle_flu"
REFERENCE_COLUMNS = {
    "timestamp_ns", "x", "y", "z", "qx", "qy", "qz", "qw",
    "reference_valid", "reason", "epoch", "time_basis",
}
KINDS = {"frontend", "optimized_keyframes", "optimized_dense_reconstructed", "localization"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_document(path: Path) -> dict:
    with path.open(encoding="utf-8-sig") as stream:
        value = yaml.safe_load(stream)
    require(isinstance(value, dict), f"{path}: expected a YAML/JSON object")
    return value


def integer(value: Any, label: str) -> int:
    require(not isinstance(value, bool), f"{label}: expected integer nanoseconds")
    if isinstance(value, int):
        return value
    require(isinstance(value, str) and value.lstrip("-").isdigit(), f"{label}: expected integer nanoseconds")
    return int(value)


def positive(value: Any, label: str) -> float:
    value = float(value)
    require(math.isfinite(value) and value > 0, f"{label} must be finite and positive")
    return value


def seconds_ns(value: Any, label: str) -> int:
    result = int(round(positive(value, label) * NS))
    require(result > 0, f"{label} must be at least 1 ns")
    return result


def unit_quaternion(values: Any, label: str) -> np.ndarray:
    q = np.asarray(values, dtype=float)
    require(q.shape == (4,) and np.isfinite(q).all(), f"{label}: quaternion must be finite xyzw")
    norm = np.linalg.norm(q)
    require(abs(norm - 1.0) <= 1e-5, f"{label}: quaternion is not unit length")
    return q / norm


def quat_matrix(q: np.ndarray) -> np.ndarray:
    x, y, z, w = q
    return np.array([
        [1 - 2*(y*y + z*z), 2*(x*y - z*w), 2*(x*z + y*w)],
        [2*(x*y + z*w), 1 - 2*(x*x + z*z), 2*(y*z - x*w)],
        [2*(x*z - y*w), 2*(y*z + x*w), 1 - 2*(x*x + y*y)],
    ])


def matrix_quat(r: np.ndarray) -> np.ndarray:
    # Eigenvector form is well conditioned near rotations of pi.
    k = np.array([
        [r[0, 0]-r[1, 1]-r[2, 2], r[0, 1]+r[1, 0], r[0, 2]+r[2, 0], r[2, 1]-r[1, 2]],
        [r[0, 1]+r[1, 0], r[1, 1]-r[0, 0]-r[2, 2], r[1, 2]+r[2, 1], r[0, 2]-r[2, 0]],
        [r[0, 2]+r[2, 0], r[1, 2]+r[2, 1], r[2, 2]-r[0, 0]-r[1, 1], r[1, 0]-r[0, 1]],
        [r[2, 1]-r[1, 2], r[0, 2]-r[2, 0], r[1, 0]-r[0, 1], np.trace(r)],
    ]) / 3
    _, vectors = np.linalg.eigh(k)
    q = vectors[:, -1]
    return q if q[3] >= 0 else -q


def slerp(a: np.ndarray, b: np.ndarray, ratio: float) -> np.ndarray:
    dot = float(a @ b)
    if dot < 0:
        b, dot = -b, -dot
    if dot > 0.9995:
        result = (1-ratio)*a + ratio*b
        return result / np.linalg.norm(result)
    theta = math.acos(np.clip(dot, -1, 1))
    return (math.sin((1-ratio)*theta)*a + math.sin(ratio*theta)*b) / math.sin(theta)


def euler_xyz(r: np.ndarray) -> np.ndarray:
    # Fixed-axis XYZ / Rz(yaw) Ry(pitch) Rx(roll); diagnostic near gimbal lock.
    pitch = math.asin(np.clip(-r[2, 0], -1, 1))
    if abs(math.cos(pitch)) > 1e-8:
        roll, yaw = math.atan2(r[2, 1], r[2, 2]), math.atan2(r[1, 0], r[0, 0])
    else:
        roll, yaw = 0.0, math.atan2(-r[0, 1], r[1, 1])
    return np.degrees([roll, pitch, yaw])


def rotation_angle_deg(r: np.ndarray) -> float:
    # atan2 remains accurate for very small angles and includes pi rotations.
    skew = np.array([r[2, 1]-r[1, 2], r[0, 2]-r[2, 0], r[1, 0]-r[0, 1]]) / 2
    return math.degrees(math.atan2(np.linalg.norm(skew), np.clip((np.trace(r)-1)/2, -1, 1)))


def wrap_degrees(value: np.ndarray) -> np.ndarray:
    return (value + 180.0) % 360.0 - 180.0


def parse_intervals(value: Any, label: str) -> list[dict]:
    require(isinstance(value, list) and bool(value), f"{label}: require recorded/analysis intervals")
    result, last_end = [], {}
    for raw in value:
        epoch = str(raw["epoch"])
        require(bool(epoch), f"{label}: empty epoch")
        start, end = integer(raw["start_ns"], label), integer(raw["end_ns"], label)
        require(end > start, f"{label}: interval must have positive duration")
        require(epoch not in last_end or start >= last_end[epoch], f"{label}: overlapping/unsorted intervals")
        result.append({"start_ns": start, "end_ns": end, "epoch": epoch})
        last_end[epoch] = end
    return result


def interval_for(intervals: list[dict], epoch: str, timestamp: int, endpoint: bool = False) -> int | None:
    for i, item in enumerate(intervals):
        if item["epoch"] == epoch and item["start_ns"] <= timestamp and (
                timestamp < item["end_ns"] or (endpoint and timestamp == item["end_ns"])):
            return i
    return None


def subset_intervals(small: list[dict], large: list[dict]) -> bool:
    return all(any(a["epoch"] == b["epoch"] and b["start_ns"] <= a["start_ns"]
                   and a["end_ns"] <= b["end_ns"] for b in large) for a in small)


def intervals_overlap(a: list[dict], b: list[dict]) -> bool:
    return any(x["epoch"] == y["epoch"] and max(x["start_ns"], y["start_ns"]) < min(x["end_ns"], y["end_ns"])
               for x in a for y in b)


def shifted(intervals: list[dict], offset: int, epoch: str | None = None) -> list[dict]:
    return [{"start_ns": i["start_ns"]+offset, "end_ns": i["end_ns"]+offset,
             "epoch": epoch if epoch is not None else i["epoch"]} for i in intervals]


def read_metadata(path: Path, reference: bool) -> dict:
    meta = read_document(path)
    require(meta.get("schema_version") == 1, f"{path}: unsupported metadata schema")
    for key in ("source_recording_id", "world_frame_id", "body_frame_id", "time_basis"):
        require(isinstance(meta.get(key), str) and bool(meta[key].strip()), f"{path}: missing {key}")
    require(meta["body_frame_id"] == BODY_FRAME,
            f"{path}: pose must already be {BODY_FRAME}; apply the measured lever arm/extrinsic in the exporter once")
    require(isinstance(meta.get("conditions", []), list) and all(isinstance(c, str) for c in meta.get("conditions", [])),
            f"{path}: conditions must be a list of evidence/limitation strings")
    meta["conditions"] = list(meta.get("conditions", []))
    if any(word in meta["time_basis"].lower() for word in ("receive", "arrival", "proxy")):
        meta["conditions"].append("Timestamp is a receive/arrival proxy; measurement epoch and transport latency require independent evidence.")
    meta["recorded_intervals"] = parse_intervals(meta.get("recorded_intervals"), str(path))
    if not reference:
        require(meta.get("trajectory_kind") in KINDS, f"{path}: unsupported/missing trajectory_kind")
        require(meta["world_frame_id"].startswith(("map:", "run:")),
                f"{path}: world_frame_id must identify map:<asset hash> or run:<unique run id>")
        epochs = {i["epoch"] for i in meta["recorded_intervals"]}
        require(len(epochs) == 1, f"{path}: split TUM by clock epoch; TUM has no epoch column")
        meta["clock_epoch"] = str(meta.get("clock_epoch", next(iter(epochs))))
        require(epochs == {meta["clock_epoch"]}, f"{path}: clock_epoch disagrees with recording intervals")
    return meta


@dataclass
class Pose:
    timestamp: int
    position: np.ndarray
    quaternion: np.ndarray
    epoch: str
    valid: bool = True
    reason: str = ""


def read_reference(path: Path, meta: dict) -> list[Pose]:
    result, previous, previous_valid, closed_epochs, active_epoch = [], {}, {}, set(), None
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        require(REFERENCE_COLUMNS <= set(reader.fieldnames or []), f"{path}: missing canonical reference columns")
        for line, row in enumerate(reader, 2):
            label = f"{path}:{line}"
            timestamp, epoch = integer(row["timestamp_ns"], label), row["epoch"]
            require(row["time_basis"] == meta["time_basis"], f"{label}: time_basis differs from metadata")
            require(epoch and epoch not in closed_epochs, f"{label}: repeated/non-contiguous epoch")
            if active_epoch is not None and epoch != active_epoch:
                closed_epochs.add(active_epoch)
            active_epoch = epoch
            require(interval_for(meta["recorded_intervals"], epoch, timestamp, True) is not None,
                    f"{label}: reference row outside recorded intervals")
            require(row["reference_valid"].lower() in {"0", "1", "true", "false"}, f"{label}: invalid validity flag")
            valid = row["reference_valid"].lower() in {"1", "true"}
            require(epoch not in previous or timestamp >= previous[epoch], f"{label}: backward time requires a new epoch")
            require(epoch not in previous or timestamp != previous[epoch] or not (valid and previous_valid[epoch]),
                    f"{label}: duplicate valid reference time")
            previous[epoch], previous_valid[epoch] = timestamp, valid
            require(valid or bool(row["reason"].strip()), f"{label}: invalid reference needs a reason")
            if valid:
                position = np.array([float(row[k]) for k in ("x", "y", "z")])
                require(np.isfinite(position).all(), f"{label}: valid position must be finite")
                q = unit_quaternion([float(row[k]) for k in ("qx", "qy", "qz", "qw")], label)
            else:
                position, q = np.full(3, np.nan), np.full(4, np.nan)
            result.append(Pose(timestamp, position, q, epoch, valid, row["reason"]))
    return result


def read_tum(path: Path, meta: dict) -> list[Pose]:
    result = []
    with path.open(encoding="utf-8-sig") as stream:
        for line, raw in enumerate(stream, 1):
            fields = raw.partition("#")[0].split()
            if not fields:
                continue
            label = f"{path}:{line}"
            require(len(fields) == 8, f"{label}: expected t x y z qx qy qz qw")
            try:
                dec = Decimal(fields[0]) * NS
                require(dec.is_finite(), f"{label}: nonfinite time")
                timestamp = int(dec.to_integral_value(rounding=ROUND_HALF_EVEN))
            except InvalidOperation as exc:
                raise ValueError(f"{label}: invalid TUM timestamp") from exc
            require(not result or timestamp > result[-1].timestamp, f"{label}: duplicate/backward TUM time; split epochs explicitly")
            position = np.array([float(v) for v in fields[1:4]])
            require(np.isfinite(position).all(), f"{label}: nonfinite estimate position")
            q = unit_quaternion([float(v) for v in fields[4:]], label)
            epoch = meta["clock_epoch"]
            require(interval_for(meta["recorded_intervals"], epoch, timestamp, True) is not None,
                    f"{label}: estimate row outside recorded intervals")
            result.append(Pose(timestamp, position, q, epoch))
    return result


class Reference:
    def __init__(self, rows: list[Pose], intervals: list[dict], max_gap_ns: int):
        self.rows, self.intervals, self.max_gap_ns = rows, intervals, max_gap_ns
        self.epochs: dict[str, list[Pose]] = {}
        for row in rows:
            segment = self.epochs.setdefault(row.epoch, [])
            if segment and segment[-1].timestamp == row.timestamp:
                old = segment[-1]
                if not old.valid or not row.valid:
                    reason = ";".join(r.reason for r in (old, row) if not r.valid)
                    segment[-1] = Pose(row.timestamp, np.full(3, np.nan), np.full(4, np.nan), row.epoch, False, reason)
            else:
                segment.append(row)
        self.times = {e: [r.timestamp for r in rs] for e, rs in self.epochs.items()}

    def at(self, epoch: str, timestamp: int) -> tuple[Pose | None, str]:
        rows = self.epochs.get(epoch, [])
        index = bisect.bisect_left(self.times.get(epoch, []), timestamp)
        if index < len(rows) and rows[index].timestamp == timestamp:
            stop = bisect.bisect_right(self.times[epoch], timestamp)
            invalid = [r for r in rows[index:stop] if not r.valid]
            if invalid:
                return None, "reference_invalid:" + ";".join(r.reason for r in invalid)
            row = rows[index]
            return (row, "") if row.valid else (None, "reference_invalid:" + row.reason)
        if index == 0 or index == len(rows):
            return None, "reference_outside_samples"
        a, b = rows[index-1:index+1]
        if not (a.valid and b.valid):
            return None, "reference_invalid_gap:" + ";".join(r.reason for r in (a, b) if not r.valid)
        if b.timestamp-a.timestamp > self.max_gap_ns:
            return None, "reference_gap_exceeded"
        if interval_for(self.intervals, epoch, a.timestamp) != interval_for(self.intervals, epoch, b.timestamp, True):
            return None, "reference_recording_interruption"
        ratio = (timestamp-a.timestamp) / (b.timestamp-a.timestamp)
        return Pose(timestamp, (1-ratio)*a.position+ratio*b.position,
                    slerp(a.quaternion, b.quaternion, ratio), epoch), ""

    def valid_intervals(self) -> list[dict]:
        intervals = []
        for epoch, rows in self.epochs.items():
            for a, b in zip(rows, rows[1:]):
                if (a.valid and b.valid and 0 < b.timestamp-a.timestamp <= self.max_gap_ns
                        and interval_for(self.intervals, epoch, a.timestamp)
                        == interval_for(self.intervals, epoch, b.timestamp, True)):
                    intervals.append({"start_ns": a.timestamp, "end_ns": b.timestamp, "epoch": epoch})
        return merge_intervals(intervals)


def merge_intervals(intervals: list[dict]) -> list[dict]:
    result = []
    for value in sorted(intervals, key=lambda i: (i["epoch"], i["start_ns"])):
        item = dict(value)
        if result and result[-1]["epoch"] == item["epoch"] and item["start_ns"] <= result[-1]["end_ns"]:
            result[-1]["end_ns"] = max(result[-1]["end_ns"], item["end_ns"])
        else:
            result.append(item)
    return result


def estimate_availability(rows: list[Pose], recorded: list[dict], hold_ns: int, offset: int, reference_epoch: str) -> list[dict]:
    intervals = []
    for row in rows:
        index = interval_for(recorded, row.epoch, row.timestamp)
        if index is not None:
            end = min(row.timestamp+hold_ns, recorded[index]["end_ns"])
            if end > row.timestamp:
                intervals.append({"start_ns": row.timestamp+offset, "end_ns": end+offset, "epoch": reference_epoch})
    return merge_intervals(intervals)


def coverage(analysis: list[dict], reference_valid: list[dict], available: list[dict]) -> tuple[dict, list[dict]]:
    segments = []
    totals = {"recorded": 0, "reference_valid": 0, "estimate_available": 0, "joint": 0}
    # Binary-search disjoint spans, so frequent packet gaps do not cause a
    # quadratic scan on long recordings.
    reference_valid, available = merge_intervals(reference_valid), merge_intervals(available)
    indexes = []
    for intervals in (reference_valid, available):
        epochs = {}
        for span in intervals:
            starts, ends = epochs.setdefault(span["epoch"], ([], []))
            starts.append(span["start_ns"])
            ends.append(span["end_ns"])
        indexes.append(epochs)

    def contains(index, epoch, timestamp):
        starts, ends = index.get(epoch, ([], []))
        offset = bisect.bisect_right(starts, timestamp)-1
        return offset >= 0 and timestamp < ends[offset]

    for interval in analysis:
        start, end, epoch = interval["start_ns"], interval["end_ns"], interval["epoch"]
        boundaries = {start, end}
        for row in reference_valid + available:
            if row["epoch"] == epoch:
                boundaries.update(t for t in (row["start_ns"], row["end_ns"]) if start < t < end)
        times = sorted(boundaries)
        for a, b in zip(times, times[1:]):
            # Use the left endpoint of half-open spans; no floating-point epoch arithmetic.
            valid = contains(indexes[0], epoch, a)
            present = contains(indexes[1], epoch, a)
            duration = b-a
            totals["recorded"] += duration
            totals["reference_valid"] += duration * valid
            totals["estimate_available"] += duration * present
            totals["joint"] += duration * (valid and present)
            state = ("reference_valid" if valid else "reference_invalid") + ("_estimate_available" if present else "_estimate_unavailable")
            if segments and segments[-1]["epoch"] == epoch and segments[-1]["end_ns"] == a and segments[-1]["state"] == state:
                segments[-1]["end_ns"] = b
            else:
                segments.append({"epoch": epoch, "start_ns": a, "end_ns": b, "state": state})
    require(totals["recorded"] > 0, "zero recorded analysis duration")
    result = {key+"_duration_sec": value/NS for key, value in totals.items()}
    result.update({
        "reference_valid_fraction": totals["reference_valid"]/totals["recorded"],
        "estimate_availability_fraction": totals["estimate_available"]/totals["recorded"],
        "joint_coverage_fraction": totals["joint"]/totals["recorded"],
        "estimate_availability_given_valid_reference": (totals["joint"]/totals["reference_valid"] if totals["reference_valid"] else None),
        "estimate_unavailable_duration_sec": (totals["recorded"]-totals["estimate_available"])/NS,
        "reference_invalid_duration_sec": (totals["recorded"]-totals["reference_valid"])/NS,
        "denominator": "Declared recorded analysis duration; exact half-open interval union, including absent estimates.",
    })
    for segment in segments:
        segment["duration_sec"] = (segment["end_ns"]-segment["start_ns"])/NS
    return result, segments


def statistics(values: list[float] | np.ndarray, signed: bool = False) -> dict:
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return {k: None for k in ("rmse", "mean", "median", "p95", "max", "bias", "p95_abs", "max_abs")}
    return {
        "rmse": float(np.sqrt(np.mean(values**2))), "mean": float(np.mean(values)),
        "median": float(np.median(values)), "p95": float(np.percentile(values, 95)),
        "max": float(np.max(values)), "bias": float(np.mean(values)) if signed else None,
        "p95_abs": float(np.percentile(np.abs(values), 95)), "max_abs": float(np.max(np.abs(values))),
    }


def disagreement(estimate: Pose, reference: Pose, rotation: np.ndarray, translation: np.ndarray) -> dict:
    p = rotation @ estimate.position + translation
    r_est, r_ref = rotation @ quat_matrix(estimate.quaternion), quat_matrix(reference.quaternion)
    delta = p-reference.position
    body = r_ref.T @ delta
    rpy_est, rpy_ref = euler_xyz(r_est), euler_xyz(r_ref)
    rpy = wrap_degrees(rpy_est-rpy_ref)
    result = dict(zip(("x_m", "y_m", "z_m"), map(float, delta)))
    result.update(dict(zip(("longitudinal_m", "lateral_m", "vertical_m"), map(float, body))))
    result.update(dict(zip(("roll_deg", "pitch_deg", "yaw_deg"), map(float, rpy))))
    result.update({
        "horizontal_m": float(np.linalg.norm(delta[:2])), "position_3d_m": float(np.linalg.norm(delta)),
        "attitude_geodesic_deg": rotation_angle_deg(r_ref.T @ r_est),
        "euler_near_gimbal_lock": bool(abs(rpy_est[1]) > 89 or abs(rpy_ref[1]) > 89),
        "estimate_enu_x": float(p[0]), "estimate_enu_y": float(p[1]), "estimate_enu_z": float(p[2]),
        "reference_x": float(reference.position[0]), "reference_y": float(reference.position[1]), "reference_z": float(reference.position[2]),
    })
    return result


METRICS = ["x_m", "y_m", "z_m", "horizontal_m", "position_3d_m", "longitudinal_m", "lateral_m", "vertical_m",
           "roll_deg", "pitch_deg", "yaw_deg", "attitude_geodesic_deg"]
UNSIGNED = {"horizontal_m", "position_3d_m", "attitude_geodesic_deg"}


def summarize(errors: list[dict]) -> dict:
    return {key: statistics([e[key] for e in errors], signed=key not in UNSIGNED) for key in METRICS}


def file_contract(path: Path, meta_path: Path, meta: dict, rows: list[Pose]) -> dict:
    canonical = [[r.timestamp, r.epoch, bool(r.valid), r.reason,
                  [format(v, ".14g") for v in r.position] if r.valid else None,
                  [format(v, ".14g") for v in quat_matrix(r.quaternion).ravel()] if r.valid else None] for r in rows]
    canonical_hash = hashlib.sha256(json.dumps(canonical, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    return {"file": str(path), "sha256": sha256(path), "metadata_file": str(meta_path),
            "metadata_sha256": sha256(meta_path), "canonical_pose_sha256": canonical_hash, "metadata": meta}


def resolve(base: Path, value: str) -> Path:
    require(isinstance(value, str) and bool(value), "expected nonempty file path")
    return (base / value).resolve()


def load_config(path: Path) -> dict:
    config = read_document(path)
    require(config.get("schema_version") == 1, "unsupported configuration schema")
    require(config.get("role") in {"calibration", "evaluation"}, "role must be calibration or evaluation")
    config["analysis_intervals"] = parse_intervals(config.get("analysis_intervals"), "analysis_intervals")
    pairing = config.get("pairing", {})
    for key in ("max_reference_gap_sec", "availability_hold_sec"):
        seconds_ns(pairing.get(key), key)
    require(isinstance(config.get("estimates"), list) and config["estimates"], "at least one estimate is required")
    names = [e.get("name") for e in config["estimates"]]
    require(all(isinstance(n, str) and n and all(c.isalnum() or c in "_-" for c in n) for n in names), "estimate names must be file-safe")
    require(len(set(names)) == len(names), "estimate names must be unique")
    groups = config.get("common_sample_groups", [])
    require(isinstance(groups, list), "common_sample_groups must be a list")
    group_names = set()
    for group in groups:
        name, members = group.get("name"), group.get("estimates")
        require(isinstance(name, str) and name and all(c.isalnum() or c in "_-" for c in name) and name not in group_names,
                "common sample group names must be unique and file-safe")
        require(isinstance(members, list) and len(members) >= 2 and len(set(members)) == len(members)
                and set(members) <= set(names), "common sample group must select at least two distinct configured estimates")
        group_names.add(name)
    require(isinstance(config.get("conditions", []), list) and all(isinstance(c, str) for c in config.get("conditions", [])),
            "conditions must be a list of evidence/limitation strings")
    config["_path"], config["_base"] = path.resolve(), path.resolve().parent
    return config


def load_inputs(config: dict) -> tuple[Reference, dict, list[dict]]:
    base, source = config["_base"], config["reference"]
    path, meta_path = resolve(base, source["csv"]), resolve(base, source["metadata"])
    meta = read_metadata(meta_path, reference=True)
    require(subset_intervals(config["analysis_intervals"], meta["recorded_intervals"]), "analysis interval extends beyond reference recording metadata")
    ref = Reference(read_reference(path, meta), meta["recorded_intervals"], seconds_ns(config["pairing"]["max_reference_gap_sec"], "max_reference_gap_sec"))
    contract = file_contract(path, meta_path, meta, ref.rows)
    estimates = []
    for item in config["estimates"]:
        path, meta_path = resolve(base, item["tum"]), resolve(base, item["metadata"])
        emeta = read_metadata(meta_path, reference=False)
        timing = {**config.get("time", {}), **item.get("time", {})}
        status = timing.get("offset_status")
        require(status in {"known", "assumed"}, "unknown time offset: require measured known offset, or explicit assumed offset with allow_conditional")
        require(status == "known" or timing.get("allow_conditional") is True, "assumed time offset requires allow_conditional: true")
        require(isinstance(timing.get("evidence"), str) and bool(timing["evidence"].strip()), "time evidence is required")
        offset = integer(timing.get("estimate_to_reference_offset_ns"), "estimate_to_reference_offset_ns")
        reference_epoch = str(timing.get("reference_epoch", emeta["clock_epoch"]))
        require(reference_epoch in {i["epoch"] for i in config["analysis_intervals"]}, "estimate epoch mapping is outside analysis intervals")
        condition = list(meta["conditions"]) + list(emeta["conditions"])
        if status == "assumed":
            condition.append("Time offset is assumed; dynamic pose disagreement is conditional on that assumption.")
        condition = list(dict.fromkeys(condition))
        rows = read_tum(path, emeta)
        estimates.append({"name": item["name"], "rows": rows, "contract": file_contract(path, meta_path, emeta, rows),
                          "metadata": emeta, "offset": offset, "reference_epoch": reference_epoch,
                          "hold_ns": seconds_ns(item.get("availability_hold_sec", config["pairing"]["availability_hold_sec"]), "availability_hold_sec"),
                          "time": timing, "conditions": condition, "alignment_path": resolve(base, item["alignment"]) if item.get("alignment") else None})
    return ref, contract, estimates


def matched_rows(config: dict, reference: Reference, estimate: dict) -> tuple[list[tuple[Pose, Pose]], list[dict]]:
    pairs, unmatched = [], []
    for source in estimate["rows"]:
        timestamp, epoch = source.timestamp+estimate["offset"], estimate["reference_epoch"]
        if interval_for(config["analysis_intervals"], epoch, timestamp) is None:
            continue
        target, reason = reference.at(epoch, timestamp)
        if target is None:
            unmatched.append({"timestamp_ns": timestamp, "estimate_timestamp_ns": source.timestamp,
                              "epoch": epoch, "status": reason})
        else:
            pairs.append((source, target))
    return pairs, unmatched


def fit_rigid(pairs: list[tuple[Pose, Pose]], options: dict) -> tuple[np.ndarray, np.ndarray, dict]:
    minimum = integer(options.get("min_pairs"), "fit.min_pairs")
    require(minimum >= 3 and len(pairs) >= minimum, f"calibration requires at least {minimum} matched poses")
    x, y = np.array([a.position for a, _ in pairs]), np.array([b.position for _, b in pairs])
    xc, yc = x-x.mean(axis=0), y-y.mean(axis=0)
    spans = [np.linalg.svd(v, compute_uv=False)/math.sqrt(len(v)) for v in (xc, yc)]
    min_span = positive(options.get("min_secondary_span_m"), "fit.min_secondary_span_m")
    require(min(spans[0][1], spans[1][1]) >= min_span,
            "calibration path is stationary/near-collinear; collect a separate route with turns")
    h = xc.T @ yc
    method = options.get("method")
    require(method in {"position", "position_attitude"}, "fit.method must be position or position_attitude")
    if method == "position_attitude":
        weight = positive(options.get("attitude_weight_m"), "fit.attitude_weight_m")
        for a, b in pairs:
            h += weight**2 * quat_matrix(a.quaternion) @ quat_matrix(b.quaternion).T
    u, _, vt = np.linalg.svd(h)
    d = np.eye(3)
    d[2, 2] = np.linalg.det(vt.T @ u.T)
    rotation = vt.T @ d @ u.T
    translation = y.mean(axis=0) - rotation @ x.mean(axis=0)
    if method == "position_attitude":
        limit = positive(options.get("max_attitude_dispersion_deg"), "fit.max_attitude_dispersion_deg")
        require(limit <= 180, "max_attitude_dispersion_deg must be <= 180")
        angles = [rotation_angle_deg(quat_matrix(b.quaternion).T @ rotation @ quat_matrix(a.quaternion)) for a, b in pairs]
        require(max(angles) <= limit, "calibration attitudes violate configured dispersion gate")
    return rotation, translation, {"source_path_rms_spans_m": spans[0].tolist(), "reference_path_rms_spans_m": spans[1].tolist(),
                                   "options": options, "scale": 1.0}


def provenance(config: dict) -> dict:
    return {"tool_version": VERSION, "tool_sha256": sha256(Path(__file__).resolve()),
            "python_version": platform.python_version(), "numpy_version": np.__version__, "pyyaml_version": yaml.__version__,
            "config_file": str(config["_path"]), "config_sha256": sha256(config["_path"]),
            "config": {key: value for key, value in config.items() if not key.startswith("_")}}


def fit(config: dict, name: str, output: Path) -> dict:
    require(config["role"] == "calibration", "fit only accepts role: calibration")
    reference, ref_contract, estimates = load_inputs(config)
    selected = [e for e in estimates if e["name"] == name]
    require(len(selected) == 1, "--estimate must select a configured estimate")
    estimate = selected[0]
    pairs, unmatched = matched_rows(config, reference, estimate)
    rotation, translation, diagnostics = fit_rigid(pairs, config.get("fit", {}))
    result = {
        "schema_version": 1, "artifact_type": "frozen_reference_alignment", "status": "FROZEN_CALIBRATION_ONLY",
        "convention": "reference_from_estimate_world", "source_world_frame_id": estimate["metadata"]["world_frame_id"],
        "target_world_frame_id": ref_contract["metadata"]["world_frame_id"], "body_frame_id": BODY_FRAME,
        "transform": {"translation_m": translation.tolist(), "quaternion_xyzw": matrix_quat(rotation).tolist(), "scale": 1.0},
        "calibration": {"reference": ref_contract, "estimate": estimate["contract"],
                        "analysis_intervals": config["analysis_intervals"], "time": estimate["time"],
                        "matched_pose_count": len(pairs), "unmatched_pose_count": len(unmatched),
                        "diagnostics": diagnostics, "training_reference_disagreement": summarize([disagreement(a, b, rotation, translation) for a, b in pairs]),
                        "claim": "Calibration residuals are training diagnostics, not independent accuracy or acceptance results."},
        "conditions": list(dict.fromkeys(estimate["conditions"] + config.get("conditions", []))), "provenance": provenance(config),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return result


def load_alignment(path: Path, config: dict, reference: dict, estimate: dict) -> tuple[dict, np.ndarray, np.ndarray]:
    value = read_document(path)
    require(value.get("schema_version") == 1 and value.get("artifact_type") == "frozen_reference_alignment", "not a supported frozen alignment")
    require(value.get("status") == "FROZEN_CALIBRATION_ONLY" and value.get("convention") == "reference_from_estimate_world", "invalid alignment state/convention")
    require(value.get("source_world_frame_id") == estimate["metadata"]["world_frame_id"], "frozen transform belongs to another map/run frame")
    require(value.get("target_world_frame_id") == reference["metadata"]["world_frame_id"], "frozen transform belongs to another ENU origin/frame")
    require(value.get("body_frame_id") == BODY_FRAME, "frozen alignment uses another body/reference point")
    transform = value["transform"]
    require(transform.get("scale") == 1.0, "scale must remain exactly 1")
    rotation = quat_matrix(unit_quaternion(transform["quaternion_xyzw"], str(path)))
    translation = np.asarray(transform["translation_m"], dtype=float)
    require(translation.shape == (3,) and np.isfinite(translation).all(), "invalid alignment translation")
    calibration = value["calibration"]
    old_contracts = [calibration["reference"], calibration["estimate"]]
    for current in (reference, estimate["contract"]):
        for old in old_contracts:
            require(current["sha256"] != old["sha256"], "calibration/evaluation file hashes overlap; export independent data, never fit and test the same batch")
            require(current["canonical_pose_sha256"] != old.get("canonical_pose_sha256"),
                    "calibration/evaluation canonical pose hashes overlap despite file formatting changes")
            same_recording = current["metadata"]["source_recording_id"] == old["metadata"]["source_recording_id"]
            old_raw_hash = old["metadata"].get("source_recording_sha256")
            if old_raw_hash and old_raw_hash == current["metadata"].get("source_recording_sha256"):
                same_recording = True
            if same_recording:
                require(not intervals_overlap(config["analysis_intervals"], calibration["analysis_intervals"]),
                        "calibration/evaluation recording intervals overlap for the same source_recording_id")
    return value, rotation, translation


def acceptance(metrics: dict, exposure: dict, count: int, thresholds: Any, conditions: list[str]) -> dict:
    if not count:
        return {"status": "NO_PAIRED_DATA", "pass": None, "checks": {}, "reason": "No reference disagreement can be assessed; coverage still includes the complete declared duration."}
    if not thresholds:
        return {"status": "NOT_ASSESSED_NO_THRESHOLDS", "pass": None, "checks": {}, "reason": "No acceptance thresholds supplied."}
    supported = {
        "position_3d_rmse_m": (metrics["position_3d_m"]["rmse"], "max"),
        "position_3d_p95_m": (metrics["position_3d_m"]["p95"], "max"),
        "attitude_geodesic_rmse_deg": (metrics["attitude_geodesic_deg"]["rmse"], "max"),
        "attitude_geodesic_p95_deg": (metrics["attitude_geodesic_deg"]["p95"], "max"),
        "min_estimate_availability_fraction": (exposure["estimate_availability_fraction"], "min"),
        "min_reference_valid_fraction": (exposure["reference_valid_fraction"], "min"),
    }
    require(isinstance(thresholds, dict) and set(thresholds) <= set(supported), "unsupported acceptance threshold")
    checks = {}
    for key, raw in thresholds.items():
        limit = positive(raw, "threshold."+key)
        if key.startswith("min_"):
            require(limit <= 1, "fraction thresholds must be <= 1")
        observed, direction = supported[key]
        checks[key] = {"observed": observed, "limit": limit, "satisfied": observed <= limit if direction == "max" else observed >= limit}
    satisfied = all(check["satisfied"] for check in checks.values())
    status = "CONDITIONAL" if conditions else ("PASS" if satisfied else "FAIL")
    return {"status": status, "pass": None if conditions else satisfied, "checks": checks,
            "reason": "Applies only to the explicitly supplied thresholds and stated reference conditions."}


def evaluate(config: dict, output: Path) -> dict:
    require(config["role"] == "evaluation", "evaluate only accepts role: evaluation")
    reference, ref_contract, estimates = load_inputs(config)
    report = {"schema_version": 1, "report_type": "reference_disagreement", "provenance": provenance(config),
              "reference": ref_contract, "analysis_intervals": config["analysis_intervals"], "estimates": [],
              "interpretation": [
                  "CGI-430 is a quality-gated reference, not error-free ground truth; all errors below are reference disagreement.",
                  "Frozen SE(3) sets the map gauge using separate calibration data; results do not measure an uncalibrated map origin.",
                  "Pose metrics are estimate-sample-weighted and conditional on valid reference pairs; duration coverage is reported independently.",
                  "Signed errors are estimate minus reference. Body errors use reference rear-axle FLU. Euler differences are wrapped diagnostics; SO(3) geodesic is the primary attitude metric.",
                  "Availability uses the explicitly configured maximum sample age, without assuming a missing estimate is a valid zero error.",
                  "Availability is coverage of pose timestamps; offline trajectories do not establish live publication latency or wall-clock uptime.",
              ]}
    reference_reasons = {}
    selected_reference = [r for r in reference.rows if interval_for(config["analysis_intervals"], r.epoch, r.timestamp) is not None]
    for row in selected_reference:
        if not row.valid:
            reference_reasons[row.reason] = reference_reasons.get(row.reason, 0)+1
    report["reference_quality"] = {"input_rows": len(reference.rows), "rows_in_analysis_intervals": len(selected_reference),
                                   "valid_rows": sum(r.valid for r in selected_reference), "invalid_reasons": reference_reasons}
    alignments, sample_outputs, segment_outputs = {}, {}, {}
    ref_valid = reference.valid_intervals()
    for estimate in estimates:
        path = estimate["alignment_path"]
        require(path is not None, f"{estimate['name']}: evaluation requires frozen alignment")
        alignment, rotation, translation = load_alignment(path, config, ref_contract, estimate)
        frame, alignment_hash = estimate["metadata"]["world_frame_id"], sha256(path)
        require(frame not in alignments or alignments[frame] == alignment_hash,
                "trajectories in the same map/run frame must share the identical frozen alignment file")
        alignments[frame] = alignment_hash
        pairs, unmatched = matched_rows(config, reference, estimate)
        errors, samples = [], list(unmatched)
        for source, target in pairs:
            error = disagreement(source, target, rotation, translation)
            errors.append(error)
            samples.append({"timestamp_ns": target.timestamp, "estimate_timestamp_ns": source.timestamp,
                            "epoch": target.epoch, "status": "paired", **error})
        samples.sort(key=lambda x: (x["epoch"], x["timestamp_ns"]))
        available = estimate_availability(estimate["rows"], estimate["metadata"]["recorded_intervals"],
                                          estimate["hold_ns"],
                                          estimate["offset"], estimate["reference_epoch"])
        exposure, segments = coverage(config["analysis_intervals"], ref_valid, available)
        conditions = list(dict.fromkeys(estimate["conditions"] + alignment.get("conditions", []) + config.get("conditions", [])))
        reason_counts = {}
        for row in unmatched:
            reason_counts[row["status"]] = reason_counts.get(row["status"], 0) + 1
        metrics = summarize(errors)
        entry = {"name": estimate["name"], "trajectory_kind": estimate["metadata"]["trajectory_kind"],
                 "estimate": estimate["contract"], "frozen_alignment": {"file": str(path), "sha256": alignment_hash, "content": alignment},
                 "time": estimate["time"], "conditions": conditions, "coverage": exposure,
                 "coverage_criterion": {"maximum_estimate_age_sec": estimate["hold_ns"]/NS,
                                        "basis": "pose timestamps; sparse/offline sample coverage is not online publication uptime"},
                 "sample_counts": {"input_estimate": len(estimate["rows"]), "in_analysis_intervals": len(samples),
                                   "outside_analysis_intervals": len(estimate["rows"])-len(samples),
                                   "paired": len(pairs), "unpaired": len(unmatched), "unpaired_reasons": reason_counts,
                                   "euler_near_gimbal_lock": sum(e["euler_near_gimbal_lock"] for e in errors)},
                 "metrics": metrics, "acceptance": acceptance(metrics, exposure, len(pairs), config.get("thresholds"), conditions),
                 "segments": segments, "sample_csv": estimate["name"]+"_samples.csv", "coverage_csv": estimate["name"]+"_coverage.csv"}
        if entry["trajectory_kind"] == "optimized_dense_reconstructed":
            entry["interpretation"] = "Dense output reconstructed by interpolating optimized keyframe corrections; it is not independently optimized at every timestamp. Use optimized_keyframes for primary backend conclusions."
        report["estimates"].append(entry)
        sample_outputs[estimate["name"]], segment_outputs[estimate["name"]] = samples, segments
    report["common_sample_groups"], group_outputs = [], {}
    for group in config.get("common_sample_groups", []):
        members = group["estimates"]
        lookup = {name: {(row["epoch"], row["timestamp_ns"]): row for row in sample_outputs[name]} for name in members}
        all_keys = set.intersection(*(set(lookup[name]) for name in members))
        paired_keys = sorted(key for key in all_keys if all(lookup[name][key]["status"] == "paired" for name in members))
        metrics = {name: summarize([lookup[name][key] for key in paired_keys]) for name in members}
        report["common_sample_groups"].append({
            "name": group["name"], "estimates": members, "common_timestamp_count": len(all_keys),
            "paired_common_count": len(paired_keys), "metrics": metrics,
            "paired_rows_excluded_from_common": {name: sum(r["status"] == "paired" for r in lookup[name].values())-len(paired_keys) for name in members},
            "sample_csv": "group_"+group["name"]+"_samples.csv",
            "contract": "Exact intersection of reference epoch and timestamp; no pose resampling, time fitting, or re-alignment. Original per-stream coverage remains authoritative for outages.",
        })
        group_outputs[group["name"]] = [{"estimate_name": name, **lookup[name][key]} for key in paired_keys for name in members]
    # Validation and computation happen before creating the output directory.
    output.mkdir(parents=True, exist_ok=False)
    for name, samples in sample_outputs.items():
        fields = ["timestamp_ns", "estimate_timestamp_ns", "epoch", "status"] + METRICS + ["euler_near_gimbal_lock", "estimate_enu_x", "estimate_enu_y", "estimate_enu_z", "reference_x", "reference_y", "reference_z"]
        write_csv(output / (name+"_samples.csv"), fields, samples)
        write_csv(output / (name+"_coverage.csv"), ["epoch", "start_ns", "end_ns", "state", "duration_sec"], segment_outputs[name])
    for name, rows in group_outputs.items():
        write_csv(output / ("group_"+name+"_samples.csv"), ["estimate_name"]+fields, rows)
    with (output/"report.json").open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    with (output/"report.md").open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(markdown_report(report))
    return report


def write_csv(path: Path, names: list[str], rows: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def formatted(value: Any) -> str:
    return "N/A" if value is None else format(value, ".6g")


def markdown_report(report: dict) -> str:
    lines = ["# 轨迹与 CGI-430 参考差异报告", "", "本报告中的误差均为 **reference disagreement（与参考的差异）**。CGI-430 不被视为无误差真值。", "",
             f"工具版本：`{report['provenance']['tool_version']}`；配置 SHA-256：`{report['provenance']['config_sha256']}`。",
             "", "完整配置、源文件及元数据哈希、坐标系、时钟、冻结标定和连续区段均见 `report.json`。", "",
             "地图坐标经独立数据冻结的 SE(3) 变换定义到 ENU；不报告未经标定的地图原点绝对精度。位置按估计减参考计算，FLU 误差在参考后轴车体系表达。", "",
             "位姿统计按实际估计样本加权；有效参考和估计可用率的分母分别保留完整声明录制时长。未配对样本与停发时段不作为零误差。", ""]
    for result in report["estimates"]:
        c = result["coverage"]
        lines += [f"## {result['name']} ({result['trajectory_kind']})", "", f"判定：**{result['acceptance']['status']}**。{result['acceptance']['reason']}", "",
                  f"录制时长 {formatted(c['recorded_duration_sec'])} s；参考有效 {formatted(100*c['reference_valid_fraction'])}%；估计可用 {formatted(100*c['estimate_availability_fraction'])}%；共同覆盖 {formatted(100*c['joint_coverage_fraction'])}%。",
                  f"配对样本 {result['sample_counts']['paired']}；未配对样本 {result['sample_counts']['unpaired']}；估计不可用 {formatted(c['estimate_unavailable_duration_sec'])} s。", "",
                  f"冻结变换 SHA-256：`{result['frozen_alignment']['sha256']}`。", ""]
        if result.get("interpretation"):
            lines += [result["interpretation"], ""]
        if result["conditions"]:
            lines += ["条件证据：" + "；".join(result["conditions"]), ""]
        lines += ["| 指标 | RMSE | mean | median | P95 | max | bias | P95 abs | max abs |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for key, metric in result["metrics"].items():
            lines += ["| " + key + " | " + " | ".join(formatted(metric[k]) for k in ("rmse", "mean", "median", "p95", "max", "bias", "p95_abs", "max_abs")) + " |"]
        lines += ["", "SO(3) geodesic 为姿态主指标。roll/pitch/yaw 是逐姿态欧拉角之差并 wrap 到 [-180°,180°)，接近万向节锁时仅作诊断。", "",
                  f"逐样本审计：`{result['sample_csv']}`；完整时段覆盖与停发：`{result['coverage_csv']}`。", ""]
    for group in report["common_sample_groups"]:
        lines += [f"## 同样本对比：{group['name']}", "",
                  f"同参考 epoch 与时间戳精确交集 {group['common_timestamp_count']} 个，其中所有轨迹均有有效参考配对 {group['paired_common_count']} 个。以下各路使用完全相同样本，缺失覆盖仍见各路原报告。", "",
                  "| 轨迹 | 3D RMSE (m) | 3D P95 (m) | 姿态 geodesic RMSE (deg) | 姿态 geodesic P95 (deg) |", "|---|---:|---:|---:|---:|"]
        for name, metrics in group["metrics"].items():
            values = [metrics["position_3d_m"]["rmse"], metrics["position_3d_m"]["p95"],
                      metrics["attitude_geodesic_deg"]["rmse"], metrics["attitude_geodesic_deg"]["p95"]]
            lines += ["| "+name+" | "+" | ".join(formatted(v) for v in values)+" |"]
        lines += ["", f"完整同样本指标见 JSON；逐样本交集见 `{group['sample_csv']}`。", ""]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ("validate", "fit", "evaluate"):
        child = commands.add_parser(command)
        child.add_argument("--config", required=True, type=Path)
        if command != "validate":
            child.add_argument("--output", required=True, type=Path)
        if command == "fit":
            child.add_argument("--estimate", required=True)
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config)
        if args.command == "validate":
            reference, _, estimates = load_inputs(config)
            result = {"valid": True, "reference_rows": len(reference.rows), "estimates": {e["name"]: len(e["rows"]) for e in estimates}}
        elif args.command == "fit":
            value = fit(config, args.estimate, args.output)
            result = {"output": str(args.output.resolve()), "status": value["status"], "matched_pose_count": value["calibration"]["matched_pose_count"]}
        else:
            value = evaluate(config, args.output)
            result = {"output": str(args.output.resolve()), "estimates": {e["name"]: e["acceptance"]["status"] for e in value["estimates"]}}
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except (ValueError, KeyError, TypeError, OSError, yaml.YAMLError) as error:
        print(json.dumps({"error": str(error)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
