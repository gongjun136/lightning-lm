#!/usr/bin/env python3
"""Offline reference-point translation and frozen-alignment lever sensitivity.

Only a rigid OUTPUT-point translation is supported, never the GNSS/IMU lever
inside a navigation filter. Quaternion axes already describe rear-body FLU in
the fixed reference world. Invalid rows remain source-point audit evidence.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import sys

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent / "accuracy_eval"))
import accuracy_eval as ae

VERSION = "1.0.0"
SENSITIVITY = "LEVER_ARM_SENSITIVITY_NOT_ACCURACY_ACCEPTANCE"
COVARIANCE = "REFERENCE_LEVER_ARM_COVARIANCE_NOT_PROPAGATED"
ATTITUDE_METRICS = ("roll_deg", "pitch_deg", "yaw_deg", "attitude_geodesic_deg")


def nonempty(value, label):
    ae.require(isinstance(value, str) and bool(value.strip()), f"{label}: explicit nonempty evidence/identity required")
    return value


def vector(value, label):
    ae.require(isinstance(value, list) and len(value) == 3
               and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in value),
               f"{label}: require three numeric metres in rear_axle_flu")
    result = np.asarray(value, dtype=float)
    ae.require(np.isfinite(result).all(), f"{label}: NaN/infinite lever arm is forbidden")
    return result


def snapshot(paths):
    return {str(Path(p).resolve()): ae.sha256(Path(p)) for p in paths}


def unchanged(hashes):
    for name, expected in hashes.items():
        ae.require(ae.sha256(Path(name)) == expected, f"input/frozen alignment changed during processing: {name}")


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def normalization_config(path):
    config = ae.read_document(path)
    ae.require(type(config.get("schema_version")) is int and config["schema_version"] == 1,
               "unsupported reference normalization schema")
    ae.require(config.get("mode") in {"fixed_point", "rear_axle_delta", "already_rear"}, "unsupported normalization mode")
    ae.require(config.get("confirmed") is True, "reference normalization must be confirmed")
    nonempty(config.get("evidence"), "evidence")
    nonempty(config.get("source_point_id"), "source_point_id")
    ae.require(config.get("rigid_to_rear_axle_confirmed") is True,
               "rigid_to_rear_axle_confirmed is required; a constant lever cannot span articulation")
    lever = vector(config.get("lever_arm_m"), "lever_arm_m")
    delta = lever.copy()
    if config["mode"] == "rear_axle_delta":
        old = vector(config.get("previous_device_lever_arm_m"), "previous_device_lever_arm_m")
        ae.require(config.get("previous_compensation_is_rigid_translation_confirmed") is True,
                   "previous compensation must be confirmed as a rigid output translation, not internal filter calibration")
        nonempty(config.get("previous_compensation_evidence"), "previous_compensation_evidence")
        delta -= old
    elif config["mode"] == "already_rear":
        ae.require(np.array_equal(lever, np.zeros(3)), "already_rear only permits zero lever_arm_m")
    conditions = config.get("conditions", [])
    ae.require(isinstance(conditions, list) and all(isinstance(c, str) for c in conditions), "conditions must be a string list")
    return config, delta


def read_source(source, metadata, config):
    """Validate raw metadata independently; reuse the canonical row validator."""
    meta = ae.read_document(metadata)
    ae.require(type(meta.get("schema_version")) is int and meta["schema_version"] == 1, "unsupported source metadata schema")
    for field in ("source_recording_id", "world_frame_id", "body_frame_id", "time_basis"):
        nonempty(meta.get(field), f"metadata.{field}")
    ae.require("reference_normalization" not in meta, "reference was already normalized offline; use the original source to avoid repeated compensation")
    conditions = meta.get("conditions", [])
    ae.require(isinstance(conditions, list) and all(isinstance(c, str) for c in conditions), "metadata.conditions must be a string list")
    intervals = meta.get("recorded_intervals")
    ae.require(isinstance(intervals, list) and bool(intervals), "recorded_intervals required")
    for interval in intervals:
        ae.require(isinstance(interval, dict), "invalid recorded interval")
        nonempty(interval.get("epoch"), "recorded_intervals.epoch")
    parsed_intervals = ae.parse_intervals(intervals, "recorded_intervals")
    if "clock_epoch" in meta:
        nonempty(meta["clock_epoch"], "clock_epoch")
        ae.require({i["epoch"] for i in parsed_intervals} == {meta["clock_epoch"]}, "clock_epoch differs from recorded intervals")
    point = meta.get("reference_point")
    if config["mode"] == "fixed_point":
        ae.require(meta["body_frame_id"] == "cgi_output_flu", "fixed_point requires cgi_output_flu source, not an already compensated rear axle")
        ae.require(isinstance(point, dict) and point.get("position") == "fixed_point", "source reference_point must declare fixed_point")
        ae.require(point.get("point_id") == config["source_point_id"], "source_point_id does not match source reference_point.point_id")
        ae.require(point.get("attitude_axes") == "FLU", "source attitude axes must already be rear-body FLU")
        ae.require(point.get("device_lever_arm_already_applied") is False,
                   "fixed_point requires O-to-rear OUTPUT translation not already applied; retain the internal GNSS-IMU calibration")
        ae.require(point.get("confirmed") is True and point.get("rigid_to_rear_axle_confirmed") is True,
                   "source fixed point/rigid-body contract must be confirmed")
        nonempty(point.get("evidence"), "reference_point.evidence")
    else:
        ae.require(meta["body_frame_id"] == ae.BODY_FRAME, "rear axle modes require rear_axle_flu source")
        # Historical exporter products have no reference_point object. Explicit
        # newer declarations must agree; never relabel a fixed-point product.
        if point is not None:
            ae.require(isinstance(point, dict) and point.get("position") in {"rear_axle_center", "rear_axle"}
                       and point.get("attitude_axes") == "FLU" and point.get("confirmed") is True,
                       "source reference_point contradicts rear-axle mode")
            ae.require(point.get("rigid_to_rear_axle_confirmed", True) is True,
                       "source denies a rigid rear-body relationship")
    with source.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        ae.require(len(fields) == len(set(fields)) and ae.REFERENCE_COLUMNS <= set(fields), "invalid/duplicate canonical CSV columns")
        rows = list(reader)
    ae.require(all(None not in row and all(v is not None for v in row.values()) for row in rows), "CSV row width differs from header")
    poses = ae.read_reference(source, {**meta, "recorded_intervals": parsed_intervals})
    ae.require(len(rows) == len(poses), "CSV row count changed while reading")
    return meta, fields, rows, poses


def normalize(source, metadata, config_path, output):
    source, metadata, config_path, output = [Path(p).resolve() for p in (source, metadata, config_path, output)]
    output_meta = output.with_suffix(".metadata.json")
    ae.require(output.suffix.lower() == ".csv", "normalization output must be a CSV file")
    ae.require(not output.exists() and not output_meta.exists(), "output CSV/metadata already exists; refusing overwrite")
    ae.require(output not in (source, metadata, config_path) and output_meta not in (source, metadata, config_path), "output must not replace an input")
    hashes = snapshot([source, metadata, config_path])
    config, delta = normalization_config(config_path)
    meta, fields, rows, poses = read_source(source, metadata, config)
    for row, pose in zip(rows, poses):
        if pose.valid and np.any(delta != 0):
            position = pose.position + ae.quat_matrix(pose.quaternion) @ delta
            ae.require(np.isfinite(position).all(), "translated position is not finite")
            row.update(zip(("x", "y", "z"), (format(float(v), ".17g") for v in position)))
    text = io.StringIO(newline="")
    writer = csv.DictWriter(text, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    csv_text = text.getvalue()
    result = copy.deepcopy(meta)
    result["body_frame_id"] = ae.BODY_FRAME
    result["conditions"] = list(meta.get("conditions", []))
    for condition in config.get("conditions", []) + ([COVARIANCE] if np.any(delta != 0) else []):
        if condition not in result["conditions"]:
            result["conditions"].append(condition)
    result["reference_point"] = {"position": "rear_axle_center", "point_id": "rear_axle_center", "attitude_axes": "FLU",
                                 "confirmed": True, "rigid_to_rear_axle_confirmed": True,
                                 "offline_translation_applied": True, "evidence": config["evidence"]}
    result["reference_normalization"] = {
        "tool_version": VERSION, "tool_sha256": ae.sha256(Path(__file__).resolve()),
        "mode": config["mode"], "config_file": str(config_path), "config_sha256": hashes[str(config_path)], "config": config,
        "source": ae.file_contract(source, metadata, meta, poses),
        "lever_arm_frame": ae.BODY_FRAME, "lever_arm_direction": "original fixed output point to physical rear axle",
        "applied_delta_m": delta.tolist(), "formula": "p_target = p_source + R_world_rear_body * applied_delta_m",
        "quaternion_time_validity_and_quality_columns": "preserved without change",
        "invalid_rows": "All fields preserved as source-point audit evidence; invalid XYZ are NOT rear-axle positions and must not be evaluated.",
        "quality_fields_reference_point": copy.deepcopy(meta.get("reference_point", {"position": "rear_axle_center", "declaration": "legacy body_frame_id"})),
        "target_covariance_status": "not_propagated" if np.any(delta != 0) else "unchanged_no_position_shift",
        "quality_gate": "source-point quality gate preserved; no target covariance or new validity mask is inferred",
        "output_csv_sha256": hashlib.sha256(csv_text.encode("utf-8")).hexdigest(),
        "rows": len(rows), "valid_rows": sum(p.valid for p in poses),
    }
    json_text = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    unchanged(hashes)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Both paths are exclusive. Do not modify or remove any existing artifact.
    with output.open("x", encoding="utf-8", newline="") as csv_stream, output_meta.open("x", encoding="utf-8", newline="\n") as meta_stream:
        csv_stream.write(csv_text)
        meta_stream.write(json_text)
    return result


def sample_signature(report, directory):
    signature = {}
    for estimate in report["estimates"]:
        with (directory / estimate["sample_csv"]).open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        signature[estimate["name"]] = {
            "samples": [[row.get(k, "") for k in ("epoch", "timestamp_ns", "estimate_timestamp_ns", "status", *ATTITUDE_METRICS)] for row in rows],
            "attitude_metrics": {key: estimate["metrics"][key] for key in ATTITUDE_METRICS},
            "sample_counts": estimate["sample_counts"], "coverage": estimate["coverage"], "segments": estimate["segments"],
        }
    return signature


def sweep(source, metadata, config_path, output):
    source, metadata, config_path, output = [Path(p).resolve() for p in (source, metadata, config_path, output)]
    ae.require(not output.exists(), "sweep output directory already exists")
    config = ae.read_document(config_path)
    ae.require(type(config.get("schema_version")) is int and config["schema_version"] == 1, "unsupported sweep schema")
    normal_path = ae.resolve(config_path.parent, config.get("normalization_config"))
    evaluation_path = ae.resolve(config_path.parent, config.get("evaluation_config"))
    normal, _ = normalization_config(normal_path)
    ae.require(normal["mode"] != "already_rear", "already_rear is not a lever sweep mode")
    read_source(source, metadata, normal)
    evaluation = ae.load_config(evaluation_path)
    ae.require(evaluation["role"] == "evaluation", "sweep only accepts role: evaluation; it never fits an alignment")
    candidates = config.get("candidates")
    ae.require(isinstance(candidates, list) and bool(candidates), "at least one explicit candidate is required")
    names = set()
    for candidate in candidates:
        ae.require(isinstance(candidate, dict), "each candidate must be an object")
        name = candidate.get("name")
        ae.require(isinstance(name, str) and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", name) is not None,
                   "candidate name must be 1-64 ASCII letters/digits/underscore/hyphen, starting with a letter or digit")
        ae.require(name.casefold() not in names and re.fullmatch(r"(?i:con|prn|aux|nul|com[1-9]|lpt[1-9])", name) is None,
                   "candidate names must be unique, case-insensitive, and not reserved device names")
        names.add(name.casefold())
        vector(candidate.get("lever_arm_m"), f"{name}.lever_arm_m")
    # Resolve against the ORIGINAL evaluation YAML, never a candidate folder.
    base = {key: copy.deepcopy(value) for key, value in evaluation.items() if not key.startswith("_")}
    inputs = [source, metadata, config_path, normal_path, evaluation_path]
    alignment_paths = []
    for estimate in base["estimates"]:
        for key in ("tum", "metadata", "alignment"):
            path = ae.resolve(evaluation_path.parent, estimate.get(key))
            estimate[key] = str(path)
            inputs.append(path)
            if key == "alignment":
                alignment_paths.append(path)
    hashes = snapshot(inputs)
    ae.require(len({hashes[str(p)] for p in alignment_paths}) == 1,
               "one lever sweep must use the same frozen alignment for all estimates")
    base["thresholds"] = None
    base["conditions"] = list(dict.fromkeys(base.get("conditions", []) + [SENSITIVITY]))
    output.mkdir(parents=True, exist_ok=False)
    result = {"schema_version": 1, "report_type": "lever_arm_sensitivity", "status": "SENSITIVITY_ONLY",
              "automatic_selection": False, "thresholds": None, "conditions": [SENSITIVITY],
              "interpretation": "仅比较输出点杆臂敏感性；不选择最优杆臂、不拟合地图或时间、不作精度验收。用本数据选择杆臂后须另取独立验证数据。",
              "calibration_caveat": "冻结地图配准可能已吸收原杆臂误差；扫描最小残差不能证明物理杆臂正确。内部GNSS-IMU杆臂不在本工具范围内。",
              "provenance": {"tool_version": VERSION, "tool_sha256": ae.sha256(Path(__file__).resolve()),
                             "evaluation_tool_sha256": ae.sha256(Path(ae.__file__)), "input_sha256": hashes,
                             "sweep_config": config, "base_thresholds_ignored": evaluation.get("thresholds")},
              "frozen_alignments": {str(p): hashes[str(p)] for p in alignment_paths}, "candidates": []}
    reference_signature = None
    summary_rows = []
    for candidate in candidates:
        unchanged(hashes)
        directory = output / candidate["name"]
        directory.mkdir()
        candidate_normal = {**normal, "lever_arm_m": candidate["lever_arm_m"]}
        candidate_normal["conditions"] = list(dict.fromkeys(normal.get("conditions", []) + [SENSITIVITY]))
        write_json(directory / "normalization.yaml", candidate_normal)
        normalize(source, metadata, directory / "normalization.yaml", directory / "reference_rear.csv")
        candidate_evaluation = copy.deepcopy(base)
        candidate_evaluation["reference"] = {"csv": str(directory / "reference_rear.csv"), "metadata": str(directory / "reference_rear.metadata.json")}
        write_json(directory / "evaluation.yaml", candidate_evaluation)
        report_directory = directory / "report"
        report = ae.evaluate(ae.load_config(directory / "evaluation.yaml"), report_directory)
        unchanged(hashes)
        signature = sample_signature(report, report_directory)
        if reference_signature is None:
            reference_signature = signature
        ae.require(signature == reference_signature, "sweep candidates changed pairing/status, attitude metrics or duration coverage")
        for estimate in report["estimates"]:
            ae.require(estimate["acceptance"]["pass"] is None and estimate["acceptance"]["status"] in {"NOT_ASSESSED_NO_THRESHOLDS", "NO_PAIRED_DATA"},
                       "sensitivity report must never make an accuracy acceptance decision")
            ae.require(estimate["frozen_alignment"]["sha256"] == hashes[estimate["frozen_alignment"]["file"]], "candidate alignment differs from frozen input")
            for condition in estimate["conditions"]:
                if condition not in result["conditions"]:
                    result["conditions"].append(condition)
            for metric, values in estimate["metrics"].items():
                summary_rows.append({"candidate": candidate["name"], "lever_x_m": candidate["lever_arm_m"][0],
                                     "lever_y_m": candidate["lever_arm_m"][1], "lever_z_m": candidate["lever_arm_m"][2],
                                     "estimate": estimate["name"], "metric": metric, **values,
                                     "paired": estimate["sample_counts"]["paired"], "status": estimate["acceptance"]["status"],
                                     "recorded_duration_sec": estimate["coverage"]["recorded_duration_sec"],
                                     "estimate_availability_fraction": estimate["coverage"]["estimate_availability_fraction"],
                                     "reference_valid_fraction": estimate["coverage"]["reference_valid_fraction"]})
        result["candidates"].append({"name": candidate["name"], "lever_arm_m": candidate["lever_arm_m"],
                                     "report": str(report_directory / "report.json"), "report_sha256": ae.sha256(report_directory / "report.json"),
                                     "estimates": [{k: e[k] for k in ("name", "trajectory_kind", "metrics", "coverage", "sample_counts", "acceptance", "conditions")} for e in report["estimates"]]})
    if not any(e["sample_counts"]["paired"] for c in result["candidates"] for e in c["estimates"]):
        result["status"] = "NO_PAIRED_DATA"
    result["invariants"] = {"input_and_alignment_hashes_unchanged": True, "same_timestamp_status_population": True,
                            "same_attitude_metrics": True, "same_duration_coverage": True}
    unchanged(hashes)
    write_json(output / "sweep.json", result)
    ae.write_csv(output / "summary.csv", list(summary_rows[0]), summary_rows)
    lines = ["# 杆臂敏感性分析", "", result["interpretation"], "", result["calibration_caveat"], "",
             "以下按配置顺序列出，未进行最优候选选择。零配对仍保留完整录制时段和缺测分母。", "",
             "| 候选 | 轨迹 | 杆臂 FLU（m） | 配对数 | 3D RMSE（m） | SO(3) RMSE（deg） | 状态 |", "|---|---|---|---:|---:|---:|---|"]
    for candidate in result["candidates"]:
        for estimate in candidate["estimates"]:
            def number(value):
                return "无配对数据" if value is None else f"{value:.6g}"
            lines.append(f"| {candidate['name']} | {estimate['name']} | {candidate['lever_arm_m']} | {estimate['sample_counts']['paired']} | "
                         f"{number(estimate['metrics']['position_3d_m']['rmse'])} | {number(estimate['metrics']['attitude_geodesic_deg']['rmse'])} | {estimate['acceptance']['status']} |")
    lines += ["", "完整指标见 summary.csv/sweep.json，各候选 report/ 保留逐样本残差与覆盖分段。", "",
              "已核验：配准与所有输入哈希未变；候选的配对时间/状态、姿态指标及覆盖分母完全相同。", "", "条件与证据：", ""]
    lines += ["- " + condition for condition in result["conditions"]]
    with (output / "sweep.md").open("x", encoding="utf-8", newline="\n") as stream:
        stream.write("\n".join(lines) + "\n")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("normalize", "sweep"):
        child = commands.add_parser(name)
        for flag in ("source", "metadata", "config", "output"):
            child.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        result = (normalize if args.command == "normalize" else sweep)(args.source, args.metadata, args.config, args.output)
        print(json.dumps({"output": str(args.output.resolve()), "status": result.get("status", "NORMALIZED"),
                          "claim": "Reference-point processing only; no accuracy acceptance or automatic lever selection."}, ensure_ascii=False))
        return 0
    except (ValueError, TypeError, KeyError, OSError, yaml.YAMLError) as error:
        print("accuracy_reference: " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
