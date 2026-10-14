#!/usr/bin/env python3
"""Generate deterministic SYNTHETIC A/B/C data; optionally execute the real CLI.

This exercises the offline assessment and report contracts only. It is neither
a CGI-430 device simulation nor a field accuracy/ROS capture validation.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import accuracy_eval as ae

MARKER = "SYNTHETIC_DATA_NOT_A_FIELD_RESULT"
BASE_NS = 1_800_000_000_000_000_000
OFFSET_NS = 35_000_000
DT_NS = 200_000_000
COUNT = 100
EPOCH = "synthetic_epoch_0"
SEED = 430


def rotation(roll, pitch, yaw):
    a, b, c = np.radians([roll, pitch, yaw])
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    return rz @ ry @ rx


def document(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")


def tum_time(ns):
    return f"{ns//ae.NS}.{ns%ae.NS:09d}"


def generate(output):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(SEED)
    # Seed determines a small fixed rotation component; no measurement noise is
    # added, so the injected pose disagreements have analytically known norms.
    map_rotation = rotation(10, -7, 42+float(rng.uniform(-1, 1)))
    map_translation = np.array([100., -75., 4.5])
    manifest = {"marker": MARKER, "seed": SEED, "asset_kind": "synthetic_map_identity_only", "scale": 1.0,
                "expected_reference_from_map": {"rotation_matrix": map_rotation.tolist(), "translation_m": map_translation.tolist()}}
    document(output/"synthetic_map_manifest.json", manifest)
    map_id = "map:"+ae.sha256(output/"synthetic_map_manifest.json")
    reference_frame = "enu:SYNTHETIC_ORIGIN_31N_121E_8M"
    (output/"SYNTHETIC_NOTICE.md").write_text(
        "# SYNTHETIC DEMONSTRATION — 不是实车结果\n\n"
        "全部位姿、地图标识、时钟和质量标志由固定 seed 430 的数学轨迹生成。"
        "本目录只验证离线代码和报告，不证明 CGI-430、CAN、rosbag2 或车上算法的实测精度。\n",
        encoding="utf-8")
    expected = {"marker": MARKER, "seed": SEED, "position_error_sign": "estimate minus reference in ENU", "streams": {}}
    settings = {
        "A": (0, .0, set(range(35, 40))),
        "B": (200, .3, set()),
        "C": (400, .7, set(range(70, 75))),
    }
    configs = {}
    for run, (start_seconds, phase_shift, invalid) in settings.items():
        directory = output/run
        directory.mkdir()
        start = BASE_NS+start_seconds*ae.NS
        end = start+COUNT*DT_NS
        interval = {"start_ns": start, "end_ns": end, "epoch": EPOCH}
        source_manifest = {"marker": MARKER, "source_recording_id": "synthetic_"+run, "seed": SEED,
                           "start_ns": start, "end_ns": end, "count": COUNT, "dt_ns": DT_NS}
        document(directory/"source_manifest.json", source_manifest)
        metadata = {"schema_version": 1, "source_recording_id": "synthetic_"+run,
                    "source_recording_sha256": ae.sha256(directory/"source_manifest.json"),
                    "world_frame_id": reference_frame, "body_frame_id": ae.BODY_FRAME,
                    "time_basis": "synthetic_unix_nanoseconds", "clock_epoch": EPOCH,
                    "recorded_intervals": [interval], "conditions": [MARKER],
                    "evidence": {"marker": MARKER, "pose_convention": "world_from_body", "origin_llh": [31.0, 121.0, 8.0],
                                 "position_unit": "m", "quaternion_order": "xyzw", "seed": SEED}}
        document(directory/"reference.metadata.json", metadata)
        truth, reference_rows = [], []
        for index in range(COUNT):
            theta = phase_shift+index*.06
            point = np.array([20*np.cos(theta), 12*np.sin(theta), .4*np.sin(2*theta)])
            attitude = rotation(2*np.sin(theta), 1.5*np.cos(theta), math.degrees(theta)+90)
            timestamp = start+index*DT_NS
            ref_position = map_rotation @ point+map_translation
            ref_quaternion = ae.matrix_quat(map_rotation @ attitude)
            valid = index not in invalid
            reference_rows.append([timestamp, *ref_position, *ref_quaternion, int(valid), "" if valid else "synthetic_RTK_FLOAT", EPOCH, metadata["time_basis"]])
            truth.append((timestamp, point, attitude))
        with (directory/"reference.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["timestamp_ns", "x", "y", "z", "qx", "qy", "qz", "qw", "reference_valid", "reason", "epoch", "time_basis"])
            writer.writerows(reference_rows)
        if run == "A":
            streams = [
                ("frontend_dense", "frontend", [0.30, -0.40, .10], 2., list(range(COUNT)), .25),
                ("frontend_keyframes", "frontend", [0.30, -0.40, .10], 2., list(range(0, COUNT, 5)), 1.1),
                ("optimized_keyframes", "optimized_keyframes", [.06, -.08, .02], .4, list(range(0, COUNT, 5)), 1.1),
                ("optimized_dense", "optimized_dense_reconstructed", [.06, -.08, .02], .4, list(range(COUNT)), .25),
            ]
        elif run == "B":
            streams = [("calibration_localization", "localization", [0., 0., 0.], 0., list(range(COUNT)), .25)]
        else:
            kept = [i for i in range(COUNT) if i not in set(range(5)) | set(range(40, 50)) | set(range(95, 100))]
            streams = [("lidar_localization", "localization", [-.12, .16, -.05], 1., kept, .25)]
        estimates = []
        for name, kind, bias, angle, kept, hold in streams:
            source_bias = map_rotation.T @ np.array(bias)
            with (directory/(name+".tum")).open("w", encoding="utf-8", newline="\n") as stream:
                stream.write("# "+MARKER+"; canonical rear axle FLU; world_from_body\n")
                for index in kept:
                    timestamp, point, attitude = truth[index]
                    q = ae.matrix_quat(attitude @ rotation(0, 0, angle))
                    stream.write(" ".join([tum_time(timestamp-OFFSET_NS), *[format(v, ".16g") for v in np.r_[point+source_bias, q]]])+"\n")
            est_meta = {**metadata, "world_frame_id": map_id, "trajectory_kind": kind,
                        "recorded_intervals": ae.shifted([interval], -OFFSET_NS),
                        "evidence": {**metadata["evidence"], "synthetic_injected_enu_bias_m": bias, "synthetic_injected_geodesic_deg": angle,
                                     "normalization": "Already canonical rear axle; no lever arm is applied by accuracy_eval."}}
            document(directory/(name+".metadata.json"), est_meta)
            item = {"name": name, "tum": run+"/"+name+".tum", "metadata": run+"/"+name+".metadata.json", "availability_hold_sec": hold}
            if run != "B":
                item["alignment"] = "frozen_map_to_enu.json"
                expected["streams"][name] = {"position_3d_rmse_m": float(np.linalg.norm(bias)), "attitude_geodesic_rmse_deg": angle}
            estimates.append(item)
        config = {"schema_version": 1, "role": "calibration" if run == "B" else "evaluation",
                  "reference": {"csv": run+"/reference.csv", "metadata": run+"/reference.metadata.json"}, "estimates": estimates,
                  "analysis_intervals": [interval], "conditions": [MARKER], "thresholds": None,
                  "time": {"offset_status": "known", "estimate_to_reference_offset_ns": OFFSET_NS, "reference_epoch": EPOCH,
                           "evidence": "SYNTHETIC fixed 35 ms clock offset; t_ref = t_est + 35000000 ns.", "allow_conditional": False},
                  "pairing": {"max_reference_gap_sec": .25, "availability_hold_sec": .25}}
        if run == "B":
            config["fit"] = {"method": "position_attitude", "min_pairs": 30, "min_secondary_span_m": .5,
                             "attitude_weight_m": 2., "max_attitude_dispersion_deg": 1.}
        if run == "A":
            config["common_sample_groups"] = [{"name": "mapping_before_after", "estimates": ["frontend_keyframes", "optimized_keyframes"]}]
        config_path = output/("config_"+run.lower()+".yaml")
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        configs[run] = config_path
    expected["localization_recorded_duration_sec"] = 20.0
    expected["localization_estimate_availability_fraction"] = .805
    expected["mapping_before_after_common_pose_count"] = 19
    expected["tolerance"] = 1e-8
    document(output/"expected.json", expected)
    script = Path(ae.__file__).resolve()
    commands = [[sys.executable, str(script), "validate", "--config", str(configs["B"])],
                [sys.executable, str(script), "fit", "--config", str(configs["B"]), "--estimate", "calibration_localization", "--output", str(output/"frozen_map_to_enu.json")]]
    for run in ("A", "C"):
        commands += [[sys.executable, str(script), "validate", "--config", str(configs[run])],
                     [sys.executable, str(script), "evaluate", "--config", str(configs[run]), "--output", str(output/("report_"+run.lower()))]]
    document(output/"run_commands.json", {"marker": MARKER, "commands": commands})
    return commands


def verify(output):
    expected = json.loads((output/"expected.json").read_text())
    for report_name in ("report_a", "report_c"):
        report = json.loads((output/report_name/"report.json").read_text(encoding="utf-8"))
        for entry in report["estimates"]:
            target = expected["streams"][entry["name"]]
            for name, metric in (("position_3d_rmse_m", "position_3d_m"), ("attitude_geodesic_rmse_deg", "attitude_geodesic_deg")):
                if abs(entry["metrics"][metric]["rmse"]-target[name]) > expected["tolerance"]:
                    raise AssertionError(f"synthetic known error mismatch for {entry['name']} {name}")
            if entry["acceptance"]["pass"] is not None:
                raise AssertionError("synthetic/no-threshold report must not claim PASS")
        if report_name == "report_a":
            assert report["common_sample_groups"][0]["paired_common_count"] == expected["mapping_before_after_common_pose_count"]
        else:
            c = report["estimates"][0]["coverage"]
            assert c["recorded_duration_sec"] == expected["localization_recorded_duration_sec"]
            assert abs(c["estimate_availability_fraction"]-expected["localization_estimate_availability_fraction"]) < expected["tolerance"]
    document(output/"synthetic_verification.json", {"marker": MARKER, "known_error_checks": "passed", "tolerance": expected["tolerance"]})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="new directory only")
    parser.add_argument("--run-core", action="store_true", help="run validate/fit/evaluate CLIs and assert the known synthetic errors")
    args = parser.parse_args()
    commands = generate(args.output)
    if args.run_core:
        for command in commands:
            subprocess.run(command, check=True)
        verify(args.output)
    print(json.dumps({"marker": MARKER, "output": str(args.output.resolve()), "core_executed": args.run_core}, indent=2))


if __name__ == "__main__":
    main()
