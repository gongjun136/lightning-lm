#!/usr/bin/env python3
"""Mapping-A-only aligned diagnostics; never independent accuracy or acceptance.

Fits one position-only SE(3) to a declared common-keyframe set in A and applies
it unchanged to every A trajectory. Does not produce a formal frozen alignment.
Unconfirmed physical contracts must be declared in metadata/config conditions.
"""
from __future__ import annotations

import argparse
import bisect
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

import accuracy_eval as ae

VERSION = "1.1.0"
STATUS = "A_ONLY_ALIGNED_DIAGNOSTIC"
LIMITATION = (
    "A_ONLY_ALIGNED_DIAGNOSTIC: alignment is fitted on the same mapping A data; "
    "residuals are in-sample reference disagreement, not independently calibrated "
    "absolute accuracy, localization accuracy, or an acceptance result."
)
SAMPLE_FIELDS = ["timestamp_ns", "estimate_timestamp_ns", "epoch", "status"] + ae.METRICS + [
    "euler_near_gimbal_lock", "estimate_enu_x", "estimate_enu_y", "estimate_enu_z",
    "reference_x", "reference_y", "reference_z",
]
RPE_METRICS = ["translation_3d_m", "error_x_m", "error_y_m", "error_z_m", "rotation_geodesic_deg"]
RPE_FIELDS = ["timestamp_ns", "end_timestamp_ns", "estimate_timestamp_ns", "epoch", "status", "detail"] + RPE_METRICS


def load_config(path: Path) -> dict:
    config = ae.read_document(path)
    ae.require(config.get("schema_version") == 1, "unsupported diagnostic schema")
    ae.require(config.get("role") == "diagnostic_mapping_a", "require role: diagnostic_mapping_a")
    ae.require(config.get("thresholds") is None, "A-only diagnostics do not accept accuracy thresholds")
    config["analysis_intervals"] = ae.parse_intervals(config.get("analysis_intervals"), "analysis_intervals")
    for key in ("max_reference_gap_sec", "availability_hold_sec"):
        ae.seconds_ns(config.get("pairing", {}).get(key), key)
    estimates = config.get("estimates")
    ae.require(isinstance(estimates, list) and estimates, "at least one estimate is required")
    names = [e.get("name") for e in estimates]
    safe_name = lambda n: isinstance(n, str) and bool(n) and all(c.isalnum() or c in "_-" for c in n)
    ae.require(all(safe_name(n) for n in names) and len(set(names)) == len(names), "unique file-safe estimate names required")
    ae.require(all(not e.get("alignment") for e in estimates), "diagnostics fit one A-only transform; do not supply formal alignments")
    groups = config.get("common_sample_groups", [])
    ae.require(isinstance(groups, list) and groups, "declare a common-keyframe comparison group")
    group_names = set()
    for group in groups:
        name, members = group.get("name"), group.get("estimates")
        ae.require(safe_name(name) and name not in group_names, "unique file-safe group names required")
        ae.require(isinstance(members, list) and len(members) >= 2 and len(set(members)) == len(members)
                   and set(members) <= set(names), "group must contain distinct configured estimates")
        group_names.add(name)
    fit = config.get("fit", {})
    ae.require(fit.get("method") == "position", "diagnostic alignment must use position only; attitude is reported without attitude fitting")
    ae.require(config.get("fit_estimate") in names, "fit_estimate must select one configured trajectory")
    ae.require(config.get("fit_common_sample_group") in group_names, "fit_common_sample_group must select a declared group")
    selected_group = next(g for g in groups if g["name"] == config["fit_common_sample_group"])
    ae.require(config["fit_estimate"] in selected_group["estimates"], "fit estimate must belong to the common fit group")
    ae.require(isinstance(config.get("conditions", []), list)
               and all(isinstance(c, str) for c in config.get("conditions", [])), "conditions must be a list of strings")
    if config.get("relative_motion") is not None:
        relative = config["relative_motion"]
        deltas = relative.get("deltas_sec")
        ae.require(isinstance(deltas, list) and deltas, "relative_motion requires declared deltas_sec")
        delta_ns = [ae.seconds_ns(d, "relative_motion.deltas_sec") for d in deltas]
        ae.require(len(delta_ns) == len(set(delta_ns)), "relative motion deltas must be distinct")
        ae.seconds_ns(relative.get("max_estimate_gap_sec"), "relative_motion.max_estimate_gap_sec")
    config["_path"], config["_base"] = path.resolve(), path.resolve().parent
    return config


def relative_motion(config: dict, reference: ae.Reference, estimate: dict) -> tuple[list[dict], dict]:
    """Fixed-time RPE without global alignment or crossing an invalid interior.

    E = inv(inv(Tref0)*Tref1) * (inv(Test0)*Test1). Its translation
    components are expressed in the reference endpoint body frame.
    """
    options = config.get("relative_motion")
    if options is None:
        return [], {}
    gap_ns = ae.seconds_ns(options["max_estimate_gap_sec"], "relative_motion.max_estimate_gap_sec")
    native = ae.Reference(estimate["rows"], estimate["metadata"]["recorded_intervals"], gap_ns)

    def support_index(intervals):
        index = {}
        for span in intervals:
            starts, ends = index.setdefault(span["epoch"], ([], []))
            starts.append(span["start_ns"])
            ends.append(span["end_ns"])
        return index

    def supported(index, epoch, start, end):
        starts, ends = index.get(epoch, ([], []))
        i = bisect.bisect_right(starts, start)-1
        return i >= 0 and end <= ends[i]

    def relative_pose(a, b):
        r0, r1 = ae.quat_matrix(a.quaternion), ae.quat_matrix(b.quaternion)
        return ae.Pose(0, r0.T @ (b.position-a.position), ae.matrix_quat(r0.T @ r1), a.epoch)

    ref_support = support_index(reference.valid_intervals())
    est_support = support_index(ae.shifted(native.valid_intervals(), estimate["offset"], estimate["reference_epoch"]))
    starts = [row for row in estimate["rows"] if ae.interval_for(config["analysis_intervals"], estimate["reference_epoch"],
                                                               row.timestamp+estimate["offset"]) is not None]
    summaries, outputs = [], {}
    for delta in options["deltas_sec"]:
        delta_ns = ae.seconds_ns(delta, "relative_motion.deltas_sec")
        samples, errors, reasons, eligible = [], [], {}, 0
        for source in starts:
            start, end, epoch = source.timestamp+estimate["offset"], source.timestamp+estimate["offset"]+delta_ns, estimate["reference_epoch"]
            row = {"timestamp_ns": start, "end_timestamp_ns": end, "estimate_timestamp_ns": source.timestamp,
                   "epoch": epoch, "status": "paired", "detail": ""}
            interval = ae.interval_for(config["analysis_intervals"], epoch, start)
            if interval != ae.interval_for(config["analysis_intervals"], epoch, end):
                row["status"] = "endpoint_outside_analysis_interval"
            else:
                eligible += 1
                final, detail = native.at(source.epoch, source.timestamp+delta_ns)
                ref0, detail0 = reference.at(epoch, start)
                ref1, detail1 = reference.at(epoch, end)
                if final is None:
                    row["status"], row["detail"] = "estimate_endpoint_unavailable", detail
                elif not supported(est_support, epoch, start, end):
                    row["status"] = "estimate_segment_discontinuous"
                elif ref0 is None:
                    row["status"], row["detail"] = "reference_start_unavailable", detail0
                elif ref1 is None:
                    row["status"], row["detail"] = "reference_endpoint_unavailable", detail1
                elif not supported(ref_support, epoch, start, end):
                    row["status"] = "reference_segment_discontinuous"
                else:
                    error = ae.disagreement(relative_pose(source, final), relative_pose(ref0, ref1), np.eye(3), np.zeros(3))
                    values = dict(zip(RPE_METRICS, [error["position_3d_m"], error["longitudinal_m"],
                        error["lateral_m"], error["vertical_m"], error["attitude_geodesic_deg"]]))
                    row.update(values)
                    errors.append(values)
            if row["status"] != "paired":
                reasons[row["status"]] = reasons.get(row["status"], 0)+1
            samples.append(row)
        filename = f"{estimate['name']}_rpe_{delta_ns}ns_samples.csv"
        summaries.append({"delta_sec": delta_ns/ae.NS, "status": "CONDITIONAL_RELATIVE_REFERENCE_DISAGREEMENT" if errors else "NO_PAIRED_DATA", "pass": None,
            "global_alignment_applied": False, "max_estimate_gap_sec": gap_ns/ae.NS,
            "start_samples_in_analysis": len(starts), "eligible_start_count": eligible,
            "paired_count": len(errors), "rejected_eligible_count": eligible-len(errors),
            "paired_fraction_of_eligible": len(errors)/eligible if eligible else None,
            "outside_analysis_endpoint_count": reasons.get("endpoint_outside_analysis_interval", 0),
            "rejected_reasons": reasons,
            "metrics": {key: ae.statistics([r[key] for r in errors], signed=key.startswith("error_")) for key in RPE_METRICS},
            "sample_csv": filename,
            "denominator": "Every estimate start sample whose full fixed delta lies within one declared analysis interval, before reference gating.",
            "contract": "E=inv(Tref0^-1*Tref1)*(Test0^-1*Test1); no global fit. Both endpoints and the entire segment require continuous qualified reference and estimate support. Position-linear/SLERP interpolation only within the declared maximum gap.",
            "interpretation": "Conditional on body/lever/time/datum contracts. Sparse keyframes may yield zero eligible pairs; reconstructed dense optimization is supplemental. Overlapping time windows are not independent trials."})
        outputs[filename] = samples
    return summaries, outputs


def diagnose(config: dict, output: Path) -> dict:
    ae.require(config.get("role") == "diagnostic_mapping_a", "diagnose only accepts diagnostic_mapping_a")
    reference, ref_contract, estimates = ae.load_inputs(config)
    recording_id = ref_contract["metadata"]["source_recording_id"]
    frame = estimates[0]["metadata"]["world_frame_id"]
    for estimate in estimates:
        ae.require(estimate["metadata"]["source_recording_id"] == recording_id,
                   "A-only inputs must declare the same physical source_recording_id")
        ae.require(estimate["metadata"]["world_frame_id"] == frame, "all A trajectories must share one map/run world frame")
        ae.require(estimate["metadata"]["trajectory_kind"] != "localization",
                   "A mapping diagnostics cannot claim localization evaluation")
    by_name = {e["name"]: e for e in estimates}
    ae.require(by_name[config["fit_estimate"]]["metadata"]["trajectory_kind"] == "optimized_keyframes",
               "predeclare optimized_keyframes as the single diagnostic alignment source")
    matched = {e["name"]: ae.matched_rows(config, reference, e) for e in estimates}
    group = next(g for g in config["common_sample_groups"] if g["name"] == config["fit_common_sample_group"])
    pair_lookup = {name: {(b.epoch, b.timestamp): (a, b) for a, b in matched[name][0]} for name in group["estimates"]}
    fit_keys = sorted(set.intersection(*(set(v) for v in pair_lookup.values())))
    fit_pairs = [pair_lookup[config["fit_estimate"]][key] for key in fit_keys]
    rotation, translation, fit_diagnostics = ae.fit_rigid(fit_pairs, config["fit"])
    fit_mask_bytes = json.dumps(fit_keys, separators=(",", ":")).encode()
    alignment = {
        "schema_version": 1, "artifact_type": "mapping_a_diagnostic_alignment", "status": STATUS,
        "usable_for_formal_evaluate": False, "independent_calibration": False,
        "convention": "reference_from_estimate_world", "source_world_frame_id": frame,
        "target_world_frame_id": ref_contract["metadata"]["world_frame_id"],
        "body_frame_id": ae.BODY_FRAME,
        "transform": {"translation_m": translation.tolist(), "quaternion_xyzw": ae.matrix_quat(rotation).tolist(), "scale": 1.0},
        "fit_estimate": config["fit_estimate"], "fit_common_sample_group": group["name"],
        "fit_pair_count": len(fit_pairs), "fit_mask_sha256": hashlib.sha256(fit_mask_bytes).hexdigest(),
        "fit_mask_csv": "alignment_fit_samples.csv", "geometry": fit_diagnostics,
        "fit_reference_disagreement": ae.summarize([ae.disagreement(a, b, rotation, translation) for a, b in fit_pairs]),
        "claim": LIMITATION,
    }
    alignment_bytes = (json.dumps(alignment, ensure_ascii=False, indent=2, allow_nan=False)+"\n").encode()
    alignment_hash = hashlib.sha256(alignment_bytes).hexdigest()
    conditions = list(dict.fromkeys([LIMITATION] + config.get("conditions", [])
                                   + [c for e in estimates for c in e["conditions"]]))
    provenance = ae.provenance(config)
    provenance.update(tool_version=VERSION, tool_sha256=ae.sha256(Path(__file__).resolve()),
                      core_tool_sha256=ae.sha256(Path(ae.__file__).resolve()))
    report = {"schema_version": 1, "report_type": "mapping_a_aligned_diagnostic", "status": STATUS,
              "independent_accuracy_assessed": False, "localization_accuracy_assessed": False,
              "provenance": provenance, "reference": ref_contract, "analysis_intervals": config["analysis_intervals"],
              "diagnostic_alignment": {"file": "diagnostic_alignment.json", "sha256": alignment_hash, "content": alignment},
              "conditions": conditions, "interpretation": [LIMITATION,
                  "One position-only SE(3), unit scale, fixed declared timestamps; no curve-specific refit, time fitting, or lever-arm tuning.",
                  "Offline pose timestamp coverage does not measure live publication uptime or latency.",
                  "Reference validity is taken from the declared reference contract, never from the residual magnitude."],
              "estimates": [], "common_sample_groups": []}
    selected_reference = [r for r in reference.rows if ae.interval_for(config["analysis_intervals"], r.epoch, r.timestamp) is not None]
    reasons = {}
    for row in selected_reference:
        if not row.valid:
            reasons[row.reason] = reasons.get(row.reason, 0)+1
    report["reference_quality"] = {"input_rows": len(reference.rows), "rows_in_analysis_intervals": len(selected_reference),
                                   "valid_rows": sum(r.valid for r in selected_reference), "invalid_reasons": reasons}
    samples_by_name, coverage_by_name, rpe_outputs = {}, {}, {}
    reference_valid = reference.valid_intervals()
    for estimate in estimates:
        name = estimate["name"]
        pairs, unmatched = matched[name]
        samples = list(unmatched)
        for source, target in pairs:
            samples.append({"timestamp_ns": target.timestamp, "estimate_timestamp_ns": source.timestamp,
                            "epoch": target.epoch, "status": "paired", **ae.disagreement(source, target, rotation, translation)})
        samples.sort(key=lambda r: (r["epoch"], r["timestamp_ns"]))
        available = ae.estimate_availability(estimate["rows"], estimate["metadata"]["recorded_intervals"],
                                             estimate["hold_ns"], estimate["offset"], estimate["reference_epoch"])
        exposure, segments = ae.coverage(config["analysis_intervals"], reference_valid, available)
        reason_counts = {}
        for row in unmatched:
            reason_counts[row["status"]] = reason_counts.get(row["status"], 0)+1
        entry = {"name": name, "trajectory_kind": estimate["metadata"]["trajectory_kind"],
                 "estimate": estimate["contract"], "time": estimate["time"], "conditions": conditions,
                 "diagnostic_alignment_sha256": alignment_hash, "coverage": exposure,
                 "sample_counts": {"input_estimate": len(estimate["rows"]), "in_analysis_intervals": len(samples),
                                   "outside_analysis_intervals": len(estimate["rows"])-len(samples),
                                   "paired": len(pairs), "unpaired": len(unmatched), "unpaired_reasons": reason_counts},
                 "metrics": ae.summarize([s for s in samples if s["status"] == "paired"]),
                 "acceptance": {"status": STATUS, "pass": None, "checks": {}, "reason": LIMITATION},
                 "coverage_criterion": {"maximum_estimate_age_sec": estimate["hold_ns"]/ae.NS,
                                        "basis": "offline pose timestamps; sparse keyframes are not live availability"},
                 "segments": segments, "sample_csv": name+"_samples.csv", "coverage_csv": name+"_coverage.csv"}
        if entry["trajectory_kind"] == "optimized_dense_reconstructed":
            entry["interpretation"] = "Dense trajectory interpolates optimized keyframe corrections; it is supplemental, not independently optimized at each frame."
        entry["relative_motion"], relative_outputs = relative_motion(config, reference, estimate)
        rpe_outputs.update(relative_outputs)
        report["estimates"].append(entry)
        samples_by_name[name], coverage_by_name[name] = samples, segments
    group_outputs = {}
    for group in config["common_sample_groups"]:
        lookup = {name: {(r["epoch"], r["timestamp_ns"]): r for r in samples_by_name[name]} for name in group["estimates"]}
        keys = set.intersection(*(set(v) for v in lookup.values()))
        valid_keys = sorted(k for k in keys if all(v[k]["status"] == "paired" for v in lookup.values()))
        report["common_sample_groups"].append({"name": group["name"], "estimates": group["estimates"],
            "common_timestamp_count": len(keys), "paired_common_count": len(valid_keys),
            "metrics": {name: ae.summarize([lookup[name][k] for k in valid_keys]) for name in group["estimates"]},
            "sample_csv": "group_"+group["name"]+"_samples.csv",
            "contract": "Exact epoch/time intersection, same reference mask and single A-only alignment; no resampling or refitting."})
        group_outputs[group["name"]] = [{"estimate_name": name, **lookup[name][k]} for k in valid_keys for name in group["estimates"]]
    output.mkdir(parents=True, exist_ok=False)
    (output/"diagnostic_alignment.json").write_bytes(alignment_bytes)
    ae.write_csv(output/"alignment_fit_samples.csv", ["epoch", "timestamp_ns"],
                 [{"epoch": k[0], "timestamp_ns": k[1]} for k in fit_keys])
    for name, samples in samples_by_name.items():
        ae.write_csv(output/(name+"_samples.csv"), SAMPLE_FIELDS, samples)
        ae.write_csv(output/(name+"_coverage.csv"), ["epoch", "start_ns", "end_ns", "state", "duration_sec"], coverage_by_name[name])
    for name, rows in group_outputs.items():
        ae.write_csv(output/("group_"+name+"_samples.csv"), ["estimate_name"]+SAMPLE_FIELDS, rows)
    for filename, rows in rpe_outputs.items():
        ae.write_csv(output/filename, RPE_FIELDS, rows)
    with (output/"report.json").open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    (output/"report.md").write_text(markdown_report(report), encoding="utf-8")
    return report


def markdown_report(report: dict) -> str:
    lines = ["# 建图 A 同数据对齐诊断", "", "**A_ONLY_ALIGNED_DIAGNOSTIC：不是独立绝对精度或纯雷达定位精度验收。**", "",
             "在 A 的共同有效关键帧上，以最终优化关键帧预先指定为拟合源，仅拟合一次位置 SE(3)，随后所有曲线共用它。"
             "因此优化轨迹的拟合残差存在同数据乐观偏差，不能将前后差值直接解释为独立验证的优化收益。", "",
             f"同数据拟合样本数：{report['diagnostic_alignment']['content']['fit_pair_count']}；变换 SHA-256：`{report['diagnostic_alignment']['sha256']}`。", "",
             "| 轨迹 | 配对数 | 3D RMSE(m) | 水平 RMSE(m) | 姿态 SO(3) RMSE(deg) | 参考有效时长(s) | 完整分析时长(s) |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for entry in report["estimates"]:
        m, c = entry["metrics"], entry["coverage"]
        values = [entry["name"], str(entry["sample_counts"]["paired"]), ae.formatted(m["position_3d_m"]["rmse"]),
                  ae.formatted(m["horizontal_m"]["rmse"]), ae.formatted(m["attitude_geodesic_deg"]["rmse"]),
                  ae.formatted(c["reference_valid_duration_sec"]), ae.formatted(c["recorded_duration_sec"])]
        lines.append("| "+" | ".join(values)+" |")
    lines += ["", "完整三轴位置、FLU位置、roll/pitch/yaw与SO(3)分布、共同关键帧统计及时间覆盖见 `report.json` 和逐样本CSV。"
              "参考无效与估计缺失均进入完整区间分母；稀疏关键帧覆盖不等于实时定位可用率。", ""]
    if any(e.get("relative_motion") for e in report["estimates"]):
        lines += ["固定时距相对位姿差（RPE）不使用A拟合变换；要求两端及整个区间连续有效。平移误差三轴按标准SE(3)残差在参考末端FLU表达。"
                  "非均匀采样只在声明最大gap内线性位置/SLERP插值；稀疏关键帧可能零对。各路按自身真实起始样本加权，不同样本集合不能直接宣称优化收益，重叠窗口不属于独立重复实验。", "",
                  "| 轨迹 | 时距(s) | 分析区间起始样本 | 完整时距候选 | 配对 | 拒配候选 | 平移RMSE(m) | SO(3)RMSE(deg) | 拒配原因 |",
                  "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
        for entry in report["estimates"]:
            for relative in entry.get("relative_motion", []):
                values = [entry["name"], ae.formatted(relative["delta_sec"]), str(relative["start_samples_in_analysis"]),
                    str(relative["eligible_start_count"]), str(relative["paired_count"]), str(relative["rejected_eligible_count"]),
                    ae.formatted(relative["metrics"]["translation_3d_m"]["rmse"]), ae.formatted(relative["metrics"]["rotation_geodesic_deg"]["rmse"]),
                    "; ".join(f"{key}:{count}" for key, count in relative["rejected_reasons"].items()) or "none"]
                lines.append("| "+" | ".join(values)+" |")
        lines.append("")
    lines += ["条件与未确认项：", ""]
    lines += ["- "+c for c in report["conditions"]]
    return "\n".join(lines)+"\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        result = diagnose(load_config(args.config), args.output)
        print(json.dumps({"status": result["status"], "output": str(args.output.resolve()), "independent_accuracy_assessed": False}))
        return 0
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"diagnostic error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
