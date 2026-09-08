#!/usr/bin/env python3
"""Compare replay artifacts with a frozen baseline (relative, not ground truth)."""
import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np


def distribution(values):
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if not len(a):
        return {"count": 0}
    return dict(count=len(a), mean=float(a.mean()), p5=float(np.percentile(a, 5)), p50=float(np.percentile(a, 50)),
                p95=float(np.percentile(a, 95)), p99=float(np.percentile(a, 99)), max=float(a.max()))


def read_run(path):
    lio = []
    ndt = []
    with (path / "stderr.log").open(errors="replace") as log:
        for line in log:
            is_lio = "LIO_BENCH_FRAME" in line
            is_ndt = "COMPUTE_BENCH_FRAME module=lidar_loc " in line
            if (is_lio or is_ndt) and "phase=tracking" in line:
                fields = dict(re.findall(r"(\w+)=([^\s]+)", line))
                numeric = {}
                for key, value in fields.items():
                    try:
                        numeric[key] = int(value) if key == "frame_id" else float(value)
                    except ValueError:
                        pass
                (lio if is_lio else ndt).append(numeric)
    with (path / "frames.csv").open() as source:
        frames = list(csv.DictReader(source))
    trajectory = np.atleast_2d(np.loadtxt(path / "ndt.tum"))
    return lio, frames, trajectory, ndt


def rotations(q):
    q = q / np.linalg.norm(q, axis=1)[:, None]
    x, y, z, w = q.T
    return np.stack((1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w),
                     2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w),
                     2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)), axis=1).reshape(-1, 3, 3)


def summarize(lio, frames, ndt):
    numeric = ("total_ms", "scan_match_ms", "downsample_ms", "output_points", "effective_features",
               "obs_evaluations", "selected_lidars", "point_budget")
    result = {key: distribution([row[key] for row in lio if key in row]) for key in numeric}
    result["lio_gt_100ms"] = sum(row["total_ms"] > 100 for row in lio)
    result["tracking_unhealthy"] = sum(row.get("tracking_healthy", 1) == 0 for row in lio)
    for key in ("ndt_align_wall_ms", "ndt_align_process_cpu_ms_concurrent", "ndt_input_points"):
        result[key] = distribution([row[key] for row in ndt if key in row])
    if frames:
        result["frames"] = len(frames)
        first_valid = next((i for i, row in enumerate(frames) if row.get("valid") == "1"), len(frames))
        result["post_init_invalid_frames"] = sum(row.get("valid") != "1" for row in frames[first_valid:])
        result["post_init_failed_matches"] = sum(row.get("match_success") != "1" for row in frames[first_valid:])
        tracking = [row for row in frames if row.get("status") == "2"
                    and row.get("relocalization_attempted") == "0"]
        result["ndt_tracking_ms"] = distribution([float(row["processing_ms"]) for row in tracking])
        for key in ("processing_ms", "confidence", "loc_valid", "ndt_success", "valid", "match_success"):
            if key in frames[0]:
                result[key] = distribution([float(row[key]) for row in frames])
    return result


def compare(reference, candidate):
    # Sensor timestamps should match without a fitted spatial transformation.
    insertion = np.searchsorted(reference[:, 0], candidate[:, 0])
    right = np.clip(insertion, 0, len(reference)-1)
    left = np.clip(insertion-1, 0, len(reference)-1)
    indices = np.where(abs(reference[left, 0]-candidate[:, 0]) <
                       abs(reference[right, 0]-candidate[:, 0]), left, right)
    mask = abs(reference[indices, 0]-candidate[:, 0]) <= 0.001
    ref, cand = reference[indices[mask]], candidate[mask]
    if not len(ref):
        raise ValueError("no trajectories associate within 1 ms")
    if np.any(np.diff(reference[:, 0]) <= 0) or np.any(np.diff(candidate[:, 0]) <= 0):
        raise ValueError("trajectory timestamps must be strictly increasing")
    if len(np.unique(indices[mask])) != len(ref):
        raise ValueError("trajectory association is not one-to-one")
    rr, rc = rotations(ref[:, 4:8]), rotations(cand[:, 4:8])
    yaw_r = np.arctan2(rr[:, 1, 0], rr[:, 0, 0])
    yaw_c = np.arctan2(rc[:, 1, 0], rc[:, 0, 0])
    yaw_error = np.abs(np.rad2deg(np.arctan2(np.sin(yaw_c-yaw_r), np.cos(yaw_c-yaw_r))))
    result = dict(associated=len(ref), reference_rows=len(reference), candidate_rows=len(candidate),
                  coverage=len(ref)/len(reference),
                  translation_m=distribution(np.linalg.norm(cand[:, 1:4]-ref[:, 1:4], axis=1)),
                  yaw_deg=distribution(yaw_error))
    for seconds in (1, 10):
        future = np.searchsorted(ref[:, 0], ref[:, 0]+seconds)
        right = np.clip(future, 0, len(ref)-1)
        left = np.clip(future-1, 0, len(ref)-1)
        future = np.where(abs(ref[left, 0]-ref[:, 0]-seconds) <
                          abs(ref[right, 0]-ref[:, 0]-seconds), left, right)
        valid = abs(ref[future, 0]-ref[:, 0]-seconds) <= 0.005
        i = np.flatnonzero(valid)
        j = future[valid]
        tr = np.einsum("nji,nj->ni", rr[i], ref[j, 1:4]-ref[i, 1:4])
        tc = np.einsum("nji,nj->ni", rc[i], cand[j, 1:4]-cand[i, 1:4])
        dr = rr[i].transpose(0, 2, 1) @ rr[j]
        dc = rc[i].transpose(0, 2, 1) @ rc[j]
        difference = dr.transpose(0, 2, 1) @ dc
        angle = np.rad2deg(np.arccos(np.clip((np.trace(difference, axis1=1, axis2=2)-1)/2, -1, 1)))
        result[f"rpe_{seconds}s"] = dict(translation_m=distribution(np.linalg.norm(tc-tr, axis=1)),
                                          rotation_deg=distribution(angle))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--enforce-relative-gates", action="store_true")
    parser.add_argument("--expect-point-budget", type=int, default=0)
    parser.add_argument("--expect-lidar-count", type=int, default=0)
    args = parser.parse_args()
    base_lio, base_frames, base_trajectory, base_ndt = read_run(args.baseline)
    cand_lio, cand_frames, cand_trajectory, cand_ndt = read_run(args.candidate)
    report = dict(baseline=str(args.baseline), candidate=str(args.candidate),
                  limitation="WSL offline relative non-regression; not Orin full-stack realtime or absolute accuracy",
                  baseline_metrics=summarize(base_lio, base_frames, base_ndt),
                  candidate_metrics=summarize(cand_lio, cand_frames, cand_ndt),
                  trajectory=compare(base_trajectory, cand_trajectory))
    t = report["trajectory"]
    gates = dict(coverage=t["coverage"] >= 0.999,
                 position_p95=t["translation_m"]["p95"] <= 0.15,
                 position_p99=t["translation_m"]["p99"] <= 0.30,
                 yaw_p95=t["yaw_deg"]["p95"] <= 0.3,
                 yaw_p99=t["yaw_deg"]["p99"] <= 0.8)
    for seconds, meters, degrees in ((1, 0.05, 0.2), (10, 0.15, 0.5)):
        for field, limit in (("translation_m", meters), ("rotation_deg", degrees)):
            value = t[f"rpe_{seconds}s"][field]
            gates[f"rpe_{seconds}s_{field}"] = value.get("count", 0) > 0 and value["p95"] <= limit
    for metric in ("post_init_invalid_frames", "post_init_failed_matches"):
        gates[f"no_new_{metric}"] = report["candidate_metrics"].get(metric, float("inf")) <= \
            report["baseline_metrics"].get(metric, -1)
    report["relative_gates"] = gates
    report["relative_nonregression_pass"] = all(gates.values())
    activation = {}
    if args.expect_point_budget:
        capped = [row for row in cand_lio if row.get("point_budget", 0) > 0]
        activation["budget_active_fraction"] = len(capped) / max(1, len(cand_lio))
        activation["hard_cap_observed"] = bool(capped) and all(
            row["output_points"] <= row["point_budget"] <= args.expect_point_budget for row in capped)
        activation["budget_path_pass"] = activation["budget_active_fraction"] >= 0.95 and activation["hard_cap_observed"]
    if args.expect_lidar_count:
        activation["expected_lidar_fraction"] = sum(
            row.get("selected_lidars") == args.expect_lidar_count for row in cand_lio) / max(1, len(cand_lio))
        activation["lidar_path_pass"] = activation["expected_lidar_fraction"] >= 0.95
    report["activation"] = activation
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if args.enforce_relative_gates and (not report["relative_nonregression_pass"] or
            not all(value for key, value in activation.items() if key.endswith("_pass"))):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
