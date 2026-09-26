#!/usr/bin/env python3
"""Summarize production health evidence without treating internal validity as ROS delivery."""
import argparse
import collections
import datetime
import hashlib
import json
import math
import pathlib
import re

ANSI = re.compile(r"\x1b\[[0-9;]*m")
HEADER = re.compile(r"[IWEF](\d{4}) (\d{2}:\d{2}:\d{2}\.\d+) .*? ([\w.]+:\d+)\] (.*)")
SCALAR = re.compile(r"(\w+)=([^\s\[]+)")
VECTOR = re.compile(r"(\w+)=\[([^\]]+)\]")


def analyze(run):
    log = run / "logs/run_loc_online.stderr.log"
    if not log.is_file():
        return {"run": run.name, "missing_algorithm_log": True}
    data = log.read_bytes()
    sources = collections.Counter()
    events, audits, pgo, warnings = [], [], [], []
    first, last = None, None
    for number, raw in enumerate(data.decode("utf-8", errors="replace").splitlines(), 1):
        line = ANSI.sub("", raw)
        match = HEADER.search(line)
        if not match:
            continue
        date, clock, source, message = match.groups()
        first = first or clock
        last = clock
        sources[source] += 1
        entry = {"line": number, "time": clock, "message": message}
        if any(term in message for term in (
            "accepted:", "tracking lost after", "outputs disabled", "localization failed!",
        )):
            events.append(entry)
        if "smoother motion is too large:" in message:
            warnings.append(entry)
        if "LIDAR_LOC_AUDIT" in message:
            values = dict(SCALAR.findall(message))
            vectors = {key: list(map(float, value.split())) for key, value in VECTOR.findall(message)}
            # A forcibly terminated run can end with a partial audit line.
            if not all(key in values for key in ("confidence", "published_valid")):
                continue
            if not all(key in vectors for key in ("guess_xyz", "ndt_xyz", "balanced_xyz")):
                continue
            audits.append({"line": number, "time": clock, **values, **vectors,
                           "raw_correction_m": math.dist(vectors["guess_xyz"], vectors["ndt_xyz"])})
        if "PGO_OUTPUT_AUDIT" in message:
            values = dict(SCALAR.findall(message))
            if "output_lidar_delta_xy" in values:
                pgo.append({**entry, "delta_xy_m": float(values["output_lidar_delta_xy"])})
    feedback = {time: next((a for a in audits if a["time"].startswith(time)), None)
                for time in ("16:15:00", "17:54:00")}
    high_frequency = run / "results/trajectory_high_frequency.tum"
    nearest = {}
    if high_frequency.exists():
        date = re.search(r"_(\d{8})_", run.name).group(1)
        for time, audit in feedback.items():
            if audit:
                wall = datetime.datetime.strptime(date + time, "%Y%m%d%H:%M:%S")
                target = wall.replace(tzinfo=datetime.timezone(datetime.timedelta(hours=8))).timestamp()
                nearest[time] = {"target_unix": target, "offset_sec": float("inf")}
        if nearest:
            with high_frequency.open() as stream:
                for line in stream:
                    if not line.strip() or line.startswith("#"):
                        continue
                    values = list(map(float, line.split()))
                    for time, best in nearest.items():
                        offset = abs(values[0] - best["target_unix"])
                        if offset < best["offset_sec"]:
                            best.update(offset_sec=offset, timestamp=values[0], xyz=values[1:4])
            for time, best in nearest.items():
                if "xyz" in best:
                    best["sampled_lidar_disagreement_xy_m"] = math.dist(
                        best["xyz"][:2], feedback[time]["balanced_xyz"][:2])
    published = run / "results/trajectory_published_rear_axle.tum"
    recorded_count = (sum(1 for line in published.open() if line.strip() and not line.startswith("#"))
                      if published.exists() else None)
    return {"run": run.name, "log": str(log), "sha256": hashlib.sha256(data).hexdigest(),
            "first_log_time": first, "last_log_time": last,
            "sampled_matches": len(audits),
            "sampled_invalid_matches": sum(a["published_valid"] == "0" for a in audits),
            "minimum_sampled_confidence": min((float(a["confidence"]) for a in audits), default=None),
            "max_raw_correction": max(audits, key=lambda a: a["raw_correction_m"], default=None),
            "max_pgo_lidar_disagreement": max(pgo, key=lambda a: a["delta_xy_m"], default=None),
            "feedback_samples": feedback, "nearest_internal_high_frequency_tum": nearest, "events": events,
            "smoother_motion_rejections": len(warnings),
            "first_smoother_motion_rejection": warnings[0] if warnings else None,
            "last_smoother_motion_rejection": warnings[-1] if warnings else None,
            "recorded_published_poses": recorded_count,
            "top_log_sources": sources.most_common(12)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("diagnostics", type=pathlib.Path)
    parser.add_argument("--date", default="20260925")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    result = {"time_basis": "Log wall clock, Asia/Shanghai; feedback time is not failure onset.",
              "limitation": "No raw bags or ground truth. published_valid is an internal LidarLoc flag, not a ROS receipt.",
              "runs": [analyze(run) for run in sorted(args.diagnostics.glob(f"*{args.date}*")) if run.is_dir()]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for run in result["runs"]:
        print(run["run"], run.get("first_log_time"), run.get("last_log_time"),
              "sampled_invalid=", run.get("sampled_invalid_matches"),
              "smoother_rejections=", run.get("smoother_motion_rejections"))


if __name__ == "__main__":
    main()
