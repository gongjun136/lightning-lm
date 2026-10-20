#!/usr/bin/env python3
"""Read-only correlation of safety transitions, clock steps and sampled load.

Log durations are wall-clock intervals, not physical elapsed durations when the
clock steps. Log rows are not uniformly sampled ROS messages or ground truth.
"""
import argparse
import bisect
import collections
import csv
import datetime
import json
from pathlib import Path
import re


def quantiles(values):
    values = sorted(values)
    return {str(q): values[int((len(values) - 1) * q)]
            for q in (0.5, 0.95, 0.99, 1)} if values else {}


def analyze(root):
    metadata = dict(line.split("=", 1) for line in
                    (root / "run_metadata.txt").read_text().splitlines() if "=" in line)
    date = datetime.datetime.fromisoformat(metadata["started_at"].replace(",", "."))
    midnight = date.replace(hour=0, minute=0, second=0, microsecond=0).timestamp()
    previous = None
    steps = []
    malformed = 0
    with (root / "results/causal_trace.csv").open() as stream:
        for row in csv.DictReader(stream):
            if row["sequence"].startswith("#"):
                continue
            try:
                wall = int(row["wall_ns"]) / 1e9
                steady = int(row["steady_ns"]) / 1e9
            except (ValueError, TypeError):
                malformed += 1
                continue
            if previous and steady > previous[1]:
                delta = (wall - steady) - (previous[0] - previous[1])
                if abs(delta) > 0.1:
                    steps.append(dict(wall=wall, steady=steady, delta_sec=delta,
                                      sequence=row["sequence"]))
            previous = wall, steady
    step_times = [step["wall"] for step in steps]
    rows = collections.Counter()
    stop_conditions = collections.Counter()
    runs = []
    with (root / "logs/run_loc_online.stderr.log").open(errors="replace") as stream:
        for line_number, line in enumerate(stream, 1):
            if "SAFETY_TIME_SYNC_SHADOW" not in line:
                continue
            fields = dict(re.findall(r"(\w+)=(-?[\d.]+)", line))
            state = int(fields["published_runstate"])
            rows[state] += 1
            if state == 4:
                stop_conditions[fields["localization_consistency"] + "/" +
                                fields["localization_uncertainty"]] += 1
            time = line.split()[1]
            h, m, s = map(float, time.split(":"))
            wall = midnight + h * 3600 + m * 60 + s
            if not runs or runs[-1]["state"] != state:
                if runs:
                    runs[-1]["wall_duration_sec"] = wall - runs[-1]["wall"]
                index = bisect.bisect_right(step_times, wall) - 1
                runs.append(dict(state=state, wall=wall, time=time, line=line_number,
                                 since_clock_step_sec=wall - step_times[index] if index >= 0 else None))
    cores = []
    host = collections.defaultdict(list)
    with (root / "results/process_resources.jsonl").open() as stream:
        for line in stream:
            record = json.loads(line)
            if record["type"] != "sample":
                continue
            if record.get("cpu_core_equivalents") is not None:
                cores.append(record["cpu_core_equivalents"])
            for cpu in record.get("host_cpus", []):
                if cpu["busy_percent"] is not None:
                    host[cpu["cpu"]].append(cpu["busy_percent"])
    return dict(
        git_commit=metadata.get("git_commit"),
        caveats=["Safety log rows are change-triggered plus periodic; counts are not time percentages.",
                 "Wall durations include clock steps; do not use them as a grace-period calibration.",
                 "Clock correlation alone does not identify the service that changed time.",
                 "Sampled CPU busy time is not proof of scheduler starvation."],
        clock_steps=steps, malformed_trace_rows=malformed,
        safety_log_rows=dict(rows), stop_conditions=dict(stop_conditions),
        state_run_counts=dict(collections.Counter(r["state"] for r in runs)),
        severe_runs=[r for r in runs if r["state"] == -9],
        process_cpu_cores=quantiles(cores),
        host_cpu_busy_percent={cpu: quantiles(values) for cpu, values in host.items()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze(args.run), indent=2, ensure_ascii=False))
