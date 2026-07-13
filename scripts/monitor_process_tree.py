#!/usr/bin/env python3
"""Sample CPU and RSS for a Linux process tree until a stop file appears."""

import argparse
import csv
import json
import os
import time


def process_snapshot(root_pid):
    stats = {}
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        pid = int(name)
        try:
            with open(f"/proc/{pid}/stat", "r", encoding="utf-8") as stream:
                content = stream.read()
            closing = content.rfind(")")
            if closing < 0:
                continue
            fields = content[closing + 2:].split()
            ppid = int(fields[1])
            ticks = int(fields[11]) + int(fields[12])
            starttime = int(fields[19])
            with open(f"/proc/{pid}/statm", "r", encoding="utf-8") as stream:
                resident_pages = int(stream.read().split()[1])
            stats[pid] = (ppid, ticks, resident_pages, starttime)
        except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError, IndexError):
            continue
    descendants = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, (ppid, _, _, _) in stats.items():
            if ppid in descendants and pid not in descendants:
                descendants.add(pid)
                changed = True
    return {pid: stats[pid] for pid in descendants if pid in stats}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pid", type=int, required=True)
    parser.add_argument("--stop-file", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--allocated-cpus", type=int, default=8)
    parser.add_argument("--interval", type=float, default=0.2)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.csv)), exist_ok=True)
    ticks_per_second = os.sysconf(os.sysconf_names["SC_CLK_TCK"])
    page_size = os.sysconf("SC_PAGE_SIZE")
    start = time.monotonic()
    previous_time = start
    previous_ticks = {}
    samples = []

    with open(args.csv, "w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=("elapsed_s", "processes", "cpu_cores", "cpu_pct_of_allocation", "rss_mb"))
        writer.writeheader()
        while not os.path.exists(args.stop_file):
            now = time.monotonic()
            snapshot = process_snapshot(args.pid)
            rss_mb = sum(value[2] for value in snapshot.values()) * page_size / (1024.0 * 1024.0)
            current_ticks = {(pid, value[3]): value[1] for pid, value in snapshot.items()}
            if previous_ticks and now > previous_time:
                delta_ticks = sum(max(0, ticks - previous_ticks.get(key, ticks)) for key, ticks in current_ticks.items())
                cpu_cores = delta_ticks / ticks_per_second / (now - previous_time)
                row = {
                    "elapsed_s": now - start,
                    "processes": len(snapshot),
                    "cpu_cores": max(0.0, cpu_cores),
                    "cpu_pct_of_allocation": max(0.0, 100.0 * cpu_cores / args.allocated_cpus),
                    "rss_mb": rss_mb,
                }
                writer.writerow(row)
                output.flush()
                samples.append(row)
            previous_ticks = current_ticks
            previous_time = now
            time.sleep(args.interval)

    def percentile(values, fraction):
        if not values:
            return 0.0
        ordered = sorted(values)
        return ordered[min(len(ordered) - 1, round((len(ordered) - 1) * fraction))]

    summary = {
        "samples": len(samples),
        "duration_s": time.monotonic() - start,
        "mean_cpu_cores": sum(row["cpu_cores"] for row in samples) / len(samples) if samples else 0.0,
        "peak_cpu_cores": max((row["cpu_cores"] for row in samples), default=0.0),
        "p95_cpu_cores": percentile([row["cpu_cores"] for row in samples], 0.95),
        "mean_cpu_pct_of_allocation": sum(row["cpu_pct_of_allocation"] for row in samples) / len(samples) if samples else 0.0,
        "peak_cpu_pct_of_allocation": max((row["cpu_pct_of_allocation"] for row in samples), default=0.0),
        "mean_rss_mb": sum(row["rss_mb"] for row in samples) / len(samples) if samples else 0.0,
        "peak_rss_mb": max((row["rss_mb"] for row in samples), default=0.0),
    }
    with open(args.summary, "w", encoding="utf-8") as output:
        json.dump(summary, output, indent=2)


if __name__ == "__main__":
    main()
