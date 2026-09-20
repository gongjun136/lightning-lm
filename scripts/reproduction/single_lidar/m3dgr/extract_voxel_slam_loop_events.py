#!/usr/bin/env python3
"""Convert Voxel-SLAM accepted-edge logs to the common loop-event CSV schema."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FIND_LOOP = re.compile(r"Find Loop in session\d+:\s+(\d+)\s+(\d+)")
SCORE = re.compile(r"score:\s+([-+0-9.eE]+)")
ADD_EDGE = re.compile(r"addedge:\s+\(\d+\s+\d+\)\s+\((\d+)\s+(\d+)\)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def load_timestamps(path: Path) -> list[float]:
    timestamps: list[float] = []
    with path.open("r", encoding="utf-8-sig") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 8:
                timestamps.append(float(fields[0]))
    if not timestamps:
        raise ValueError(f"trajectory is empty: {path}")
    return timestamps


def parse_edges(path: Path) -> dict[int, dict[str, float | int]]:
    candidates: dict[int, dict[str, float | int]] = {}
    pending_current: int | None = None
    pending_score = 0.0
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = FIND_LOOP.search(line)
            if match:
                # Voxel-SLAM prints buf_base, while its graph vertex is buf_base - 1.
                pending_current = int(match.group(1)) - 1
                pending_score = 0.0
                continue
            match = SCORE.search(line)
            if match and pending_current is not None:
                pending_score = float(match.group(1))
                continue
            match = ADD_EDGE.search(line)
            if not match:
                continue
            history = int(match.group(1))
            current = int(match.group(2))
            score = pending_score if pending_current == current else 0.0
            candidates[current] = {
                "history": history,
                "score": score,
            }
            pending_current = None
            pending_score = 0.0
    return candidates


def main() -> None:
    args = parse_args()
    timestamps = load_timestamps(args.trajectory)
    edges = parse_edges(args.log)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "current_descriptor",
        "history_descriptor",
        "current_timestamp",
        "history_timestamp",
        "candidate",
        "accepted",
        "score",
        "reason",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for current, timestamp in enumerate(timestamps):
            edge = edges.get(current)
            history = int(edge["history"]) if edge else -1
            valid = edge is not None and 0 <= history < len(timestamps)
            writer.writerow(
                {
                    "current_descriptor": current,
                    "history_descriptor": history if valid else -1,
                    "current_timestamp": f"{timestamp:.9f}",
                    "history_timestamp": f"{timestamps[history]:.9f}" if valid else "0",
                    "candidate": int(valid),
                    "accepted": int(valid),
                    "score": f"{float(edge['score']):.12g}" if valid else "0",
                    "reason": "voxel_slam_addedge" if valid else "no_accepted_edge",
                }
            )
    print(args.output.resolve())
    print(f"accepted_edges={sum(1 for edge in edges if 0 <= edge < len(timestamps))}")


if __name__ == "__main__":
    main()
