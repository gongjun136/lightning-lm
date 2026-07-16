#!/usr/bin/env python3
"""Print per-topic message counts and storage timestamp ranges for a ROS2 SQLite bag."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bag", type=Path)
    args = parser.parse_args()
    databases = sorted(args.bag.glob("*.db3")) if args.bag.is_dir() else [args.bag]
    if not databases:
        raise ValueError(f"no db3 payload found: {args.bag}")

    topics: dict[tuple[str, str], dict[str, int | str | None]] = {}
    for database in databases:
        connection = sqlite3.connect(f"file:{database}?mode=ro", uri=True)
        try:
            rows = connection.execute(
                """
                SELECT topics.name, topics.type, COUNT(messages.id),
                       MIN(messages.timestamp), MAX(messages.timestamp)
                FROM topics
                LEFT JOIN messages ON messages.topic_id = topics.id
                GROUP BY topics.id
                ORDER BY topics.name
                """
            )
            for name, message_type, count, minimum, maximum in rows:
                key = (name, message_type)
                summary = topics.setdefault(
                    key,
                    {
                        "name": name,
                        "type": message_type,
                        "message_count": 0,
                        "start_ns": None,
                        "end_ns": None,
                    },
                )
                summary["message_count"] = int(summary["message_count"]) + int(count)
                if minimum is not None:
                    summary["start_ns"] = (
                        int(minimum)
                        if summary["start_ns"] is None
                        else min(int(summary["start_ns"]), int(minimum))
                    )
                if maximum is not None:
                    summary["end_ns"] = (
                        int(maximum)
                        if summary["end_ns"] is None
                        else max(int(summary["end_ns"]), int(maximum))
                    )
        finally:
            connection.close()

    output = {"bag": str(args.bag.resolve()), "payload_files": len(databases), "topics": []}
    for summary in sorted(topics.values(), key=lambda item: str(item["name"])):
        start_ns = summary["start_ns"]
        end_ns = summary["end_ns"]
        summary["start_s"] = None if start_ns is None else int(start_ns) * 1e-9
        summary["end_s"] = None if end_ns is None else int(end_ns) * 1e-9
        summary["duration_s"] = (
            None if start_ns is None or end_ns is None else (int(end_ns) - int(start_ns)) * 1e-9
        )
        output["topics"].append(summary)
    print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
