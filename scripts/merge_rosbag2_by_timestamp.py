#!/usr/bin/env python3
"""Bulk-merge disjoint-topic ROS 2 sqlite bags without decoding CDR data."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sqlite3
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repair_metadata(output: Path, update_contract: bool = False) -> None:
    metadata = output / "metadata.yaml"
    text = metadata.read_text(encoding="utf-8")
    repaired = text.replace('  storage_identifier: ""', "  storage_identifier: sqlite3", 1)
    if repaired != text:
        metadata.write_text(repaired, encoding="utf-8")
    if update_contract:
        contract_path = output / "merge_contract.json"
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        contract["metadata_sha256"] = sha256(metadata)
        contract_path.write_text(
            json.dumps(contract, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )


def open_read_only(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)


def table_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [row[1] for row in connection.execute(f"PRAGMA table_info({table})")]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, type=Path,
                        help="input bag directory; repeat for every source")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--staging-directory", type=Path,
                        help="optional fast local directory used before one sequential copy to output")
    args = parser.parse_args()

    inputs = [path.resolve() for path in args.input]
    output = args.output.resolve()
    if len(inputs) < 2:
        raise SystemExit("at least two --input bags are required")
    if len(set(inputs)) != len(inputs):
        raise SystemExit("duplicate input bag")
    payloads: list[Path] = []
    for path in inputs:
        if not (path / "metadata.yaml").is_file():
            raise SystemExit(f"input is not a ROS 2 bag: {path}")
        files = sorted(path.glob("*.db3"))
        if len(files) != 1:
            raise SystemExit(f"exactly one sqlite payload is required in {path}; found {len(files)}")
        payloads.append(files[0])
    if output.exists():
        raise SystemExit(f"refusing to overwrite: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working_output = output
    if args.staging_directory:
        staging_root = args.staging_directory.resolve()
        staging_root.mkdir(parents=True, exist_ok=True)
        working_output = staging_root / f"{output.name}_staging"
        if working_output.exists():
            raise SystemExit(f"refusing to overwrite staging output: {working_output}")
    working_output.mkdir()
    output_db = working_output / f"{output.name}_0.db3"

    source_topics: list[list[dict[str, object]]] = []
    source_counts: list[dict[str, int]] = []
    topic_definitions: dict[str, tuple[object, ...]] = {}
    first = open_read_only(payloads[0])
    topic_columns = table_columns(first, "topics")
    message_columns = table_columns(first, "messages")
    if "id" not in topic_columns or "name" not in topic_columns:
        raise SystemExit("unsupported ROS 2 topics table")
    if not {"id", "topic_id", "timestamp", "data"}.issubset(message_columns):
        raise SystemExit("unsupported ROS 2 messages table")
    table_sql = [
        row[0] for row in first.execute(
            "SELECT sql FROM sqlite_master WHERE type='table' AND sql IS NOT NULL "
            "AND name NOT LIKE 'sqlite_%' ORDER BY name"
        )
    ]
    index_sql = [
        row[0] for row in first.execute(
            "SELECT sql FROM sqlite_master WHERE type='index' AND sql IS NOT NULL ORDER BY name"
        )
    ]
    page_size = int(first.execute("PRAGMA page_size").fetchone()[0])
    first.close()

    for input_index, payload in enumerate(payloads):
        connection = open_read_only(payload)
        if table_columns(connection, "topics") != topic_columns or \
           table_columns(connection, "messages") != message_columns:
            raise SystemExit(f"sqlite schema differs in {payload}")
        connection.row_factory = sqlite3.Row
        topics = [dict(row) for row in connection.execute("SELECT * FROM topics ORDER BY id")]
        counts_by_id = {
            int(topic_id): int(count)
            for topic_id, count in connection.execute(
                "SELECT topic_id, COUNT(*) FROM messages GROUP BY topic_id"
            )
        }
        named_counts: dict[str, int] = {}
        for topic in topics:
            name = str(topic["name"])
            signature = tuple(topic[column] for column in topic_columns if column not in {"id", "name"})
            if name in topic_definitions:
                raise SystemExit(
                    f"topic {name} occurs in more than one input; this merger requires disjoint topics"
                )
            topic_definitions[name] = signature
            named_counts[name] = counts_by_id.get(int(topic["id"]), 0)
            topic["source_index"] = input_index
        source_topics.append(topics)
        source_counts.append(named_counts)
        connection.close()

    output_connection = sqlite3.connect(output_db)
    output_connection.execute(f"PRAGMA page_size={page_size}")
    output_connection.execute("PRAGMA journal_mode=OFF")
    output_connection.execute("PRAGMA synchronous=OFF")
    output_connection.execute("PRAGMA temp_store=MEMORY")
    output_connection.execute("PRAGMA cache_size=-262144")
    for statement in table_sql:
        output_connection.execute(statement)

    # Copy the storage schema version, if present. Internal metadata is not
    # copied because ros2 bag reindex regenerates authoritative bag metadata.
    source = open_read_only(payloads[0])
    existing_tables = {
        row[0] for row in source.execute("SELECT name FROM sqlite_master WHERE type='table'")
    }
    if "schema" in existing_tables:
        schema_columns = table_columns(source, "schema")
        schema_rows = list(source.execute("SELECT * FROM schema"))
        placeholders = ",".join("?" for _ in schema_columns)
        output_connection.executemany(
            f"INSERT INTO schema ({','.join(schema_columns)}) VALUES ({placeholders})", schema_rows
        )
    source.close()

    next_topic_id = 1
    id_maps: list[dict[int, int]] = []
    insert_topic_sql = (
        f"INSERT INTO topics ({','.join(topic_columns)}) VALUES "
        f"({','.join('?' for _ in topic_columns)})"
    )
    for topics in source_topics:
        id_map: dict[int, int] = {}
        for topic in topics:
            old_id = int(topic["id"])
            new_id = next_topic_id
            next_topic_id += 1
            values = [new_id if column == "id" else topic[column] for column in topic_columns]
            output_connection.execute(insert_topic_sql, values)
            id_map[old_id] = new_id
        id_maps.append(id_map)
    output_connection.commit()

    for source_index, payload in enumerate(payloads):
        id_map = id_maps[source_index]
        case_expression = "CASE topic_id " + " ".join(
            f"WHEN {old_id} THEN {new_id}" for old_id, new_id in sorted(id_map.items())
        ) + " END"
        output_connection.execute("ATTACH DATABASE ? AS source", (str(payload),))
        output_connection.execute("BEGIN IMMEDIATE")
        output_connection.execute(
            "INSERT INTO messages (topic_id,timestamp,data) "
            f"SELECT {case_expression}, timestamp, data FROM source.messages"
        )
        output_connection.commit()
        output_connection.execute("DETACH DATABASE source")

    for statement in index_sql:
        output_connection.execute(statement)
    output_connection.commit()

    output_counts = {
        str(name): int(count)
        for name, count in output_connection.execute(
            "SELECT topics.name, COUNT(messages.id) FROM topics "
            "LEFT JOIN messages ON messages.topic_id=topics.id GROUP BY topics.id"
        )
    }
    first_timestamp_ns, last_timestamp_ns = output_connection.execute(
        "SELECT MIN(timestamp), MAX(timestamp) FROM messages"
    ).fetchone()
    output_connection.close()
    expected_counts = {name: count for counts in source_counts for name, count in counts.items()}
    if output_counts != expected_counts:
        raise SystemExit(f"merged counts differ: expected={expected_counts}, actual={output_counts}")

    subprocess.run(["ros2", "bag", "reindex", str(working_output)], check=True)
    metadata = working_output / "metadata.yaml"
    if not metadata.is_file():
        raise SystemExit("ros2 bag reindex did not create metadata.yaml")
    repair_metadata(working_output)

    contract = {
        "created_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "method": "SQLite bulk CDR copy; ros2 bag reindex metadata; reader orders by timestamp index",
        "storage_id": "sqlite3",
        "sources": [
            {
                "path": str(path),
                "metadata_sha256": sha256(path / "metadata.yaml"),
                "payload_bytes": payload.stat().st_size,
                "message_counts": counts,
            }
            for path, payload, counts in zip(inputs, payloads, source_counts)
        ],
        "output": str(output),
        "output_counts": output_counts,
        "total_messages": sum(output_counts.values()),
        "first_record_timestamp_ns": int(first_timestamp_ns),
        "last_record_timestamp_ns": int(last_timestamp_ns),
        "duration_s": (int(last_timestamp_ns) - int(first_timestamp_ns)) * 1e-9,
        "metadata_sha256": sha256(metadata),
    }
    (working_output / "merge_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if working_output != output:
        shutil.copytree(working_output, output)
        shutil.rmtree(working_output)
    print(json.dumps({
        "output": str(output),
        "topics": len(output_counts),
        "messages": contract["total_messages"],
        "duration_s": contract["duration_s"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
