#!/usr/bin/env python3
"""Run formal SANY 20260701 single/four-LiDAR frontend and mapping trials."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import subprocess
import sys
from datetime import datetime
from pathlib import Path


VARIANTS = ("single114", "four_lidar", "four_no_noise", "drop114", "drop127", "drop187", "drop195")
PIPELINES = ("frontend", "slam")
SLAM_VARIANTS = ("single114", "four_lidar")
SEED = 20260720
DEFAULT_COMPLETION_TOLERANCE_S = 0.25
DROP114_COMPLETION_TOLERANCE_S = 2.0


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bag", type=Path,
        default=Path(r"F:\datasets\SANY\mid360\data_20260701\rosbag2_2026_07_01-17_30_42_merged"),
    )
    parser.add_argument(
        "--output-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "sany_20260701_mapping",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--play-rate", type=float, default=1.0)
    parser.add_argument("--cpu-set", default="0-7")
    parser.add_argument("--cpu-count", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def to_wsl(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        raise ValueError(f"expected a Windows drive path: {resolved}")
    return f"/mnt/{drive}/{resolved.as_posix()[2:].lstrip('/')}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def parse_metadata(path: Path) -> dict[str, str]:
    values = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    return values


def passed(path: Path, pipeline: str, variant: str, repeat: int, config: Path) -> bool:
    values = parse_metadata(path)
    base = (
        values.get("sequence") == f"SANY_20260701_{variant}"
        and values.get("repeat") == str(repeat)
        and values.get("completion") == "reached_final_lidar"
        and values.get("algorithm_rc") == "0"
        and values.get("watchdog_status") == "completed"
        and values.get("invalid_count") == "0"
        and values.get("nonmonotonic_count") == "0"
        and values.get("excessive_output_gap_count") == "0"
        and values.get("config_sha256") == sha256(config)
    )
    if pipeline == "frontend":
        return base and values.get("method") == "lightning_lm"
    return (
        base
        and values.get("method") == "lightning_lm_offline_slam_map_export"
        and int(values.get("map_chunk_count", "0")) > 0
    )


def main() -> int:
    args = parse_args()
    repo = Path(__file__).resolve().parents[3]
    state_dir = args.output_root / "_state"
    config_dir = state_dir / "configs"
    generator = repo / "scripts" / "reproduction" / "formal_report" / "prepare_sany_mapping_configs.py"
    configs = {variant: config_dir / f"sany_20260701_{variant}.yaml" for variant in VARIANTS}
    if not args.dry_run and not all(path.is_file() for path in configs.values()):
        subprocess.run([sys.executable, str(generator), "--output-dir", str(config_dir)], check=True)
    runners = {
        "frontend": repo / "scripts" / "run_frontend_offline.sh",
        "slam": repo / "scripts" / "run_slam_offline.sh",
    }
    required = [args.bag / "metadata.yaml", *runners.values()]
    if not args.dry_run:
        required.extend(configs.values())
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing inputs:\n" + "\n".join(missing))
    if args.repeats < 1 or args.play_rate <= 0 or args.cpu_count < 1:
        raise SystemExit("repeats, play-rate and cpu-count must be positive")

    identity = {
        "schema_version": 1,
        "seed": SEED,
        "variants": VARIANTS,
        "pipelines": PIPELINES,
        "cells": [f"{variant}/{pipeline}" for variant, pipeline in
                  ([(name, "frontend") for name in VARIANTS] + [(name, "slam") for name in SLAM_VARIANTS])],
        "repeats": args.repeats,
        "play_rate": args.play_rate,
        "cpu_set": args.cpu_set,
        "cpu_count": args.cpu_count,
        "completion_tolerances_s": {
            "default": DEFAULT_COMPLETION_TOLERANCE_S,
            "drop114": DROP114_COMPLETION_TOLERANCE_S,
        },
        "bag": {"path": str(args.bag.resolve()), "metadata_sha256": sha256(args.bag / "metadata.yaml")},
        "runners": {name: {"path": str(path.resolve()), "sha256": sha256(path)} for name, path in runners.items()},
        "configs": {name: {"path": str(path.resolve()), "sha256": sha256(path)}
                    for name, path in configs.items() if path.is_file()},
    }
    identity["experiment_fingerprint"] = hashlib.sha256(
        json.dumps(identity, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    cells = [(variant, "frontend") for variant in VARIANTS] + [(variant, "slam") for variant in SLAM_VARIANTS]
    schedule = []
    for repeat in range(1, args.repeats + 1):
        block = [{"variant": variant, "pipeline": pipeline, "repeat": repeat} for variant, pipeline in cells]
        random.Random(SEED + repeat).shuffle(block)
        schedule.extend(block)

    def command_for(item: dict[str, object], run_dir: Path) -> list[str]:
        variant = str(item["variant"])
        pipeline = str(item["pipeline"])
        command = [
            "wsl", "-d", "Ubuntu-22.04", "--", "bash", to_wsl(runners[pipeline]),
            "--bag", to_wsl(args.bag), "--config", to_wsl(configs[variant]),
            "--output-dir", to_wsl(run_dir),
            "--sequence", f"SANY_20260701_{variant}", "--repeat", str(item["repeat"]),
            "--playback-rate", str(args.play_rate),
            "--cpu-set", args.cpu_set, "--cpu-count", str(args.cpu_count),
            "--wait-ui", "false",
        ]
        if pipeline == "frontend":
            command.extend([
                "--output-lidar-tum", to_wsl(run_dir / "results" / "trajectory_lidar114.tum"),
                "--development-exposed", "true", "--confirmatory", "true",
            ])
            if variant == "drop114":
                # Lidar 114 is intentionally absent, so bag inspection cannot
                # use its final stamp and falls back to the later bag end.  The
                # observed valid trajectory tail is 1.828 s before that bound.
                command.extend(["--completion-tolerance", str(DROP114_COMPLETION_TOLERANCE_S)])
        return command

    if args.dry_run:
        for item in schedule:
            run_dir = args.output_root / str(item["variant"]) / str(item["pipeline"]) \
                / f"repeat_{int(item['repeat']):02d}"
            print(subprocess.list2cmdline(command_for(item, run_dir)))
        return 0

    state_dir.mkdir(parents=True, exist_ok=True)
    manifest = state_dir / "experiment_manifest.json"
    if manifest.exists():
        if json.loads(manifest.read_text(encoding="utf-8")) != identity:
            raise SystemExit(f"existing experiment identity differs: {manifest}")
    else:
        write_json(manifest, identity)
    write_json(state_dir / "schedule.json", schedule)

    progress = []
    failures = 0
    for index, item in enumerate(schedule, 1):
        variant = str(item["variant"])
        pipeline = str(item["pipeline"])
        repeat = int(item["repeat"])
        run_dir = args.output_root / variant / pipeline / f"repeat_{repeat:02d}"
        log = state_dir / "controller_logs" / variant / pipeline / f"repeat_{repeat:02d}.log"
        returncode: int | None = None
        if passed(run_dir / "run_metadata.txt", pipeline, variant, repeat, configs[variant]):
            status = "passed_existing"
        elif run_dir.exists() and any(run_dir.iterdir()):
            status = "failed_existing"
            failures += 1
        else:
            command = command_for(item, run_dir)
            log.parent.mkdir(parents=True, exist_ok=True)
            with log.open("w", encoding="utf-8", errors="replace") as stream:
                stream.write(
                    f"started_at={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
                    f"experiment_fingerprint={identity['experiment_fingerprint']}\n"
                    f"command={subprocess.list2cmdline(command)}\n"
                )
                stream.flush()
                result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
            returncode = result.returncode
            status = "passed" if passed(
                run_dir / "run_metadata.txt", pipeline, variant, repeat, configs[variant]
            ) else "failed"
            failures += int(status == "failed")
        progress.append({**item, "status": status, "returncode": returncode,
                         "run_dir": str(run_dir), "controller_log": str(log)})
        write_json(state_dir / "progress.json", {"items": progress, "failures": failures})
        print(f"[{index}/{len(schedule)}] {variant} {pipeline} repeat={repeat}: {status}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
