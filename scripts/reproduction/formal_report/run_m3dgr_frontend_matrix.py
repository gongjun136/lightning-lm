#!/usr/bin/env python3
"""Run the formal four-sequence M3DGR frontend comparison matrix."""

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


SEQUENCES = ("Grass02", "Outdoor04", "Z-Rough-Road01", "Dark01")
METHODS = ("lightning_lm", "fastlio", "fastlivo2_lio", "voxel_slam_frontend")
SEED = 20260717


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    code = repo.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "m3dgr_frontend",
    )
    parser.add_argument(
        "--inventory", type=Path,
        default=code / "_m3dgr_work" / "bench" / "inventory" / "bag_inventory.json",
    )
    parser.add_argument(
        "--ros1-root", type=Path,
        default=code / "_m3dgr_work" / "ros1_normalized",
    )
    parser.add_argument(
        "--ros2-root", type=Path,
        default=code / "_m3dgr_work" / "ros2_bags",
    )
    parser.add_argument("--sequences", default=",".join(SEQUENCES))
    parser.add_argument("--methods", default=",".join(METHODS))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--play-rate", type=float, default=1.0)
    parser.add_argument("--cpu-set", default="0-7")
    parser.add_argument("--cpu-count", type=int, default=8)
    parser.add_argument("--ros-port-base", type=int, default=11420)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def to_wsl(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        raise ValueError(f"expected a Windows drive path: {resolved}")
    suffix = resolved.as_posix()[2:].lstrip("/")
    return f"/mnt/{drive}/{suffix}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_values(value: str, allowed: tuple[str, ...], label: str) -> tuple[str, ...]:
    selected = tuple(item.strip() for item in value.split(",") if item.strip())
    unknown = sorted(set(selected) - set(allowed))
    if not selected or unknown:
        raise SystemExit(f"invalid {label}: selected={selected}, unknown={unknown}")
    return selected


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def metadata_passed(path: Path, sequence: str, method: str, repeat: int) -> bool:
    if not path.is_file():
        return False
    values = {}
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    base = (
        values.get("sequence") == sequence
        and values.get("method") == method
        and values.get("repeat") == str(repeat)
        and values.get("completion") == "reached_final_lidar"
        and values.get("invalid_count") == "0"
        and values.get("nonmonotonic_count") == "0"
    )
    if method == "lightning_lm":
        return base and values.get("algorithm_rc") == "0" and values.get("watchdog_status") == "completed"
    return (
        base
        and values.get("rosbag_play_rc") == "0"
        and values.get("algorithm_contract") == "passed"
        and values.get("subscriber_contract") == "passed"
        and values.get("algorithm_launch_rc") == "0"
    )


def method_paths(code: Path, repo: Path) -> dict[str, Path]:
    return {
        "lightning_lm": repo / "scripts" / "run_frontend_offline.sh",
        "fastlio": code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_fastlio" / "src" / "FAST_LIO" / "reproduction" / "m3dgr" / "scripts" / "run_frontend_offline.sh",
        "fastlivo2_lio": code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_fastlivo2" / "src" / "FAST-LIVO2" / "reproduction" / "m3dgr" / "scripts" / "run_frontend_offline.sh",
        "voxel_slam_frontend": code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_voxel_slam" / "src" / "Voxel-SLAM" / "reproduction" / "m3dgr" / "scripts" / "run_frontend_offline.sh",
    }


def command_for(
    args: argparse.Namespace,
    code: Path,
    repo: Path,
    runners: dict[str, Path],
    sequence: str,
    method: str,
    repeat: int,
    run_dir: Path,
    fingerprint: str,
    ros_port: int,
) -> list[str]:
    if method == "lightning_lm":
        bag = args.ros2_root / sequence
        config = repo / "config" / "reproduction" / "single_lidar" / "m3dgr" / "lightning_m3dgr_mid360_benchmark.yaml"
        return [
            "wsl", "-d", "Ubuntu-22.04", "--", "bash", to_wsl(runners[method]),
            "--bag", to_wsl(bag), "--config", to_wsl(config),
            "--output-dir", to_wsl(run_dir),
            "--output-lidar-tum", to_wsl(run_dir / "trajectory_mid360.tum"),
            "--sequence", sequence, "--repeat", str(repeat),
            "--inventory-json", to_wsl(args.inventory),
            "--playback-rate", str(args.play_rate),
            "--cpu-set", args.cpu_set, "--cpu-count", str(args.cpu_count),
            "--experiment-fingerprint", fingerprint,
            "--development-exposed", "true", "--confirmatory", "true",
            "--wait-ui", "false",
        ]
    bag = args.ros1_root / sequence / f"{sequence}.bag"
    environment = [
        f"BENCH_PLAY_RATE={args.play_rate}",
        f"BENCH_CPUSET={args.cpu_set}",
        f"BENCH_CPU_COUNT={args.cpu_count}",
        f"BENCH_INVENTORY_JSON={to_wsl(args.inventory)}",
        f"BENCH_EXPERIMENT_FINGERPRINT={fingerprint}",
        "BENCH_DEVELOPMENT_EXPOSED=false",
        "BENCH_CONFIRMATORY=true",
        f"BENCH_ROS_PORT={ros_port}",
    ]
    return [
        "wsl", "-d", "Ubuntu-20.04", "--", "/usr/bin/env", *environment,
        "bash", to_wsl(runners[method]), to_wsl(bag), sequence, str(repeat), to_wsl(run_dir),
    ]


def main() -> int:
    args = parse_args()
    repo = Path(__file__).resolve().parents[3]
    code = repo.parents[1]
    sequences = split_values(args.sequences, SEQUENCES, "sequences")
    methods = split_values(args.methods, METHODS, "methods")
    if args.repeats < 1 or args.play_rate <= 0 or args.cpu_count < 1:
        raise SystemExit("repeats, play-rate and cpu-count must be positive")
    runners = method_paths(code, repo)
    required = [args.inventory, *runners.values()]
    for sequence in sequences:
        required.extend([
            args.ros1_root / sequence / f"{sequence}.bag",
            args.ros2_root / sequence / "metadata.yaml",
        ])
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing inputs:\n" + "\n".join(missing))

    identity = {
        "schema_version": 1,
        "seed": SEED,
        "sequences": sequences,
        "methods": methods,
        "repeats": args.repeats,
        "play_rate": args.play_rate,
        "cpu_set": args.cpu_set,
        "cpu_count": args.cpu_count,
        "inventory": {"path": str(args.inventory.resolve()), "sha256": sha256(args.inventory)},
        "runners": {method: {"path": str(runners[method].resolve()), "sha256": sha256(runners[method])} for method in methods},
    }
    canonical = json.dumps(identity, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    fingerprint = hashlib.sha256(canonical).hexdigest()
    identity["experiment_fingerprint"] = fingerprint

    blocks = [(sequence, repeat) for sequence in sequences for repeat in range(1, args.repeats + 1)]
    rng = random.Random(SEED)
    rng.shuffle(blocks)
    schedule = []
    for sequence, repeat in blocks:
        order = list(methods)
        rng.shuffle(order)
        for method in order:
            schedule.append({"sequence": sequence, "method": method, "repeat": repeat})

    if args.dry_run:
        for index, item in enumerate(schedule, 1):
            run_dir = args.output_root / item["sequence"] / item["method"] / f"repeat_{item['repeat']:02d}"
            command = command_for(
                args, code, repo, runners, item["sequence"], item["method"], item["repeat"],
                run_dir, fingerprint, args.ros_port_base + index,
            )
            print(subprocess.list2cmdline(command))
        return 0

    args.output_root.mkdir(parents=True, exist_ok=True)
    state_dir = args.output_root / "_state"
    state_dir.mkdir(exist_ok=True)
    manifest_path = state_dir / "experiment_manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        if previous != identity:
            raise SystemExit(f"existing experiment identity differs: {manifest_path}")
    else:
        write_json(manifest_path, identity)
    write_json(state_dir / "schedule.json", schedule)

    progress = []
    failures = 0
    for index, item in enumerate(schedule, 1):
        sequence, method, repeat = item["sequence"], item["method"], item["repeat"]
        run_dir = args.output_root / sequence / method / f"repeat_{repeat:02d}"
        log_path = state_dir / "controller_logs" / sequence / method / f"repeat_{repeat:02d}.log"
        status = "pending"
        returncode: int | None = None
        if metadata_passed(run_dir / "run_metadata.txt", sequence, method, repeat):
            status = "passed_existing"
        elif run_dir.exists() and any(run_dir.iterdir()):
            status = "failed_existing"
            failures += 1
        else:
            command = command_for(
                args, code, repo, runners, sequence, method, repeat, run_dir,
                fingerprint, args.ros_port_base + index,
            )
            log_path.parent.mkdir(parents=True, exist_ok=True)
            started = datetime.now().astimezone().isoformat(timespec="seconds")
            with log_path.open("w", encoding="utf-8", errors="replace") as stream:
                stream.write(f"started_at={started}\ncommand={subprocess.list2cmdline(command)}\n")
                stream.flush()
                result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
            returncode = result.returncode
            status = "passed" if metadata_passed(run_dir / "run_metadata.txt", sequence, method, repeat) else "failed"
            if status == "failed":
                failures += 1
        progress.append({**item, "status": status, "returncode": returncode, "run_dir": str(run_dir), "controller_log": str(log_path)})
        write_json(state_dir / "progress.json", {"items": progress, "failures": failures})
        print(f"[{index}/{len(schedule)}] {sequence} {method} repeat={repeat}: {status}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
