#!/usr/bin/env python3
"""Run the formal M3DGR legacy/new/Voxel-SLAM backend matrix."""

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
METHODS = ("lightning_legacy", "lightning_new_backend", "voxel_slam_full")
SEED = 20260718


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    code = repo.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "m3dgr_backend",
    )
    parser.add_argument(
        "--inventory", type=Path,
        default=code / "_m3dgr_work" / "bench" / "inventory" / "bag_inventory.json",
    )
    parser.add_argument("--ros1-root", type=Path, default=code / "_m3dgr_work" / "ros1_normalized")
    parser.add_argument("--ros2-root", type=Path, default=code / "_m3dgr_work" / "ros2_bags")
    parser.add_argument("--sequences", default=",".join(SEQUENCES))
    parser.add_argument("--methods", default=",".join(METHODS))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--play-rate", type=float, default=1.0)
    parser.add_argument("--cpu-set", default="0-7")
    parser.add_argument("--cpu-count", type=int, default=8)
    parser.add_argument("--ros-port-base", type=int, default=11620)
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


def split_values(value: str, allowed: tuple[str, ...], label: str) -> tuple[str, ...]:
    selected = tuple(item.strip() for item in value.split(",") if item.strip())
    unknown = sorted(set(selected) - set(allowed))
    if not selected or unknown:
        raise SystemExit(f"invalid {label}: selected={selected}, unknown={unknown}")
    return selected


def parse_metadata(path: Path) -> dict[str, str]:
    values = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    return values


def metadata_passed(path: Path, sequence: str, method: str, repeat: int, config: Path | None) -> bool:
    values = parse_metadata(path)
    base = (
        values.get("sequence") == sequence
        and values.get("repeat") == str(repeat)
        and values.get("completion") == "reached_final_lidar"
        and values.get("invalid_count") == "0"
        and values.get("nonmonotonic_count") == "0"
        and values.get("excessive_output_gap_count") == "0"
    )
    if method == "voxel_slam_full":
        return (
            base
            and values.get("method") == "voxel_slam_full_backend"
            and int(values.get("pcd_count", "0")) == int(values.get("trajectory_lines", "-1"))
        )
    return (
        base
        and values.get("method") == "lightning_lm_offline_slam_map_export"
        and values.get("algorithm_rc") == "0"
        and values.get("watchdog_status") == "completed"
        and int(values.get("map_chunk_count", "0")) > 0
        and config is not None
        and values.get("config_sha256") == sha256(config)
    )


def generate_configs(repo: Path, state_dir: Path, dry_run: bool) -> dict[str, Path]:
    config_dir = state_dir / "configs"
    base = repo / "config" / "reproduction" / "single_lidar" / "m3dgr" / "lightning_m3dgr_mid360_benchmark.yaml"
    configs = {
        "lightning_legacy": config_dir / "lightning_m3dgr_legacy.yaml",
        "lightning_new_backend": config_dir / "lightning_m3dgr_ba_btc_hba.yaml",
    }
    if dry_run:
        return configs
    generator = repo / "scripts" / "reproduction" / "single_lidar" / "m3dgr" / "generate_backend_ablation_configs.py"
    if not all(path.is_file() for path in configs.values()):
        subprocess.run(
            [sys.executable, str(generator), "--base", str(base), "--output-dir", str(config_dir),
             "--variants", "legacy", "ba_btc_hba"],
            check=True,
        )
    return configs


def command_for(
    args: argparse.Namespace,
    repo: Path,
    runners: dict[str, Path],
    configs: dict[str, Path],
    sequence: str,
    method: str,
    repeat: int,
    run_dir: Path,
    ros_port: int,
) -> list[str]:
    if method.startswith("lightning_"):
        return [
            "wsl", "-d", "Ubuntu-22.04", "--", "bash", to_wsl(runners[method]),
            "--bag", to_wsl(args.ros2_root / sequence),
            "--config", to_wsl(configs[method]),
            "--output-dir", to_wsl(run_dir),
            "--sequence", sequence, "--repeat", str(repeat),
            "--playback-rate", str(args.play_rate),
            "--cpu-set", args.cpu_set, "--cpu-count", str(args.cpu_count),
            "--wait-ui", "false",
        ]
    environment = [
        f"BENCH_PLAY_RATE={args.play_rate}",
        f"BENCH_CPUSET={args.cpu_set}",
        f"BENCH_CPU_COUNT={args.cpu_count}",
        f"BENCH_INVENTORY_JSON={to_wsl(args.inventory)}",
        f"BENCH_ROS_PORT={ros_port}",
    ]
    return [
        "wsl", "-d", "Ubuntu-20.04", "--", "env", *environment,
        "bash", to_wsl(runners[method]),
        to_wsl(args.ros1_root / sequence / f"{sequence}.bag"), sequence, to_wsl(run_dir), str(repeat),
    ]


def main() -> int:
    args = parse_args()
    repo = Path(__file__).resolve().parents[3]
    code = repo.parents[1]
    sequences = split_values(args.sequences, SEQUENCES, "sequences")
    methods = split_values(args.methods, METHODS, "methods")
    if args.repeats < 1 or args.play_rate <= 0 or args.cpu_count < 1:
        raise SystemExit("repeats, play-rate and cpu-count must be positive")

    state_dir = args.output_root / "_state"
    configs = generate_configs(repo, state_dir, args.dry_run)
    slam_runner = repo / "scripts" / "run_slam_offline.sh"
    voxel_runner = repo / "scripts" / "reproduction" / "single_lidar" / "m3dgr" / "run_voxel_slam_full_backend.sh"
    runners = {
        "lightning_legacy": slam_runner,
        "lightning_new_backend": slam_runner,
        "voxel_slam_full": voxel_runner,
    }
    required = [args.inventory, *[runners[method] for method in methods]]
    if not args.dry_run:
        required.extend(configs[method] for method in methods if method in configs)
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
        "sequences": list(sequences),
        "methods": list(methods),
        "repeats": args.repeats,
        "play_rate": args.play_rate,
        "cpu_set": args.cpu_set,
        "cpu_count": args.cpu_count,
        "inventory": {"path": str(args.inventory.resolve()), "sha256": sha256(args.inventory)},
        "runners": {method: {"path": str(runners[method].resolve()), "sha256": sha256(runners[method])}
                    for method in methods},
        "configs": {method: {"path": str(configs[method].resolve()), "sha256": sha256(configs[method])}
                    for method in methods if method in configs and configs[method].is_file()},
    }
    fingerprint = hashlib.sha256(
        json.dumps(identity, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    identity["experiment_fingerprint"] = fingerprint

    blocks = [(sequence, repeat) for sequence in sequences for repeat in range(1, args.repeats + 1)]
    rng = random.Random(SEED)
    rng.shuffle(blocks)
    schedule = []
    for sequence, repeat in blocks:
        order = list(methods)
        rng.shuffle(order)
        schedule.extend({"sequence": sequence, "method": method, "repeat": repeat} for method in order)

    if args.dry_run:
        for index, item in enumerate(schedule, 1):
            run_dir = args.output_root / item["sequence"] / item["method"] / f"repeat_{item['repeat']:02d}"
            print(subprocess.list2cmdline(command_for(
                args, repo, runners, configs, item["sequence"], item["method"], item["repeat"], run_dir,
                args.ros_port_base + index,
            )))
        return 0

    state_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = state_dir / "experiment_manifest.json"
    if manifest_path.exists():
        if json.loads(manifest_path.read_text(encoding="utf-8")) != identity:
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
        config = configs.get(method)
        returncode: int | None = None
        if metadata_passed(run_dir / "run_metadata.txt", sequence, method, repeat, config):
            status = "passed_existing"
        elif run_dir.exists() and any(run_dir.iterdir()):
            status = "failed_existing"
            failures += 1
        else:
            command = command_for(
                args, repo, runners, configs, sequence, method, repeat, run_dir, args.ros_port_base + index,
            )
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open("w", encoding="utf-8", errors="replace") as stream:
                stream.write(
                    f"started_at={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
                    f"experiment_fingerprint={fingerprint}\n"
                    f"command={subprocess.list2cmdline(command)}\n"
                )
                stream.flush()
                result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
            returncode = result.returncode
            status = "passed" if metadata_passed(
                run_dir / "run_metadata.txt", sequence, method, repeat, config
            ) else "failed"
            failures += int(status == "failed")
        progress.append({**item, "status": status, "returncode": returncode,
                         "run_dir": str(run_dir), "controller_log": str(log_path)})
        write_json(state_dir / "progress.json", {"items": progress, "failures": failures})
        print(f"[{index}/{len(schedule)}] {sequence} {method} repeat={repeat}: {status}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
