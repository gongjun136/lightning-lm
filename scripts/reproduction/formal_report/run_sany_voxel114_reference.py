#!/usr/bin/env python3
"""Generate traceable Voxel-SLAM 114 proxy references for formal SANY experiments."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime
from pathlib import Path


DATASETS = ("data_20260701", "data1", "data2")


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root", type=Path,
        default=Path(r"F:\datasets\SANY\mid360\data_20260716"),
    )
    parser.add_argument(
        "--mapping-bag", type=Path,
        default=Path(
            r"F:\datasets\SANY\mid360\data_20260701\voxel_slam_114_20260701"
            r"\sany_20260701_livox_114_ros1.bag"
        ),
    )
    parser.add_argument(
        "--output-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "sany_voxel114_reference",
    )
    parser.add_argument("--datasets", default=",".join(DATASETS))
    parser.add_argument("--play-rate", type=float, default=1.0)
    parser.add_argument("--cpu-set", default="0-7")
    parser.add_argument("--cpu-count", type=int, default=8)
    parser.add_argument("--ros-port-base", type=int, default=11920)
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
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    return values


def sequence_name(dataset: str) -> str:
    return "SANY_20260701_voxel114" if dataset == "data_20260701" \
        else f"SANY_20260716_{dataset}_voxel114"


def passed(path: Path, dataset: str, runner: Path, config: Path) -> bool:
    values = parse_metadata(path)
    try:
        return (
            values.get("method") == "voxel_slam_114_reference"
            and values.get("sequence") == sequence_name(dataset)
            and values.get("repeat") == "1"
            and values.get("runner_sha256") == sha256(runner)
            and values.get("config_file_sha256") == sha256(config)
            and values.get("completion") == "reached_final_lidar"
            and int(values.get("valid_lines", "0")) >= 3
            and int(values.get("invalid_lines", "-1")) == 0
            and int(values.get("nonmonotonic_lines", "-1")) == 0
            and 0.80 <= float(values.get("output_ratio", "0")) <= 1.01
        )
    except ValueError:
        return False


def main() -> int:
    args = parse_args()
    repo = Path(__file__).resolve().parents[3]
    code = repo.parents[1]
    runner = repo / "scripts" / "run_voxel_slam_114_reference.sh"
    config = code / "WSL_Ubuntu_20.04" / "ros1_ws" / "ws_voxel_slam" / "src" / "Voxel-SLAM" \
        / "VoxelSLAM" / "config" / "sany_20260701_livox_pc2_114.yaml"
    selected = tuple(item.strip() for item in args.datasets.split(",") if item.strip())
    unknown = sorted(set(selected) - set(DATASETS))
    if not selected or unknown:
        raise SystemExit(f"invalid datasets: selected={selected}, unknown={unknown}")
    if args.play_rate <= 0 or args.cpu_count < 1:
        raise SystemExit("play-rate and cpu-count must be positive")
    bags = {
        dataset: (
            args.mapping_bag if dataset == "data_20260701" else
            args.data_root / dataset / "codex_relocalization_eval" / "voxel_114_input.bag"
        )
        for dataset in selected
    }
    missing = [str(path) for path in (runner, config, *bags.values()) if not path.is_file()]
    if missing:
        raise SystemExit("missing inputs:\n" + "\n".join(missing))

    identity = {
        "schema_version": 1,
        "datasets": selected,
        "repeat_policy": "one deterministic proxy-reference generation per SANY bag",
        "play_rate": args.play_rate,
        "cpu_set": args.cpu_set,
        "cpu_count": args.cpu_count,
        "runner": {"path": str(runner.resolve()), "sha256": sha256(runner)},
        "config": {"path": str(config.resolve()), "sha256": sha256(config)},
        "bags": {
            dataset: {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for dataset, path in bags.items()
        },
    }
    identity["experiment_fingerprint"] = hashlib.sha256(
        json.dumps(identity, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    schedule = [{"dataset": dataset, "repeat": 1} for dataset in selected]

    def command_for(dataset: str, run_dir: Path, index: int) -> list[str]:
        return [
            "wsl", "-d", "Ubuntu-20.04", "--", "bash", to_wsl(runner),
            "--bag", to_wsl(bags[dataset]),
            "--output-dir", to_wsl(run_dir),
            "--sequence", sequence_name(dataset),
            "--ros-port", str(args.ros_port_base + index),
            "--cpu-set", args.cpu_set,
            "--cpu-count", str(args.cpu_count),
            "--play-rate", str(args.play_rate),
            "--repeat", "1",
        ]

    if args.dry_run:
        for index, item in enumerate(schedule, 1):
            run_dir = args.output_root / item["dataset"]
            print(subprocess.list2cmdline(command_for(item["dataset"], run_dir, index)))
        return 0

    state_dir = args.output_root / "_state"
    state_dir.mkdir(parents=True, exist_ok=True)
    manifest = state_dir / "experiment_manifest.json"
    if manifest.exists():
        if json.loads(manifest.read_text(encoding="utf-8")) != identity:
            raise SystemExit(f"existing experiment identity differs: {manifest}")
    else:
        write_json(manifest, identity)
    write_json(state_dir / "schedule.json", schedule)

    progress: list[dict[str, object]] = []
    failures = 0
    for index, item in enumerate(schedule, 1):
        dataset = item["dataset"]
        run_dir = args.output_root / dataset
        log = state_dir / "controller_logs" / f"{dataset}.log"
        returncode: int | None = None
        if passed(run_dir / "run_metadata.txt", dataset, runner, config):
            status = "passed_existing"
        elif run_dir.exists() and any(run_dir.iterdir()):
            status = "failed_existing"
            failures += 1
        else:
            command = command_for(dataset, run_dir, index)
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
            status = "passed" if passed(run_dir / "run_metadata.txt", dataset, runner, config) else "failed"
            failures += int(status == "failed")
        progress.append({
            **item, "status": status, "returncode": returncode,
            "run_dir": str(run_dir), "controller_log": str(log),
        })
        write_json(state_dir / "progress.json", {"items": progress, "failures": failures})
        print(f"[{index}/{len(schedule)}] {dataset}: {status}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
