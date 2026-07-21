#!/usr/bin/env python3
"""Run three-repeat SANY 20260716 startup and forced-loss relocalization."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import subprocess
from datetime import datetime
from pathlib import Path


DATASETS = ("data1", "data2")
MODES = ("startup", "forced_loss")
FORCED_FRAMES = {"data1": 600, "data2": 350}
SEED = 20260721


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    formal = repo / "runs" / "formal_report_20260717"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root", type=Path,
        default=Path(r"F:\datasets\SANY\mid360\data_20260716"),
    )
    parser.add_argument("--mapping-root", type=Path, default=formal / "sany_20260701_mapping")
    parser.add_argument(
        "--reference-root", type=Path, default=formal / "sany_voxel114_reference",
    )
    parser.add_argument("--output-root", type=Path, default=formal / "sany_20260716_relocalization")
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


def passed(path: Path, dataset: str, repeat: int, map_index: Path, config: Path) -> bool:
    values = parse_metadata(path)
    return (
        values.get("method") == "lightning_lm_offline_localization"
        and values.get("sequence") == f"SANY_20260716_{dataset}"
        and values.get("repeat") == str(repeat)
        and values.get("completion") == "reached_final_lidar"
        and values.get("algorithm_rc") == "0"
        and values.get("watchdog_status") == "completed"
        and values.get("invalid_count") == "0"
        and values.get("nonmonotonic_count") == "0"
        and values.get("excessive_output_gap_count") == "0"
        and values.get("map_index_sha256") == sha256(map_index)
        and values.get("config_sha256") == sha256(config)
    )


def main() -> int:
    args = parse_args()
    repo = Path(__file__).resolve().parents[3]
    runner = repo / "scripts" / "run_loc_offline.sh"
    localization_binary = repo / "bin" / "run_loc_offline"
    map_dir = args.mapping_root / "four_lidar" / "slam" / "repeat_01" / "data" / "new_map"
    config = args.mapping_root / "_state" / "configs" / "sany_20260701_four_lidar.yaml"
    bags = {
        dataset: args.data_root / dataset / "codex_relocalization_eval" / "merged_bag_sqlite_v2"
        for dataset in DATASETS
    }
    truth = {
        dataset: args.reference_root / dataset / "results" / "trajectory_voxel_opt.tum"
        for dataset in DATASETS
    }
    required = [runner, localization_binary]
    if not args.dry_run:
        required.extend([map_dir / "index.txt", map_dir / "global.pcd", config])
    required.extend(bag / "metadata.yaml" for bag in bags.values())
    required.extend(truth.values())
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing inputs:\n" + "\n".join(missing))
    if args.repeats < 1 or args.play_rate <= 0 or args.cpu_count < 1:
        raise SystemExit("repeats, play-rate and cpu-count must be positive")

    identity = {
        "schema_version": 1,
        "seed": SEED,
        "datasets": DATASETS,
        "modes": MODES,
        "forced_frames": FORCED_FRAMES,
        "repeats": args.repeats,
        "play_rate": args.play_rate,
        "cpu_set": args.cpu_set,
        "cpu_count": args.cpu_count,
        "runner": {"path": str(runner.resolve()), "sha256": sha256(runner)},
        "localization_binary": {
            "path": str(localization_binary.resolve()),
            "sha256": sha256(localization_binary),
        },
        "map": {
            "path": str(map_dir.resolve()),
            "index_sha256": sha256(map_dir / "index.txt") if (map_dir / "index.txt").is_file() else None,
            "global_pcd_sha256": sha256(map_dir / "global.pcd") if (map_dir / "global.pcd").is_file() else None,
        },
        "config": {"path": str(config.resolve()), "sha256": sha256(config) if config.is_file() else None},
        "bags": {dataset: {"path": str(bags[dataset].resolve()),
                            "metadata_sha256": sha256(bags[dataset] / "metadata.yaml")}
                 for dataset in DATASETS},
        "proxy_truth": {dataset: {"path": str(truth[dataset].resolve()), "sha256": sha256(truth[dataset])}
                        for dataset in DATASETS},
    }
    identity["experiment_fingerprint"] = hashlib.sha256(
        json.dumps(identity, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    schedule = []
    for repeat in range(1, args.repeats + 1):
        block = [{"dataset": dataset, "mode": mode, "repeat": repeat}
                 for dataset in DATASETS for mode in MODES]
        random.Random(SEED + repeat).shuffle(block)
        schedule.extend(block)

    def command_for(item: dict[str, object], run_dir: Path) -> list[str]:
        dataset = str(item["dataset"])
        mode = str(item["mode"])
        command = [
            "wsl", "-d", "Ubuntu-22.04", "--", "bash", to_wsl(runner),
            "--bag", to_wsl(bags[dataset]), "--config", to_wsl(config), "--map", to_wsl(map_dir),
            "--output-dir", to_wsl(run_dir),
            "--sequence", f"SANY_20260716_{dataset}", "--repeat", str(item["repeat"]),
            "--playback-rate", str(args.play_rate),
            "--cpu-set", args.cpu_set, "--cpu-count", str(args.cpu_count),
            "--wait-ui", "false", "--publish-topics", "false", "--use-config-initial-pose", "false",
        ]
        if mode == "forced_loss":
            command.extend(["--", f"--force_relocalization_frame={FORCED_FRAMES[dataset]}"])
        return command

    if args.dry_run:
        for item in schedule:
            run_dir = args.output_root / str(item["dataset"]) / str(item["mode"]) \
                / f"repeat_{int(item['repeat']):02d}"
            print(subprocess.list2cmdline(command_for(item, run_dir)))
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

    progress = []
    failures = 0
    for index, item in enumerate(schedule, 1):
        dataset = str(item["dataset"])
        mode = str(item["mode"])
        repeat = int(item["repeat"])
        run_dir = args.output_root / dataset / mode / f"repeat_{repeat:02d}"
        log = state_dir / "controller_logs" / dataset / mode / f"repeat_{repeat:02d}.log"
        returncode: int | None = None
        if passed(run_dir / "run_metadata.txt", dataset, repeat, map_dir / "index.txt", config):
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
                run_dir / "run_metadata.txt", dataset, repeat, map_dir / "index.txt", config
            ) else "failed"
            failures += int(status == "failed")
        progress.append({**item, "status": status, "returncode": returncode,
                         "run_dir": str(run_dir), "controller_log": str(log)})
        write_json(state_dir / "progress.json", {"items": progress, "failures": failures})
        print(f"[{index}/{len(schedule)}] {dataset} {mode} repeat={repeat}: {status}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
