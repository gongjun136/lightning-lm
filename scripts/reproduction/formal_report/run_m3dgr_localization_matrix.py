#!/usr/bin/env python3
"""Run M3DGR localization against Lightning and Voxel-SLAM maps."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path


SEQUENCES = ("Grass02", "Outdoor04", "Z-Rough-Road01", "Dark01")
MAP_VARIANTS = ("lightning_new_backend", "voxel_slam_full")
SEED = 20260719


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    code = repo.parents[1]
    formal = repo / "runs" / "formal_report_20260717"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend-root", type=Path, default=formal / "m3dgr_backend")
    parser.add_argument("--output-root", type=Path, default=formal / "m3dgr_localization")
    parser.add_argument("--ros2-root", type=Path, default=code / "_m3dgr_work" / "ros2_bags")
    parser.add_argument("--sequences", default=",".join(SEQUENCES))
    parser.add_argument("--map-variants", default=",".join(MAP_VARIANTS))
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


def passed(path: Path, sequence: str, repeat: int, map_index: Path, config: Path) -> bool:
    values = parse_metadata(path)
    return (
        values.get("method") == "lightning_lm_offline_localization"
        and values.get("sequence") == sequence
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


def source_map(args: argparse.Namespace, sequence: str, variant: str) -> Path:
    if variant == "lightning_new_backend":
        return args.backend_root / sequence / variant / "repeat_01" / "data" / "new_map"
    return args.output_root / "maps" / sequence / variant / "data" / "new_map"


def prepare_inputs(
    args: argparse.Namespace,
    repo: Path,
    sequences: tuple[str, ...],
    variants: tuple[str, ...],
) -> tuple[dict[tuple[str, str], Path], dict[tuple[str, str], Path], Path]:
    state_dir = args.output_root / "_state"
    # Use the real ELF build output here: Windows cannot stat the Linux symlink in install/.
    converter = repo / "bin" / "convert_voxel_slam_map"
    config_generator = repo / "scripts" / "reproduction" / "formal_report" / "prepare_m3dgr_localization_config.py"
    base_config = repo / "config" / "reproduction" / "single_lidar" / "m3dgr" / "lightning_m3dgr_mid360_benchmark.yaml"
    if args.dry_run:
        maps = {(sequence, variant): source_map(args, sequence, variant)
                for sequence in sequences for variant in variants}
        configs = {(sequence, variant): state_dir / "configs" / sequence / f"{variant}.yaml"
                   for sequence in sequences for variant in variants}
        return maps, configs, converter

    if not converter.is_file():
        raise SystemExit(f"map converter is not built: {converter}")
    maps: dict[tuple[str, str], Path] = {}
    configs: dict[tuple[str, str], Path] = {}
    for sequence in sequences:
        for variant in variants:
            map_dir = source_map(args, sequence, variant)
            if variant == "voxel_slam_full" and not (map_dir / "index.txt").is_file():
                voxel_source = args.backend_root / sequence / variant / "repeat_01" / "voxel_output" / sequence
                log = state_dir / "map_conversion_logs" / f"{sequence}.log"
                log.parent.mkdir(parents=True, exist_ok=True)
                converter_args = [
                    to_wsl(converter),
                    f"--input_dir={to_wsl(voxel_source)}",
                    f"--output_map_dir={to_wsl(map_dir)}",
                    "--voxel_size=0.1", "--frame_stride=1",
                ]
                install_setup = to_wsl(repo / "install" / "setup.bash")
                command = [
                    "wsl", "-d", "Ubuntu-22.04", "--", "bash", "-lc",
                    "source /opt/ros/humble/setup.bash && "
                    f"source {shlex.quote(install_setup)} && exec "
                    + " ".join(shlex.quote(value) for value in converter_args),
                ]
                with log.open("w", encoding="utf-8", errors="replace") as stream:
                    stream.write(f"command={subprocess.list2cmdline(command)}\n")
                    stream.flush()
                    subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
            for required in (map_dir / "index.txt", map_dir / "global.pcd"):
                if not required.is_file():
                    raise SystemExit(f"missing map artifact: {required}")
            maps[(sequence, variant)] = map_dir
            config = state_dir / "configs" / sequence / f"{variant}.yaml"
            if not config.is_file():
                subprocess.run(
                    [sys.executable, str(config_generator), "--base", str(base_config),
                     "--map-index", str(map_dir / "index.txt"), "--output", str(config)],
                    check=True,
                )
            configs[(sequence, variant)] = config
    return maps, configs, converter


def main() -> int:
    args = parse_args()
    repo = Path(__file__).resolve().parents[3]
    sequences = split_values(args.sequences, SEQUENCES, "sequences")
    variants = split_values(args.map_variants, MAP_VARIANTS, "map variants")
    if args.repeats < 1 or args.play_rate <= 0 or args.cpu_count < 1:
        raise SystemExit("repeats, play-rate and cpu-count must be positive")
    runner = repo / "scripts" / "run_loc_offline.sh"
    localization_binary = repo / "bin" / "run_loc_offline"
    maps, configs, converter = prepare_inputs(args, repo, sequences, variants)

    required = [runner, localization_binary]
    if not args.dry_run:
        required.append(converter)
        required.extend(path / "index.txt" for path in maps.values())
        required.extend(configs.values())
    required.extend(args.ros2_root / sequence / "metadata.yaml" for sequence in sequences)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing inputs:\n" + "\n".join(missing))

    identity = {
        "schema_version": 1,
        "seed": SEED,
        "sequences": list(sequences),
        "map_variants": list(variants),
        "repeats": args.repeats,
        "play_rate": args.play_rate,
        "cpu_set": args.cpu_set,
        "cpu_count": args.cpu_count,
        "runner": {"path": str(runner.resolve()), "sha256": sha256(runner)},
        "localization_binary": {
            "path": str(localization_binary.resolve()),
            "sha256": sha256(localization_binary),
        },
        "maps": {
            f"{sequence}/{variant}": {
                "path": str(maps[(sequence, variant)].resolve()),
                "index_sha256": sha256(maps[(sequence, variant)] / "index.txt")
                if (maps[(sequence, variant)] / "index.txt").is_file() else None,
                "global_pcd_sha256": sha256(maps[(sequence, variant)] / "global.pcd")
                if (maps[(sequence, variant)] / "global.pcd").is_file() else None,
                "config_sha256": sha256(configs[(sequence, variant)])
                if configs[(sequence, variant)].is_file() else None,
            }
            for sequence in sequences for variant in variants
        },
    }
    identity["experiment_fingerprint"] = hashlib.sha256(
        json.dumps(identity, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    blocks = [(sequence, repeat) for sequence in sequences for repeat in range(1, args.repeats + 1)]
    rng = random.Random(SEED)
    rng.shuffle(blocks)
    schedule = []
    for sequence, repeat in blocks:
        order = list(variants)
        rng.shuffle(order)
        schedule.extend({"sequence": sequence, "map_variant": variant, "repeat": repeat} for variant in order)

    def make_command(item: dict[str, object], run_dir: Path) -> list[str]:
        sequence = str(item["sequence"])
        variant = str(item["map_variant"])
        return [
            "wsl", "-d", "Ubuntu-22.04", "--", "bash", to_wsl(runner),
            "--bag", to_wsl(args.ros2_root / sequence),
            "--config", to_wsl(configs[(sequence, variant)]),
            "--map", to_wsl(maps[(sequence, variant)]),
            "--output-dir", to_wsl(run_dir),
            "--sequence", sequence, "--repeat", str(item["repeat"]),
            "--playback-rate", str(args.play_rate),
            "--cpu-set", args.cpu_set, "--cpu-count", str(args.cpu_count),
            "--wait-ui", "false", "--publish-topics", "false", "--use-config-initial-pose", "true",
        ]

    if args.dry_run:
        for item in schedule:
            run_dir = args.output_root / str(item["sequence"]) / str(item["map_variant"]) \
                / f"repeat_{int(item['repeat']):02d}"
            print(subprocess.list2cmdline(make_command(item, run_dir)))
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
        sequence = str(item["sequence"])
        variant = str(item["map_variant"])
        repeat = int(item["repeat"])
        run_dir = args.output_root / sequence / variant / f"repeat_{repeat:02d}"
        log = state_dir / "controller_logs" / sequence / variant / f"repeat_{repeat:02d}.log"
        returncode: int | None = None
        if passed(run_dir / "run_metadata.txt", sequence, repeat, maps[(sequence, variant)] / "index.txt",
                  configs[(sequence, variant)]):
            status = "passed_existing"
        elif run_dir.exists() and any(run_dir.iterdir()):
            status = "failed_existing"
            failures += 1
        else:
            command = make_command(item, run_dir)
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
                run_dir / "run_metadata.txt", sequence, repeat, maps[(sequence, variant)] / "index.txt",
                configs[(sequence, variant)],
            ) else "failed"
            failures += int(status == "failed")
        progress.append({**item, "status": status, "returncode": returncode,
                         "run_dir": str(run_dir), "controller_log": str(log)})
        write_json(state_dir / "progress.json", {"items": progress, "failures": failures})
        print(f"[{index}/{len(schedule)}] {sequence} {variant} repeat={repeat}: {status}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
