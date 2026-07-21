#!/usr/bin/env python3
"""Run the controlled Lightning-LM M3DGR backend ablation matrix."""

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

import yaml


SEQUENCES = ("Grass02", "Outdoor04", "Z-Rough-Road01", "Dark01")
METHOD_VARIANTS = {
    "frontend_only": "frontend_only",
    "legacy_controlled": "legacy_controlled",
    "local_ba_only": "local_ba_only",
    "loop_pgo": "loop_pgo",
    "full_backend": "ba_btc_hba",
}
SEED = 20260721


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[3]
    code = repo.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "m3dgr_backend_controlled_20260721",
    )
    parser.add_argument(
        "--inventory", type=Path,
        default=code / "_m3dgr_work" / "bench" / "inventory" / "bag_inventory.json",
    )
    parser.add_argument("--ros2-root", type=Path, default=code / "_m3dgr_work" / "ros2_bags")
    parser.add_argument("--sequences", default=",".join(SEQUENCES))
    parser.add_argument("--methods", default=",".join(METHOD_VARIANTS))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--play-rate", type=float, default=0.0)
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


def select_values(raw: str, allowed: tuple[str, ...], label: str) -> tuple[str, ...]:
    values = tuple(item.strip() for item in raw.split(",") if item.strip())
    unknown = sorted(set(values) - set(allowed))
    if not values or unknown:
        raise SystemExit(f"invalid {label}: selected={values}, unknown={unknown}")
    return values


def parse_metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    return values


def run_passed(run_dir: Path, sequence: str, repeat: int, config: Path) -> bool:
    values = parse_metadata(run_dir / "run_metadata.txt")
    return (
        values.get("sequence") == sequence
        and values.get("repeat") == str(repeat)
        and values.get("method") == "lightning_lm_offline_slam_map_export"
        and values.get("completion") == "reached_final_lidar"
        and values.get("algorithm_rc") == "0"
        and values.get("watchdog_status") == "completed"
        and values.get("invalid_count") == "0"
        and values.get("nonmonotonic_count") == "0"
        and values.get("excessive_output_gap_count") == "0"
        and values.get("backend_evaluation_only") == "true"
        and values.get("config_sha256") == sha256(config)
        and (run_dir / "results" / "trajectory_slam_keyframes_lio.tum").is_file()
    )


def frontend_contract(config: Path) -> dict[str, object]:
    payload = yaml.safe_load(config.read_text(encoding="utf-8"))
    payload.pop("backend", None)
    payload.get("common", {}).pop("benchmark_backend_variant", None)
    payload.get("system", {}).pop("with_loop_closing", None)
    return payload


def generate_configs(repo: Path, state_dir: Path, methods: tuple[str, ...]) -> dict[str, Path]:
    config_dir = state_dir / "configs"
    base = repo / "config" / "reproduction" / "single_lidar" / "m3dgr" / "lightning_m3dgr_mid360_benchmark.yaml"
    generator = repo / "scripts" / "reproduction" / "single_lidar" / "m3dgr" / "generate_backend_ablation_configs.py"
    variants = [METHOD_VARIANTS[method] for method in methods]
    subprocess.run(
        [sys.executable, str(generator), "--base", str(base), "--output-dir", str(config_dir),
         "--variants", *variants],
        check=True,
    )
    configs = {
        method: config_dir / f"lightning_m3dgr_{METHOD_VARIANTS[method]}.yaml"
        for method in methods
    }
    contracts = {method: frontend_contract(path) for method, path in configs.items()}
    reference_method = methods[0]
    mismatched = [method for method, value in contracts.items() if value != contracts[reference_method]]
    if mismatched:
        raise SystemExit(f"frontend configuration mismatch: reference={reference_method}, mismatched={mismatched}")
    return configs


def main() -> int:
    args = parse_args()
    repo = Path(__file__).resolve().parents[3]
    sequences = select_values(args.sequences, SEQUENCES, "sequences")
    methods = select_values(args.methods, tuple(METHOD_VARIANTS), "methods")
    if args.repeats < 1 or args.play_rate < 0 or args.cpu_count < 1:
        raise SystemExit("repeats/cpu-count must be positive and play-rate non-negative")

    state_dir = args.output_root / "_state"
    configs = generate_configs(repo, state_dir, methods)
    runner = repo / "scripts" / "run_slam_offline.sh"
    # The install entry is a WSL symlink that Windows pathlib cannot stat.
    # Record the real build output to fingerprint the exact executable.
    binary = repo / "bin" / "run_slam_offline"
    required = [args.inventory, runner, binary]
    required.extend(args.ros2_root / sequence / "metadata.yaml" for sequence in sequences)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing inputs:\n" + "\n".join(missing))

    identity = {
        "schema_version": 1,
        "seed": SEED,
        "design": "randomized complete blocks by sequence and repeat",
        "sequences": list(sequences),
        "methods": list(methods),
        "repeats": args.repeats,
        "play_rate": args.play_rate,
        "cpu_set": args.cpu_set,
        "cpu_count": args.cpu_count,
        "frontend_control": "all non-backend YAML fields identical; post-run raw LIO SHA-256 equality required",
        "backend_evaluation_only": True,
        "inventory": {"path": str(args.inventory.resolve()), "sha256": sha256(args.inventory)},
        "runner": {"path": str(runner.resolve()), "sha256": sha256(runner)},
        "binary": {"path": str(binary.resolve()), "sha256": sha256(binary)},
        "configs": {
            method: {"path": str(config.resolve()), "sha256": sha256(config)}
            for method, config in configs.items()
        },
    }
    identity["experiment_fingerprint"] = hashlib.sha256(
        json.dumps(identity, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    blocks = [(sequence, repeat) for sequence in sequences for repeat in range(1, args.repeats + 1)]
    rng = random.Random(SEED)
    rng.shuffle(blocks)
    schedule: list[dict[str, object]] = []
    for sequence, repeat in blocks:
        order = list(methods)
        rng.shuffle(order)
        schedule.extend({"sequence": sequence, "method": method, "repeat": repeat} for method in order)

    if args.dry_run:
        print(json.dumps({"identity": identity, "schedule": schedule}, ensure_ascii=False, indent=2))
        return 0

    state_dir.mkdir(parents=True, exist_ok=True)
    manifest = state_dir / "experiment_manifest.json"
    if manifest.exists() and json.loads(manifest.read_text(encoding="utf-8")) != identity:
        raise SystemExit(f"existing experiment identity differs: {manifest}")
    write_json(manifest, identity)
    write_json(state_dir / "schedule.json", schedule)

    progress: list[dict[str, object]] = []
    failures = 0
    for index, item in enumerate(schedule, 1):
        sequence = str(item["sequence"])
        method = str(item["method"])
        repeat = int(item["repeat"])
        run_dir = args.output_root / sequence / method / f"repeat_{repeat:02d}"
        log_path = state_dir / "controller_logs" / sequence / method / f"repeat_{repeat:02d}.log"
        config = configs[method]
        returncode: int | None = None
        if run_passed(run_dir, sequence, repeat, config):
            status = "passed_existing"
        elif run_dir.exists() and any(run_dir.iterdir()):
            status = "failed_existing"
            failures += 1
        else:
            command = [
                "wsl", "-d", "Ubuntu-22.04", "--", "bash", to_wsl(runner),
                "--bag", to_wsl(args.ros2_root / sequence),
                "--config", to_wsl(config),
                "--output-dir", to_wsl(run_dir),
                "--sequence", sequence, "--repeat", str(repeat),
                "--playback-rate", str(args.play_rate),
                "--cpu-set", args.cpu_set, "--cpu-count", str(args.cpu_count),
                "--wait-ui", "false", "--backend-evaluation-only", "true",
            ]
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open("w", encoding="utf-8", errors="replace") as stream:
                stream.write(
                    f"started_at={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
                    f"experiment_fingerprint={identity['experiment_fingerprint']}\n"
                    f"command={subprocess.list2cmdline(command)}\n"
                )
                stream.flush()
                result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
            returncode = result.returncode
            status = "passed" if run_passed(run_dir, sequence, repeat, config) else "failed"
            failures += int(status == "failed")
        lio_path = run_dir / "results" / "trajectory_slam_keyframes_lio.tum"
        progress.append({
            **item, "status": status, "returncode": returncode,
            "run_dir": str(run_dir), "controller_log": str(log_path),
            "frontend_lio_sha256": sha256(lio_path) if lio_path.is_file() else None,
        })
        write_json(state_dir / "progress.json", {"items": progress, "failures": failures})
        print(f"[{index}/{len(schedule)}] {sequence} {method} repeat={repeat}: {status}", flush=True)

    frontend_hashes: dict[str, list[str]] = {}
    for sequence in sequences:
        hashes = sorted({
            str(row["frontend_lio_sha256"])
            for row in progress
            if row["sequence"] == sequence and row["frontend_lio_sha256"]
        })
        frontend_hashes[sequence] = hashes
        if len(hashes) != 1:
            failures += 1
    control = {
        "status": "passed" if not failures else "failed",
        "expected_hashes_per_sequence": 1,
        "frontend_lio_hashes": frontend_hashes,
        "failures": failures,
    }
    write_json(state_dir / "frontend_control_validation.json", control)
    print(json.dumps(control, ensure_ascii=False, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
