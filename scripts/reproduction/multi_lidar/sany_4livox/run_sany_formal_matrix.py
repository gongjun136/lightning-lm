#!/usr/bin/env python3
"""Generate external SANY configs and run the frozen 7 x 3 offline matrix."""

from __future__ import annotations

import argparse
import atexit
import copy
import csv
import hashlib
import json
import math
import os
import random
import socket
import subprocess
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath
from typing import Any

import yaml


SEED = 20260711
VARIANTS = (
    ("C0_single114_noise", "single", None),
    ("C1_four_noise", "multi", None),
    ("C2_four_no_noise", "multi_no_noise", None),
    ("C3_drop114_lidar", "drop", 0),
    ("C4_drop127", "drop", 1),
    ("C5_drop187", "drop", 2),
    ("C6_drop195", "drop", 3),
)


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def to_wsl(path: Path | str) -> str:
    text = str(path)
    if text.startswith("/mnt/"):
        return text
    pure = PureWindowsPath(text)
    if not pure.drive:
        raise ValueError(f"drive-letter path required: {text}")
    return f"/mnt/{pure.drive[0].lower()}/" + "/".join(pure.parts[1:])


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def build_config(base: dict[str, Any], mode: str, dropped_id: int | None) -> dict[str, Any]:
    config = copy.deepcopy(base)
    if mode == "single":
        config["multi_lidar"]["enabled"] = False
        config["lidar_noise_model"]["enabled"] = True
    else:
        config["multi_lidar"]["enabled"] = True
        config["lidar_noise_model"]["enabled"] = mode != "multi_no_noise"
    if dropped_id is not None:
        config["multi_lidar"]["topics"][f"lidar_{dropped_id}"] = f"/fault_injection/missing_lidar_{dropped_id}"
    return config


def read_tum(path: Path) -> tuple[list[float], list[str]]:
    timestamps: list[float] = []
    errors: list[str] = []
    previous = -math.inf
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 8:
                errors.append(f"{path.name}:{line_number}: expected 8 columns")
                continue
            try:
                values = [float(value) for value in fields]
            except ValueError:
                errors.append(f"{path.name}:{line_number}: non-numeric value")
                continue
            if not all(math.isfinite(value) for value in values):
                errors.append(f"{path.name}:{line_number}: non-finite value")
                continue
            timestamp = values[0]
            if timestamp <= previous:
                errors.append(f"{path.name}:{line_number}: timestamp is not strictly increasing")
            previous = timestamp
            quaternion_norm = math.sqrt(sum(value * value for value in values[4:8]))
            if abs(quaternion_norm - 1.0) > 1e-3:
                errors.append(f"{path.name}:{line_number}: invalid quaternion norm")
            timestamps.append(timestamp)
    return timestamps, errors


def pcd_point_count(path: Path) -> int:
    with path.open("rb") as stream:
        for _ in range(64):
            line = stream.readline()
            if not line:
                break
            text = line.decode("ascii", errors="strict").strip()
            if text.startswith("POINTS "):
                return int(text.split()[1])
            if text.startswith("DATA "):
                break
    raise ValueError(f"PCD POINTS header missing: {path}")


def validate_outputs(run_dir: Path, variant: str) -> dict[str, Any]:
    metadata_path = run_dir / "controller_metadata.json"
    del metadata_path
    lidar_tum = run_dir / "results" / "trajectory_lidar114.tum"
    rear_tum = run_dir / "results" / "trajectory_rear_axle.tum"
    frame_stats = run_dir / "results" / "frame_stats.csv"
    map_pcd = run_dir / "results" / "map_lio.pcd"
    required = (lidar_tum, rear_tum, frame_stats, map_pcd)
    errors = [f"missing or empty output: {path.name}" for path in required if not path.is_file() or path.stat().st_size == 0]
    if errors:
        return {"passed": False, "errors": errors}

    lidar_times, lidar_errors = read_tum(lidar_tum)
    rear_times, rear_errors = read_tum(rear_tum)
    errors.extend(lidar_errors)
    errors.extend(rear_errors)
    if len(lidar_times) < 1000 or len(rear_times) < 1000:
        errors.append("trajectory has fewer than 1000 poses")
    if len(lidar_times) != len(rear_times) or any(abs(left - right) > 1e-6 for left, right in zip(lidar_times, rear_times)):
        errors.append("lidar/rear-axle trajectory timestamps differ")

    frame_rows: list[dict[str, str]] = []
    with frame_stats.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        expected_header = {
            "begin_time", "end_time", "partial", "merged_points", "present_lidar_ids", "missing_lidar_ids",
        }
        if reader.fieldnames is None or not expected_header.issubset(reader.fieldnames):
            errors.append("frame_stats header is incomplete")
        else:
            frame_rows = list(reader)

    frame_contract_ratio: float | None = None
    dropped = next((item[2] for item in VARIANTS if item[0] == variant), None)
    if variant.startswith("C0_"):
        if frame_rows:
            errors.append("single-lidar run unexpectedly contains multi-lidar frame rows")
    else:
        if len(frame_rows) < 1100:
            errors.append("multi-lidar frame_stats has fewer than 1100 frames")
        expected_present = {0, 1, 2, 3} - ({dropped} if dropped is not None else set())
        expected_missing = {dropped} if dropped is not None else set()
        valid_rows = 0
        previous_end = -math.inf
        for row in frame_rows:
            try:
                begin = float(row["begin_time"])
                end = float(row["end_time"])
                merged = int(row["merged_points"])
                partial = int(row["partial"])
                present = {int(value) for value in row["present_lidar_ids"].split(";") if value}
                missing = {int(value) for value in row["missing_lidar_ids"].split(";") if value}
                point_counts = {lidar_id: int(row[f"points_lidar_{lidar_id}"]) for lidar_id in range(4)}
            except (KeyError, TypeError, ValueError):
                continue
            row_valid = (
                math.isfinite(begin) and math.isfinite(end) and begin < end and end > previous_end
                and merged >= 1000 and present == expected_present and missing == expected_missing
                and partial == int(dropped is not None)
                and all(point_counts[lidar_id] >= 100 for lidar_id in expected_present)
                and all(point_counts[lidar_id] == 0 for lidar_id in expected_missing)
            )
            valid_rows += int(row_valid)
            previous_end = end
        frame_contract_ratio = valid_rows / len(frame_rows) if frame_rows else 0.0
        if frame_contract_ratio < 0.99:
            errors.append("fewer than 99% of multi-lidar frames satisfy the configured source/content contract")
        if lidar_times and frame_rows:
            try:
                last_frame_end = float(frame_rows[-1]["end_time"])
                if abs(lidar_times[-1] - last_frame_end) > 0.20:
                    errors.append("trajectory does not reach the final assembled lidar frame")
            except (KeyError, ValueError):
                errors.append("invalid final frame timestamp")

    try:
        map_points = pcd_point_count(map_pcd)
    except (OSError, UnicodeError, ValueError) as error:
        map_points = 0
        errors.append(str(error))
    if map_points < 10_000:
        errors.append("map contains fewer than 10000 points")
    return {
        "passed": not errors,
        "errors": errors[:100],
        "trajectory_pose_count": len(lidar_times),
        "trajectory_start_s": lidar_times[0] if lidar_times else None,
        "trajectory_end_s": lidar_times[-1] if lidar_times else None,
        "frame_stats_count": len(frame_rows),
        "frame_contract_ratio": frame_contract_ratio,
        "map_point_count": map_points,
        "output_sha256": {path.name: sha256(path) for path in required},
    }


def valid_completed_run(run_dir: Path, fingerprint: str, variant: str) -> bool:
    metadata_path = run_dir / "controller_metadata.json"
    manifest_path = run_dir / "completion_manifest.json"
    if not metadata_path.is_file() or not manifest_path.is_file():
        return False
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    validation = validate_outputs(run_dir, variant)
    return (
        metadata.get("returncode") == 0
        and metadata.get("experiment_fingerprint") == fingerprint
        and manifest.get("experiment_fingerprint") == fingerprint
        and manifest.get("variant") == variant
        and manifest.get("output_validation", {}).get("passed") is True
        and manifest.get("output_validation", {}).get("output_sha256") == validation.get("output_sha256")
        and validation["passed"]
    )


def first_valid_attempt(cell_dir: Path, fingerprint: str, variant: str) -> Path | None:
    for attempt_dir in sorted(cell_dir.glob("attempt_[0-9][0-9]")):
        if valid_completed_run(attempt_dir, fingerprint, variant):
            return attempt_dir
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--experiment-root", type=Path, required=True)
    parser.add_argument("--bag", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--wsl", default=str(Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / "wsl.exe"))
    parser.add_argument("--distro", default="Ubuntu-22.04")
    parser.add_argument("--cpu-set", default="0-7")
    parser.add_argument("--cell-timeout-s", type=float, default=1800.0)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    args.base_config = args.base_config.resolve()
    args.experiment_root = args.experiment_root.resolve()
    args.bag = args.bag.resolve()
    args.repo = args.repo.resolve()
    if args.cell_timeout_s <= 0.0 or args.max_attempts < 1:
        raise SystemExit("timeout and max attempts must be positive")
    if is_within(args.experiment_root, args.repo):
        raise SystemExit("experiment root must be outside the code repository")
    if is_within(args.experiment_root, args.bag) or is_within(args.bag, args.experiment_root):
        raise SystemExit("experiment root and input bag must not overlap")

    base = yaml.safe_load(args.base_config.read_text(encoding="utf-8"))
    configs_dir = args.experiment_root / "configs"
    runs_root = args.experiment_root / "runs"
    state_dir = args.experiment_root / "formal_matrix"
    configs_dir.mkdir(parents=True, exist_ok=True)
    runs_root.mkdir(parents=True, exist_ok=True)
    state_dir.mkdir(parents=True, exist_ok=True)
    lock_path = state_dir / "controller.lock"
    try:
        lock_fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as error:
        raise SystemExit(f"experiment root is locked: {lock_path}") from error
    os.write(
        lock_fd,
        json.dumps({"pid": os.getpid(), "host": socket.gethostname(), "started_at": now()}).encode("utf-8"),
    )
    os.close(lock_fd)
    def release_lock() -> None:
        try:
            lock_path.unlink()
        except FileNotFoundError:
            pass
    atexit.register(release_lock)

    config_records = []
    for name, mode, dropped_id in VARIANTS:
        config_path = configs_dir / f"{name}.yaml"
        rendered = yaml.safe_dump(
            build_config(base, mode, dropped_id), allow_unicode=True, sort_keys=False, width=120
        )
        if config_path.exists() and config_path.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"refusing to change frozen config: {config_path}")
        if not config_path.exists():
            config_path.write_text(rendered, encoding="utf-8")
        config_records.append({"variant": name, "path": str(config_path), "sha256": sha256(config_path)})

    # The standard install entry is a Linux symlink to repo/bin and cannot be
    # stat'ed reliably by Windows Python. Fingerprint the real executable.
    binary = args.repo / "bin" / "run_frontend_offline"
    template = args.repo / "scripts" / "run_frontend_offline.sh"
    if not binary.is_file() or not template.is_file() or not args.bag.is_dir():
        raise SystemExit("binary, run template, or bag is missing")
    git = ["git", "-c", f"safe.directory={args.repo.as_posix()}", "-C", str(args.repo)]
    git_head = subprocess.run(
        [*git, "rev-parse", "HEAD"], check=True, text=True, capture_output=True
    ).stdout.strip()
    git_status = subprocess.run(
        [*git, "status", "--porcelain"], check=True, text=True, capture_output=True
    ).stdout
    if git_status.strip():
        raise SystemExit("formal matrix requires a clean Lightning-LM worktree")
    git_diff = subprocess.run(
        [*git, "diff", "--binary", "HEAD"], check=True, text=True, capture_output=True
    ).stdout
    bag_files = sorted(
        path for path in args.bag.iterdir()
        if path.is_file() and (path.name == "metadata.yaml" or path.suffix in {".db3", ".mcap"})
    )
    if not bag_files or not any(path.name == "metadata.yaml" for path in bag_files):
        raise SystemExit("input bag metadata/storage files are missing")
    bag_identity = [
        {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)} for path in bag_files
    ]
    environment_script = (
        "set -e; source /opt/ros/humble/setup.bash; "
        f"source {to_wsl(args.repo / 'install' / 'setup.bash')}; "
        "printf 'ros_distro=%s\\n' \"${ROS_DISTRO:-}\"; "
        "printf 'os='; . /etc/os-release; printf '%s %s\\n' \"$ID\" \"$VERSION_ID\"; "
        "printf 'prefix='; ros2 pkg prefix lightning; "
        f"ldd {to_wsl(binary)} | sed -E 's/ \\(0x[0-9a-f]+\\)//g'"
    )
    environment_identity = subprocess.run(
        [args.wsl, "-d", args.distro, "--", "bash", "-lc", environment_script],
        check=True, text=True, capture_output=True,
    ).stdout
    expected_prefix = to_wsl(args.repo / "install" / "lightning")
    if f"prefix={expected_prefix}" not in environment_identity:
        raise SystemExit("ros2 package prefix does not resolve to the frozen install tree")
    controller_path = Path(__file__).resolve()
    fingerprint_payload = {
        "seed": SEED,
        "base_config": {"path": str(args.base_config), "sha256": sha256(args.base_config)},
        "configs": config_records,
        "binary": {"path": str(binary), "sha256": sha256(binary)},
        "template": {"path": str(template), "sha256": sha256(template)},
        "controller": {"path": str(controller_path), "sha256": sha256(controller_path)},
        "install_setup": {
            "path": str(args.repo / "install" / "setup.bash"),
            "sha256": sha256(args.repo / "install" / "setup.bash"),
        },
        "git_head": git_head,
        "git_status_porcelain": git_status.splitlines(),
        "git_diff_binary_sha256": text_sha256(git_diff),
        "bag": {"path": str(args.bag), "files": bag_identity},
        "environment": {
            "distro": args.distro,
            "identity": environment_identity.splitlines(),
            "identity_sha256": text_sha256(environment_identity),
        },
        "cpu_set": args.cpu_set,
        "repeats": 3,
        "cell_timeout_s": args.cell_timeout_s,
    }
    fingerprint = hashlib.sha256(
        json.dumps(fingerprint_payload, ensure_ascii=False, sort_keys=True).encode("utf-8")
    ).hexdigest()
    fingerprint_payload["experiment_fingerprint"] = fingerprint
    fingerprint_path = state_dir / "fingerprint.json"
    if fingerprint_path.exists() and json.loads(fingerprint_path.read_text(encoding="utf-8")) != fingerprint_payload:
        raise SystemExit("experiment fingerprint changed; use a new experiment root")
    atomic_json(fingerprint_path, fingerprint_payload)

    schedule = []
    for repeat in range(1, 4):
        block = [{"variant": variant, "repeat": repeat, "block": repeat} for variant, _, _ in VARIANTS]
        random.Random(SEED + repeat).shuffle(block)
        schedule.extend(block)
    for index, item in enumerate(schedule, start=1):
        item["order_index"] = index
    schedule_payload = {"seed": SEED, "design": "three independently randomized complete blocks", "items": schedule}
    schedule_path = state_dir / "schedule.json"
    if schedule_path.exists() and json.loads(schedule_path.read_text(encoding="utf-8")) != schedule_payload:
        raise SystemExit("frozen schedule changed; use a new experiment root")
    atomic_json(schedule_path, schedule_payload)
    if args.dry_run:
        print(json.dumps({"fingerprint": fingerprint, "schedule": schedule}, ensure_ascii=False, indent=2))
        return 0

    failures = 0
    progress = []
    config_by_variant = {record["variant"]: Path(record["path"]) for record in config_records}
    for item in schedule:
        variant, repeat = item["variant"], item["repeat"]
        cell_dir = runs_root / variant / f"repeat_{repeat:02d}"
        cell_dir.mkdir(parents=True, exist_ok=True)
        selected = first_valid_attempt(cell_dir, fingerprint, variant)
        if selected is not None:
            selection = {
                **item, "status": "skipped_success", "cell_dir": str(cell_dir),
                "selected_attempt": str(selected),
            }
            atomic_json(cell_dir / "selected_attempt.json", selection)
            progress.append(selection)
            atomic_json(state_dir / "progress.json", {"updated_at": now(), "items": progress})
            continue
        prior_attempts = sorted(cell_dir.glob("attempt_[0-9][0-9]"))
        if len(prior_attempts) >= args.max_attempts:
            failures += 1
            progress.append({
                **item, "status": "attempt_limit_reached", "cell_dir": str(cell_dir),
                "prior_attempts": [str(path) for path in prior_attempts],
            })
            atomic_json(state_dir / "progress.json", {"updated_at": now(), "items": progress})
            continue
        attempt_index = len(prior_attempts) + 1
        run_dir = cell_dir / f"attempt_{attempt_index:02d}"
        if run_dir.exists():
            raise SystemExit(f"attempt allocation collision: {run_dir}")
        run_dir.mkdir(parents=True)
        config_path = config_by_variant[variant]
        time_log = run_dir / "gnu_time.txt"
        command = [
            args.wsl, "-d", args.distro, "--",
            "timeout", "--signal=INT", "--kill-after=30s", f"{math.ceil(args.cell_timeout_s)}s",
            "/usr/bin/time", "-v", "-o", to_wsl(time_log),
            "bash", to_wsl(template),
            "--bag", to_wsl(args.bag),
            "--config", to_wsl(config_path),
            "--output-dir", to_wsl(run_dir),
            "--sequence", variant,
            "--repeat", str(repeat),
            "--cpu-set", args.cpu_set,
            "--cpu-count", "8",
            "--experiment-fingerprint", fingerprint,
        ]
        started = now()
        running_metadata = {
            **item,
            "variant": variant,
            "repeat": repeat,
            "attempt": attempt_index,
            "status": "running",
            "started_at": started,
            "experiment_fingerprint": fingerprint,
            "config": str(config_path),
            "config_sha256": sha256(config_path),
            "run_dir": str(run_dir),
            "command": command,
            "cell_timeout_s": args.cell_timeout_s,
        }
        atomic_json(run_dir / "controller_metadata.json", running_metadata)
        controller_stdout = run_dir / "controller.stdout.log"
        controller_stderr = run_dir / "controller.stderr.log"
        timed_out = False
        with controller_stdout.open("wb") as stdout, controller_stderr.open("wb") as stderr:
            try:
                result = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=args.cell_timeout_s + 90.0)
                returncode = result.returncode
            except subprocess.TimeoutExpired:
                timed_out = True
                returncode = 124
        output_validation = validate_outputs(run_dir, variant)
        metadata = {
            **running_metadata,
            "status": "finished",
            "finished_at": now(),
            "returncode": returncode,
            "timed_out": timed_out,
            "output_validation": output_validation,
        }
        atomic_json(run_dir / "controller_metadata.json", metadata)
        if returncode == 0 and output_validation["passed"]:
            atomic_json(
                run_dir / "completion_manifest.json",
                {
                    "created_at": now(),
                    "experiment_fingerprint": fingerprint,
                    "variant": variant,
                    "repeat": repeat,
                    "attempt": attempt_index,
                    "output_validation": output_validation,
                },
            )
        success = returncode == 0 and valid_completed_run(run_dir, fingerprint, variant)
        status = "success" if success else "failed"
        if not success:
            failures += 1
        progress_item = {**metadata, "status": status}
        if success:
            atomic_json(
                cell_dir / "selected_attempt.json",
                {**item, "status": "success", "cell_dir": str(cell_dir), "selected_attempt": str(run_dir)},
            )
        progress.append(progress_item)
        atomic_json(state_dir / "progress.json", {"updated_at": now(), "items": progress})

    print(json.dumps({"runs": len(schedule), "failures": failures, "fingerprint": fingerprint}, ensure_ascii=False))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
