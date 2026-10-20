#!/usr/bin/env python3
"""Read-only check that workspace dependency sources match the delivery manifest."""
import argparse
from pathlib import Path
import re
import subprocess
import sys

import yaml


def check_interfaces(repo, workspace):
    """Check cross-repository source contracts, not just matching Git hashes."""
    try:
        loc = yaml.safe_load((repo / "config/sany/localization.yaml").read_text())["diagnostics"]
        driver = workspace / "src/livox_ros_driver2"
        heartbeat = (workspace / "src/common_msgs/diagnostic_monitor_interfaces/msg/NodeHeartbeat.msg").read_text()
        lidar_state = (workspace / "src/common_msgs/livox_ros_driver/msg/LidarDiagnosticState.msg").read_text()
        errors = []
        for field in ("bool bounded_time_sync_warning", "int64 time_sync_offset_ns",
                      "uint64 time_sync_evidence_age_us"):
            if not re.search(r"^" + field.replace(" ", r"\s+") + r"\s*(?:#.*)?$", lidar_state, re.MULTILINE):
                errors.append(f"LidarDiagnosticState.msg lacks {field}")
        source = (driver / "src/livox_ros_driver2.cpp").read_text()
        if '"livox-dtc-v3-bounded-sync"' not in source:
            errors.append("radar bounded-warning evidence requires livox-dtc-v3-bounded-sync")
        if not re.search(r"^uint32\[\]\s+active_fault\s*(?:#.*)?$", heartbeat, re.MULTILINE):
            errors.append("NodeHeartbeat.msg lacks uint32[] active_fault")
        ids = []
        for role in ("front", "rear", "left", "right"):
            launch = (driver / f"launch_ROS2/lidar_{role}.launch.py").read_text()
            match = re.search(r"^heartbeat_node_id\s*=\s*(\d+)\s*$", launch, re.MULTILINE)
            ids.append(int(match[1]) if match else None)
            if f"'/diagnostics/heartbeat/livox_{role}'" not in launch:
                errors.append(f"{role}: heartbeat topic mismatch")
        if None in ids or len(set(ids)) != 4 or ids != loc.get("radar_role_source_ids"):
            errors.append("radar role source IDs mismatch")
        if "dtc::PublishDiagnosticSnapshot(" not in source:
            errors.append("driver does not update its own heartbeat DTCs")
        if "lidar_dtc_aggregator" in (driver / "launch_ROS2/master_lidar.launch.py").read_text():
            errors.append("master launch still starts the unregistered aggregate node")
        if "radar_source_id" in loc or "internal_radar_topic" in loc:
            errors.append("localization still depends on the retired aggregate source")
    except (OSError, KeyError, TypeError, yaml.YAMLError) as exc:
        errors = [f"cannot inspect interface contract: {exc}"]
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    if not errors:
        print("PASS: active_fault, bounded PTP schema/version and four per-device heartbeat identities match (runtime still requires verification)")
    return 1 if errors else 0


def check(repo, workspace):
    dependencies = yaml.safe_load((repo / "workspace.repos").read_text())["repositories"]
    errors = []
    for name, spec in dependencies.items():
        version = str(spec.get("version", ""))
        if spec.get("type") != "git" or not re.fullmatch(r"[0-9a-f]{40}", version):
            errors.append(f"{name}: expected an exact 40-character Git commit")
            continue
        path = workspace / "src" / name
        def git(*args):
            return subprocess.check_output(
                ["git", "-C", str(path), *args], text=True, stderr=subprocess.STDOUT
            ).strip()
        try:
            # Reject a missing checkout accidentally resolving to a parent repository.
            if Path(git("rev-parse", "--show-toplevel")).resolve() != path.resolve():
                raise ValueError("dependency path is not a repository root")
            head = git("rev-parse", "HEAD")
            if head != version:
                errors.append(f"{name}: HEAD={head}, manifest={version}")
            if git("status", "--porcelain"):
                errors.append(f"{name}: dirty checkout; commit/review changes before validation")
            print(f"{name}: HEAD={head}, manifest={version}")
        except (subprocess.CalledProcessError, ValueError) as exc:
            detail = exc.output.strip() if isinstance(exc, subprocess.CalledProcessError) else str(exc)
            errors.append(f"{name}: cannot inspect checkout: {detail}")
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    if not errors:
        print("PASS: dependency sources match workspace.repos (installed overlays not checked)")
    return 1 if errors else 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    repo, workspace = Path(__file__).resolve().parents[1], args.workspace.resolve()
    versions = check(repo, workspace)
    interfaces = check_interfaces(repo, workspace)
    sys.exit(versions or interfaces)
