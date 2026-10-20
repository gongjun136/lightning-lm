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
        radar = yaml.safe_load((driver / "config/dtc_aggregator.yaml").read_text())["lidar_dtc_aggregator"]["ros__parameters"]
        heartbeat = (workspace / "src/common_msgs/diagnostic_monitor_interfaces/msg/NodeHeartbeat.msg").read_text()
        errors = []
        if not re.search(r"^uint32\[\]\s+active_fault\s*(?:#.*)?$", heartbeat, re.MULTILINE):
            errors.append("NodeHeartbeat.msg lacks uint32[] active_fault")
        topic = radar.get("internal_output_topic")
        if not topic or topic != loc.get("internal_radar_topic"):
            errors.append(f"radar topic mismatch: producer={topic!r}, consumer={loc.get('internal_radar_topic')!r}")
        if radar.get("functional_safety.node_id") != loc.get("radar_source_id"):
            errors.append("radar aggregator source ID mismatch")
        if radar.get("expected_source_ids") != loc.get("radar_role_source_ids"):
            errors.append("radar role source IDs mismatch")
        # A matching YAML key alone is insufficient if an old producer ignores it.
        source = (driver / "src/lidar_dtc_aggregator.cpp").read_text()
        if not re.search(r'declare_parameter<std::string>\(\s*"internal_output_topic"', source):
            errors.append("radar producer does not declare internal_output_topic")
    except (OSError, KeyError, TypeError, yaml.YAMLError) as exc:
        errors = [f"cannot inspect interface contract: {exc}"]
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    if not errors:
        print("PASS: heartbeat field and radar topic/source-ID contracts match (runtime still requires verification)")
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
