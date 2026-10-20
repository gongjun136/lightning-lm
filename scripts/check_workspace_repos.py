#!/usr/bin/env python3
"""Read-only check that workspace dependency sources match the delivery manifests."""
import argparse
from pathlib import Path
import re
import subprocess
import sys

import yaml


def check(repo, workspace):
    dependencies = yaml.safe_load((repo / "workspace.repos").read_text())["repositories"]
    messages = yaml.safe_load((repo / "common_msgs.repos").read_text())["repositories"]
    errors = []
    for name, spec in messages.items():
        if dependencies.get(name) != spec:
            errors.append(f"{name}: common_msgs.repos differs from workspace.repos")
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
        print("PASS: dependency sources match both manifests (installed overlays not checked)")
    return 1 if errors else 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    sys.exit(check(Path(__file__).resolve().parents[1], args.workspace.resolve()))
