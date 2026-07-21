#!/usr/bin/env python3
"""Capture the host/WSL environment used by the formal experiment suite."""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
from datetime import datetime
from pathlib import Path


def run(command: list[str]) -> dict[str, object]:
    result = subprocess.run(command, text=True, encoding="utf-8", errors="replace", capture_output=True)
    return {
        "command": subprocess.list2cmdline(command),
        "returncode": result.returncode,
        "stdout": result.stdout.strip(),
        "stderr": result.stderr.strip(),
    }


def wsl(distro: str, command: str) -> dict[str, object]:
    return run(["wsl", "-d", distro, "--", "bash", "-lc", command])


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=repo / "runs" / "formal_report_20260717" / "_environment" / "environment.json",
    )
    args = parser.parse_args()
    payload = {
        "schema_version": 1,
        "captured_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "windows_python": platform.platform(),
        "host": {
            "cpu": run([
                "powershell", "-NoProfile", "-Command",
                "Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors | ConvertTo-Json -Compress",
            ]),
            "memory": run([
                "powershell", "-NoProfile", "-Command",
                "Get-CimInstance Win32_ComputerSystem | Select-Object TotalPhysicalMemory | ConvertTo-Json -Compress",
            ]),
            "gpu": run([
                "nvidia-smi", "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader,nounits",
            ]),
        },
        "ubuntu_22_04": wsl(
            "Ubuntu-22.04",
            "source /opt/ros/humble/setup.bash; "
            "printf 'ROS_DISTRO='; printenv ROS_DISTRO || true; "
            "uname -a; lscpu; free -b; "
            "printf 'cmake='; cmake --version | head -n1; "
            "printf 'gcc='; gcc --version | head -n1; "
            "printf 'ros2='; ros2 --version 2>&1 || true",
        ),
        "ubuntu_20_04": wsl(
            "Ubuntu-20.04",
            "source /opt/ros/noetic/setup.bash; "
            "printf 'ROS_DISTRO='; printenv ROS_DISTRO || true; "
            "uname -a; "
            "printf 'cmake='; cmake --version | head -n1; "
            "printf 'gcc='; gcc --version | head -n1; "
            "printf 'rosversion='; rosversion -d",
        ),
        "repository": wsl(
            "Ubuntu-22.04",
            f"cd '{repo.as_posix().replace('F:/', '/mnt/f/')}' && "
            "git rev-parse HEAD; git branch --show-current; git status --short",
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
