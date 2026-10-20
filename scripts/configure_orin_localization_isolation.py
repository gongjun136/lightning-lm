#!/usr/bin/env python3
"""Preview/apply Orin CPU 10-11 scheduler and managed IRQ isolation.

Changes only the primary extlinux APPEND line, backs up exact original bytes,
and never reboots. Root is required for --apply. No nohz_full is requested:
the inspected field kernel does not support CONFIG_NO_HZ_FULL.
"""
import argparse
import datetime
import os
from pathlib import Path
import re
import shutil
import tempfile


def proposed(text):
    lines = text.splitlines(keepends=True)
    primary = False
    matched = 0
    for index, line in enumerate(lines):
        label = re.match(r"\s*LABEL\s+(\S+)", line)
        if label:
            primary = label[1] == "primary"
        if not primary or not re.match(r"\s*APPEND\s+", line):
            continue
        matched += 1
        for key, value in (("isolcpus", "domain,managed_irq,10-11"), ("irqaffinity", "0-9")):
            existing = re.search(rf"\b{key}=(\S+)", line)
            if existing and existing[1] != value:
                raise ValueError(f"conflicting existing {key}; review instead of overwriting")
            if not existing:
                line = line.rstrip("\r\n") + f" {key}={value}\n"
        lines[index] = line
    if matched != 1:
        raise ValueError("expected exactly one primary APPEND line")
    return "".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    model = Path("/proc/device-tree/model").read_bytes().decode().strip("\0")
    if "Jetson AGX Orin" not in model or not Path("/sys/devices/system/cpu/cpu11").exists():
        raise SystemExit("This deployment profile requires the inspected 12-core AGX Orin.")
    path = Path("/boot/extlinux/extlinux.conf")
    original = path.read_text()
    updated = proposed(original)
    print(updated)
    if args.apply and updated != original:
        if os.geteuid() != 0:
            raise SystemExit("--apply requires sudo")
        backup = path.with_name(path.name + ".before-localization-" +
                                datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        shutil.copy2(path, backup)
        fd, temporary = tempfile.mkstemp(prefix=".extlinux-localization-", dir=path.parent)
        try:
            with os.fdopen(fd, "w") as stream:
                stream.write(updated)
                stream.flush()
                os.fsync(stream.fileno())
            shutil.copymode(path, temporary)
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        print(f"BACKUP={backup}\nReboot required. Roll back by restoring this backup and rebooting.")
