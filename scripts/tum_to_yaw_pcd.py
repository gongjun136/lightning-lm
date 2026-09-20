#!/usr/bin/env python3
"""Convert a TUM pose trajectory to an ASCII PCD with intensity equal to yaw in degrees."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def yaw_degrees(qx: float, qy: float, qz: float, qw: float) -> float:
    yaw = math.atan2(2.0 * (qw * qz + qx * qy), 1.0 - 2.0 * (qy * qy + qz * qz))
    return math.degrees(yaw) % 360.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="input TUM trajectory")
    parser.add_argument("--output", type=Path, required=True, help="output ASCII PCD")
    args = parser.parse_args()

    points: list[tuple[float, float, float, float]] = []
    with args.input.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            values = [float(value) for value in line.split()]
            if len(values) != 8 or not all(math.isfinite(value) for value in values):
                raise ValueError(f"invalid TUM pose at line {line_number}")
            _, x, y, z, qx, qy, qz, qw = values
            norm = math.sqrt(qx * qx + qy * qy + qz * qz + qw * qw)
            if norm <= 1e-12:
                raise ValueError(f"zero quaternion at line {line_number}")
            qx, qy, qz, qw = qx / norm, qy / norm, qz / norm, qw / norm
            points.append((x, y, z, yaw_degrees(qx, qy, qz, qw)))

    if not points:
        raise ValueError("input trajectory contains no valid poses")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# .PCD v0.7 - trajectory with intensity=yaw_deg_0_360\n")
        stream.write("VERSION 0.7\n")
        stream.write("FIELDS x y z intensity\n")
        stream.write("SIZE 4 4 4 4\n")
        stream.write("TYPE F F F F\n")
        stream.write("COUNT 1 1 1 1\n")
        stream.write(f"WIDTH {len(points)}\n")
        stream.write("HEIGHT 1\n")
        stream.write("VIEWPOINT 0 0 0 1 0 0 0\n")
        stream.write(f"POINTS {len(points)}\n")
        stream.write("DATA ascii\n")
        for x, y, z, yaw in points:
            stream.write(f"{x:.12f} {y:.12f} {z:.12f} {yaw:.12f}\n")

    print(f"points={len(points)} output={args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
