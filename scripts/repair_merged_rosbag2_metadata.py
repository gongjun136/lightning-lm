#!/usr/bin/env python3
"""Normalize metadata produced by the bulk ROS 2 bag merger."""

import argparse
from pathlib import Path

from merge_rosbag2_by_timestamp import repair_metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bag", nargs="+", type=Path)
    args = parser.parse_args()
    for bag in args.bag:
        repair_metadata(bag.resolve(), update_contract=True)
        print(bag.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
