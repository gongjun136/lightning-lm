#!/usr/bin/env python3
"""Create a reviewable SANY compute candidate from the frozen deployment YAML."""
import argparse
from pathlib import Path

import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", choices=("conservative", "balanced"), default="conservative")
    parser.add_argument("--allow-experimental-rotation", action="store_true",
                        help="balanced failed the first SANY accuracy gate; research use only")
    parser.add_argument("--lio-points", type=int, default=2200)
    parser.add_argument("--minimum-lio-points", type=int, default=1500)
    parser.add_argument("--ndt-points", type=int, default=3500)
    args = parser.parse_args()
    if args.profile == "balanced" and not args.allow_experimental_rotation:
        parser.error("balanced failed SANY yaw/RPE gates; use conservative for deployment, or explicitly "
                     "pass --allow-experimental-rotation for further experiments")
    if args.base.resolve() == args.output.resolve():
        parser.error("output must differ from the frozen base YAML")
    if args.lio_points < 100 or not 100 <= args.ndt_points <= 1000000:
        parser.error("point budgets must be >=100; NDT cap must be <=1000000")
    if not 100 <= args.minimum_lio_points <= args.lio_points:
        parser.error("minimum LIO points must be between 100 and the top LIO budget")
    config = yaml.safe_load(args.base.read_text(encoding="utf-8"))
    multi = config["multi_lidar"]
    if not multi.get("enabled") or len([k for k in multi["topics"] if k.startswith("lidar_")]) != 3:
        parser.error("the candidate requires the three-LiDAR configuration")
    load = multi["adaptive_load"]
    balanced = args.profile == "balanced"
    # Never remove observations to satisfy timing. Conservative keeps all three
    # sources; balanced alternates the two sides at every 10 Hz state update.
    load.update(enabled=True, tracking_min_lidars=2 if balanced else 3,
                tracking_lidar_count=2 if balanced else 3,
                relocalization_min_lidars=3, cloud_publish_min_lidars=3,
                rotate_secondary_lidars=balanced,
                point_strides=[1, 1, 1],
                lio_point_budgets=[args.lio_points, max(args.minimum_lio_points, round(args.lio_points * 0.82)),
                                   args.minimum_lio_points],
                target_latency_sec=0.09, hard_latency_sec=0.3,
                degrade_processing_ratio=0.72, recover_processing_ratio=0.5,
                degrade_consecutive_frames=1, recover_consecutive_frames=75,
                predictive=True)
    config.setdefault("compute_budget", {})["ndt_max_points"] = args.ndt_points
    config["fasterlio"]["skip_lidar_num"] = 0
    config["system"]["enable_lidar_loc_skip"] = False
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "# Experimental compute candidate; requires target-Orin full-stack acceptance.\n"
        f"# Base: {args.base.name}; profile: {args.profile}\n"
        + yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
