#!/usr/bin/env python3
"""Resolve shared SANY parameters using the vehicle's explicit VIN and LIDAR layout."""
import argparse
import copy
import hashlib
from pathlib import Path

import numpy as np
import yaml


def load(path):
    value = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected YAML mapping: {path}")
    return value


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def resolve(base_path, vehicle_path, calibration_path=None):
    base_path, vehicle_path = Path(base_path).resolve(), Path(vehicle_path).resolve()
    base, vehicle = load(base_path), load(vehicle_path)
    vin = vehicle.get("VIN")
    layout = vehicle.get("LIDAR")
    if not isinstance(vin, str) or not vin.strip():
        raise ValueError("vehicle_config.yaml requires a non-empty VIN")
    vin = vin.strip()
    layouts = {"MASTER AND SLAVE": 3, "MASTER": 4}
    if layout not in layouts:
        raise ValueError("LIDAR must be 'MASTER AND SLAVE' or 'MASTER'")
    count = layouts[layout]
    declaration = base.get("sany_vehicle_config", {})
    calibration_path = Path(calibration_path) if calibration_path else (
        base_path.parent / declaration["calibrations"])
    calibration_path = calibration_path.resolve()
    registry = load(calibration_path)
    if (registry.get("transform_convention") != "T_primary_from_lidar_row_major"
            or registry.get("translation_unit") != "m" or registry.get("primary_lidar_id") != 0):
        raise ValueError("Unsupported SANY calibration convention, unit or primary LiDAR")
    if vin not in registry.get("vehicles", {}):
        raise ValueError(f"Unknown VIN: {vin}; add its verified calibration before launch")
    profile = registry["vehicles"][vin]
    extrinsics = profile["extrinsics"]
    for index in range(4):
        matrix = np.asarray(extrinsics[f"lidar{index}"]["T"], dtype=float).reshape(4, 4)
        rotation = matrix[:3, :3]
        if (not np.isfinite(matrix).all()
                or not np.allclose(matrix[3], [0, 0, 0, 1], rtol=0, atol=1e-9)
                or not np.allclose(rotation.T @ rotation, np.eye(3), rtol=0, atol=1e-6)
                or abs(np.linalg.det(rotation) - 1) > 1e-6):
            raise ValueError(f"Invalid calibration transform: {vin}/lidar{index}")
        if index == 0 and not np.allclose(matrix, np.eye(4), rtol=0, atol=1e-9):
            raise ValueError("Primary LiDAR calibration must be identity")
    result = copy.deepcopy(base)
    result.pop("sany_vehicle_config", None)
    multi = result["multi_lidar"]
    if not multi.get("enabled") or multi.get("primary_lidar_id") != 0:
        raise ValueError("SANY requires enabled multi_lidar and primary_lidar_id: 0")
    # Topics describe the hardware once; only active sensors survive in the run snapshot.
    topics = multi["topics"]
    multi["topics"] = {key: topics[key] for index in range(count)
                       for key in (f"lidar_{index}", f"imu_{index}")}
    multi["extrinsics"] = {f"lidar{index}": copy.deepcopy(extrinsics[f"lidar{index}"])
                            for index in range(count)}
    if multi.get("adaptive_load", {}).get("enabled"):
        multi["adaptive_load"]["cloud_publish_min_lidars"] = count
        multi["adaptive_load"]["cloud_publish_require_primary"] = True
    result["vehicle_configuration"] = {
        "VIN": vin, "DEVICE_CODE": vehicle.get("DEVICE_CODE"), "LIDAR": layout,
        "lidar_count": count, "calibration_id": profile["calibration_id"],
        "vehicle_config_path": str(vehicle_path), "vehicle_config_sha256": digest(vehicle_path),
        "base_config_path": str(base_path), "base_config_sha256": digest(base_path),
        "calibrations_path": str(calibration_path), "calibrations_sha256": digest(calibration_path),
        "transform_convention": registry["transform_convention"], "translation_unit": "m",
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", required=True, type=Path)
    parser.add_argument("--vehicle-config", type=Path, default=Path("/opt/conf/vehicle_config.yaml"))
    parser.add_argument("--calibrations", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ins-site", action="store_true", help="Resolve an INS site's shared LiDAR cloud configuration")
    args = parser.parse_args()
    try:
        if args.ins_site:
            result = load(args.base_config)
            cloud = result.get("cloud", {})
            if cloud.get("enabled") and cloud.get("lidar_config"):
                lidar_path = (args.base_config.parent / cloud["lidar_config"]).resolve()
                if load(lidar_path).get("sany_vehicle_config"):
                    lidar_result = resolve(lidar_path, args.vehicle_config, args.calibrations)
                    lidar_output = args.output.parent / "lidar_config.yaml"
                    lidar_output.parent.mkdir(parents=True, exist_ok=True)
                    with lidar_output.open("x", encoding="utf-8") as stream:
                        yaml.safe_dump(lidar_result, stream, sort_keys=False, allow_unicode=True)
                    cloud["lidar_config"] = str(lidar_output.resolve())
                else:
                    cloud["lidar_config"] = str(lidar_path)
        else:
            result = resolve(args.base_config, args.vehicle_config, args.calibrations)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", encoding="utf-8") as stream:
            yaml.safe_dump(result, stream, sort_keys=False, allow_unicode=True)
    except (ValueError, KeyError, TypeError, OSError, yaml.YAMLError) as error:
        parser.exit(2, f"SANY vehicle configuration rejected: {error}\n")
    print(str(args.output.resolve()))


if __name__ == "__main__":
    main()
