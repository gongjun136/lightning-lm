#!/usr/bin/env python3
"""Pin one vehicle/calibration revision into a complete Lightning-LM YAML.

Requires PyYAML and NumPy (provided by the ROS development environment).
This resolves relative LiDAR extrinsics only; it does not infer IMU/body/CAN calibration.
"""
import argparse
import copy
import hashlib
from pathlib import Path
import re

import numpy as np
import yaml


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    value = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected YAML mapping: {path}")
    return value


def identifier(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", value):
        raise ValueError(f"Invalid vehicle/calibration identifier: {value!r}")
    return value


def resolve(registry_root, vehicle_id, calibration_id, base_config):
    root = Path(registry_root).resolve()
    vehicle_id, calibration_id = identifier(vehicle_id), identifier(calibration_id)
    vehicle_dir = (root / vehicle_id).resolve()
    revision = (vehicle_dir / "calibrations" / calibration_id).resolve()
    if root not in vehicle_dir.parents or vehicle_dir not in revision.parents:
        raise ValueError("Vehicle/calibration path escapes registry")
    vehicle_path = vehicle_dir / "vehicle.yaml"
    metadata_path = revision / "metadata.yaml"
    extrinsics_path = revision / "extrinsics.yaml"
    for path in (vehicle_path, metadata_path, extrinsics_path):
        if root not in path.resolve().parents:
            raise ValueError("Registry file escapes registry")
    vehicle, metadata = load(vehicle_path), load(metadata_path)
    if vehicle.get("vehicle_id") != vehicle_id:
        raise ValueError("Vehicle identity mismatch")
    if metadata.get("vehicle_id") != vehicle_id or metadata.get("calibration_id") != calibration_id:
        raise ValueError("Calibration identity/revision mismatch")
    if metadata.get("transform_convention") != "T_primary_from_lidar_row_major" or metadata.get("translation_unit") != "m":
        raise ValueError("Unsupported transform convention or translation unit")
    digest = sha256(extrinsics_path)
    if digest != metadata.get("extrinsics_sha256"):
        raise ValueError("Calibration SHA-256 mismatch; create a new revision instead of editing one")
    base = load(base_config)
    multi = base.get("multi_lidar", {})
    primary = vehicle.get("primary_lidar_id")
    if not multi.get("enabled") or multi.get("primary_lidar_id") != primary or metadata.get("primary_lidar_id") != primary:
        raise ValueError("Primary LiDAR/configuration mismatch")
    sensors = vehicle.get("sensors", {})
    extrinsics = load(extrinsics_path).get("extrinsics", {})
    if not sensors or set(sensors) != set(extrinsics):
        raise ValueError("Sensor/extrinsic key mismatch")
    expected_topics = {}
    for key, sensor in sensors.items():
        if key != f"lidar{sensor['id']}":
            raise ValueError("Sensor key/ID mismatch")
        index = sensor["id"]
        expected_topics[f"lidar_{index}"] = sensor["lidar_topic"]
        expected_topics[f"imu_{index}"] = sensor["imu_topic"]
        values = extrinsics[key].get("T")
        if not isinstance(values, list) or len(values) != 16:
            raise ValueError(f"{key}: expected a row-major 4x4 matrix")
        matrix = np.asarray(values, dtype=float).reshape(4, 4)
        rotation = matrix[:3, :3]
        if not np.isfinite(matrix).all() or not np.allclose(matrix[3], [0, 0, 0, 1], rtol=0, atol=1e-9):
            raise ValueError(f"{key}: invalid homogeneous transform")
        if not np.allclose(rotation.T @ rotation, np.eye(3), rtol=0, atol=1e-6) or abs(np.linalg.det(rotation) - 1) > 1e-6:
            raise ValueError(f"{key}: rotation is not SO(3)")
        if index == primary and not np.allclose(matrix, np.eye(4), rtol=0, atol=1e-9):
            raise ValueError("Primary LiDAR transform must be identity")
    if multi.get("topics") != expected_topics:
        raise ValueError("Topic/ID mapping mismatch; IP addresses alone do not identify a vehicle")
    if base.get("common", {}).get("lidar_topic") != sensors[f"lidar{primary}"]["lidar_topic"] or base.get("common", {}).get("imu_topic") != sensors[f"lidar{primary}"]["imu_topic"]:
        raise ValueError("Common primary LiDAR/IMU topic mismatch")
    result = copy.deepcopy(base)
    result["multi_lidar"]["extrinsics"] = copy.deepcopy(extrinsics)
    result["calibration_reference"] = {
        "vehicle_id": vehicle_id, "vehicle_identity_status": vehicle.get("identity_status"),
        "calibration_id": calibration_id, "extrinsics_sha256": digest,
        "vehicle_profile_sha256": sha256(vehicle_path), "metadata_sha256": sha256(metadata_path),
        "base_config_sha256": sha256(base_config), "registry_revision": str(revision),
        "transform_convention": metadata["transform_convention"], "translation_unit": "m",
        "scope": "relative_lidar_extrinsics_only", "other_vehicle_parameters": "inherited_unverified",
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry-root", required=True, type=Path)
    parser.add_argument("--vehicle-id", required=True)
    parser.add_argument("--calibration-id", required=True)
    parser.add_argument("--base-config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = resolve(args.registry_root, args.vehicle_id, args.calibration_id, args.base_config)
    except (ValueError, KeyError, TypeError, OSError, yaml.YAMLError) as error:
        parser.exit(2, f"Calibration configuration rejected: {error}\n")
    if args.output.exists():
        parser.exit(2, "Output already exists; keep archived run configurations immutable\n")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(result, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(f"Resolved {args.vehicle_id}/{args.calibration_id} -> {args.output}")


if __name__ == "__main__":
    main()
