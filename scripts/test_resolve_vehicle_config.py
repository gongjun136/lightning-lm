#!/usr/bin/env python3
"""Checks for rejecting mismatched/tampered vehicle calibration archives."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

import yaml

spec = importlib.util.spec_from_file_location("resolver", Path(__file__).with_name("resolve_vehicle_config.py"))
resolver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(resolver)


class ResolveVehicleTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.revision = self.root / "vehicle_a/calibrations/v1"
        self.revision.mkdir(parents=True)
        self.identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        sensor = {"id": 0, "lidar_topic": "/lidar", "imu_topic": "/imu"}
        self.vehicle = {"vehicle_id": "vehicle_a", "primary_lidar_id": 0, "sensors": {"lidar0": sensor}}
        self.metadata = {"vehicle_id": "vehicle_a", "calibration_id": "v1", "primary_lidar_id": 0,
                         "transform_convention": "T_primary_from_lidar_row_major", "translation_unit": "m"}
        self.base = {"common": {"lidar_topic": "/lidar", "imu_topic": "/imu"},
                     "backend": {"btc": {"max_drift_ratio": 0.03}},
                     "multi_lidar": {"enabled": True, "primary_lidar_id": 0,
                                     "topics": {"lidar_0": "/lidar", "imu_0": "/imu"},
                                     "extrinsics": {"lidar0": {"T": []}}}}
        self.write(self.root / "vehicle_a/vehicle.yaml", self.vehicle)
        self.write(self.root / "base.yaml", self.base)
        self.set_extrinsics({"lidar0": {"T": self.identity}})

    def write(self, path, data):
        path.write_text(yaml.safe_dump(data), encoding="utf-8")

    def set_extrinsics(self, value):
        path = self.revision / "extrinsics.yaml"
        self.write(path, {"extrinsics": value})
        self.metadata["extrinsics_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        self.write(self.revision / "metadata.yaml", self.metadata)

    def resolve(self):
        return resolver.resolve(self.root, "vehicle_a", "v1", self.root / "base.yaml")

    def test_overlay_preserves_runtime_parameters(self):
        result = self.resolve()
        self.assertEqual(result["multi_lidar"]["extrinsics"]["lidar0"]["T"], self.identity)
        result.pop("calibration_reference")
        result["multi_lidar"]["extrinsics"] = copy.deepcopy(self.base["multi_lidar"]["extrinsics"])
        self.assertEqual(result, self.base)

    def test_rejects_wrong_vehicle_revision_and_topic_mapping(self):
        for field, value in [("vehicle_id", "vehicle_b"), ("calibration_id", "v2")]:
            original = self.metadata[field]
            self.metadata[field] = value
            self.write(self.revision / "metadata.yaml", self.metadata)
            with self.assertRaisesRegex(ValueError, "identity/revision"):
                self.resolve()
            self.metadata[field] = original
        self.write(self.revision / "metadata.yaml", self.metadata)
        self.base["multi_lidar"]["topics"]["lidar_0"] = "/other"
        self.write(self.root / "base.yaml", self.base)
        with self.assertRaisesRegex(ValueError, "Topic/ID"):
            self.resolve()

    def test_rejects_tampered_matrix_bytes(self):
        with (self.revision / "extrinsics.yaml").open("a") as stream:
            stream.write("# edited\n")
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            self.resolve()

    def test_rejects_reflection_nonfinite_bad_bottom_and_nonidentity_primary(self):
        for index, value, message in [(0, -1, "SO\\(3\\)"), (3, float("nan"), "homogeneous"),
                                      (15, 0, "homogeneous"), (3, 0.2, "identity")]:
            matrix = self.identity.copy(); matrix[index] = value
            self.set_extrinsics({"lidar0": {"T": matrix}})
            with self.assertRaisesRegex(ValueError, message):
                self.resolve()

    def test_rejects_missing_sensors_and_path_traversal(self):
        self.set_extrinsics({})
        with self.assertRaisesRegex(ValueError, "key mismatch"):
            self.resolve()
        with self.assertRaisesRegex(ValueError, "identifier"):
            resolver.resolve(self.root, "../vehicle_a", "v1", self.root / "base.yaml")


if __name__ == "__main__":
    unittest.main()
