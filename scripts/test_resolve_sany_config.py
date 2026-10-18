#!/usr/bin/env python3
"""Verify VIN/layout selection, shared parameters and fail-closed vehicle startup."""
import copy
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import yaml
from resolve_sany_config import resolve

REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "config/sany"


class SanyConfigTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.vehicle = Path(self.temp.name) / "vehicle.yaml"

    def vehicle_config(self, vin="SANLV77EXSH000858", layout="MASTER AND SLAVE"):
        self.vehicle.write_text(yaml.safe_dump({"VIN": vin, "LIDAR": layout}))

    def test_each_vin_and_layout_for_both_modes(self):
        registry = yaml.safe_load((CONFIG / "calibrations.yaml").read_text())
        for mode in ("mapping", "localization"):
            for vin, profile in registry["vehicles"].items():
                for layout, count in (("MASTER AND SLAVE", 3), ("MASTER", 4)):
                    with self.subTest(mode=mode, vin=vin, layout=layout):
                        self.vehicle_config(vin, layout)
                        result = resolve(CONFIG / f"{mode}.yaml", self.vehicle)
                        multi = result["multi_lidar"]
                        self.assertEqual(len(multi["topics"]), count * 2)
                        self.assertEqual(multi["extrinsics"], {f"lidar{i}": profile["extrinsics"][f"lidar{i}"] for i in range(count)})
                        self.assertEqual(multi["match_tolerance"], 0.02)
                        self.assertEqual(multi["min_lidars"], 1)
                        self.assertEqual(result["fasterlio"]["livox_point_time_scale"], 1e-9)
                        self.assertNotIn("sany_vehicle_config", result)
                        expected_scale = {"SANLV77EXSH000858": 0.00120639253574024,
                                          "SANLU77E8TH000006": 0.0011687}
                        self.assertEqual(result["system"]["wheel_speed_scale_mps_per_rpm"],
                                         expected_scale[vin])
                        if mode == "localization":
                            load = multi["adaptive_load"]
                            self.assertEqual(load["cloud_publish_min_lidars"], count)
                            self.assertEqual(load["tracking_min_lidars"], 1)
                            self.assertEqual(load["relocalization_min_lidars"], 1)
                            self.assertEqual(load["tracking_lidar_count"], 0)
                            self.assertEqual(load["lio_point_budgets"], [2200, 1804, 1500])

    def test_algorithm_parameters_are_shared_across_vehicles(self):
        self.vehicle_config()
        old = resolve(CONFIG / "localization.yaml", self.vehicle)
        self.vehicle_config("SANLU77E8TH000006", "MASTER")
        new = resolve(CONFIG / "localization.yaml", self.vehicle)
        for result in (old, new):
            result.pop("vehicle_configuration")
            result["system"].pop("wheel_speed_scale_mps_per_rpm")
            result["multi_lidar"].pop("extrinsics")
            result["multi_lidar"].pop("topics")
            result["multi_lidar"]["adaptive_load"].pop("cloud_publish_min_lidars")
        self.assertEqual(old, new)
        self.assertEqual(old["relocalization"]["backend"], "solid_kiss")
        self.assertEqual(old["relocalization"]["refinement_backend"], "plane_icp")
        self.assertFalse(old["output"]["fixed_map_transform"]["enabled"])

    def test_real_new_vehicle_vin_selects_confirmed_calibration(self):
        self.vehicle_config("SANLU77E8TH000006", "MASTER")
        result = resolve(CONFIG / "localization.yaml", self.vehicle)
        self.assertEqual(result["vehicle_configuration"]["VIN"], "SANLU77E8TH000006")
        self.assertEqual(result["vehicle_configuration"]["calibration_id"], "user_provided_20261010_v1")
        self.assertEqual(result["vehicle_configuration"]["lidar_count"], 4)
        self.assertIn("lidar_3", result["multi_lidar"]["topics"])

    def test_missing_unknown_or_illegal_vehicle_fails(self):
        with self.assertRaises(OSError):
            resolve(CONFIG / "localization.yaml", self.vehicle)
        for vin, layout in (("UNKNOWN", "MASTER"), (None, "MASTER"),
                            ("SW956ELCG0858", "MASTER"), ("SANLU77E8TH000006", "SLAVE")):
            self.vehicle_config(vin, layout)
            with self.assertRaises(ValueError):
                resolve(CONFIG / "localization.yaml", self.vehicle)

    def test_site_output_transform_override_preserves_shared_default(self):
        self.vehicle_config()
        original = (CONFIG / "localization.yaml").read_bytes()
        default = resolve(CONFIG / "localization.yaml", self.vehicle)
        enabled = resolve(CONFIG / "localization.yaml", self.vehicle, fixed_map_transform=1)
        disabled = resolve(CONFIG / "localization.yaml", self.vehicle, fixed_map_transform=0)
        self.assertFalse(default["output"]["fixed_map_transform"]["enabled"])
        self.assertTrue(enabled["output"]["fixed_map_transform"]["enabled"])
        self.assertFalse(disabled["output"]["fixed_map_transform"]["enabled"])
        self.assertTrue(enabled["vehicle_configuration"]["fixed_map_transform_override"])
        self.assertEqual(enabled["output"]["fixed_map_transform"]["translation_xyz"],
                         default["output"]["fixed_map_transform"]["translation_xyz"])
        self.assertEqual((CONFIG / "localization.yaml").read_bytes(), original)

    def test_invalid_calibration_is_rejected(self):
        registry = yaml.safe_load((CONFIG / "calibrations.yaml").read_text())
        broken = copy.deepcopy(registry)
        broken["vehicles"]["SANLU77E8TH000006"]["extrinsics"]["lidar1"]["T"][0] = 2
        path = Path(self.temp.name) / "broken.yaml"
        path.write_text(yaml.safe_dump(broken))
        self.vehicle_config("SANLU77E8TH000006", "MASTER")
        with self.assertRaises(ValueError):
            resolve(CONFIG / "localization.yaml", self.vehicle, path)

    def test_invalid_wheel_speed_scale_is_rejected(self):
        registry = yaml.safe_load((CONFIG / "calibrations.yaml").read_text())
        path = Path(self.temp.name) / "broken.yaml"
        self.vehicle_config("SANLU77E8TH000006", "MASTER")
        for value in (None, 0, -1, float("nan"), float("inf"), True, "0.0011687"):
            with self.subTest(value=value):
                broken = copy.deepcopy(registry)
                profile = broken["vehicles"]["SANLU77E8TH000006"]
                if value is None:
                    profile.pop("wheel_speed_scale_mps_per_rpm")
                else:
                    profile["wheel_speed_scale_mps_per_rpm"] = value
                path.write_text(yaml.safe_dump(broken))
                with self.assertRaisesRegex(ValueError, "wheel_speed_scale_mps_per_rpm"):
                    resolve(CONFIG / "localization.yaml", self.vehicle, path)

    def test_run_snapshot_is_immutable(self):
        self.vehicle_config()
        output = Path(self.temp.name) / "runtime.yaml"
        command = [sys.executable, str(REPO / "scripts/resolve_sany_config.py"), "--base-config", str(CONFIG / "mapping.yaml"),
                   "--vehicle-config", str(self.vehicle), "--output", str(output)]
        self.assertEqual(subprocess.run(command, capture_output=True).returncode, 0)
        original = output.read_bytes()
        self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
        self.assertEqual(output.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
