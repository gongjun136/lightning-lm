#!/usr/bin/env python3
"""Exercise launcher failure paths without starting ROS nodes or hardware."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import yaml


REPO = Path(__file__).resolve().parents[1]
MANIFEST_TOOL = Path(os.environ.get("LIGHTNING_TEST_MAP_MANIFEST_TOOL",
                     str(REPO.parent / "lightning_lm_ws/build/lightning_lm/bin/export_map_manifest")))


class LauncherFailureTest(unittest.TestCase):
    def run_launcher(self, probe_exit, algorithm_exit, fixed_map_transform=None, missing_manifest=False, aligned_map=False):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "setup.bash").write_text(":\n")
            (root / "vehicle.yaml").write_text('VIN: "SANLV77EXSH000858"\nLIDAR: "MASTER AND SLAVE"\n')
            (root / "map").mkdir()
            (root / "map/index.txt").write_text("test index\n")
            (root / "map/global.pcd").write_text("test geometry\n")
            if not missing_manifest:
                subprocess.run([str(MANIFEST_TOOL), "--map=" + str(root / "map"),
                                "--status=not_required"], check=True, capture_output=True)
                if fixed_map_transform == "1" or aligned_map:
                    # Deliberately differs from the shared Bauma transform.
                    manifest_path = root / "map/map_manifest.yaml"
                    manifest = yaml.safe_load(manifest_path.read_text())
                    alignment = manifest["legacy_alignment"]
                    alignment.update(status="aligned", needs_transform_to_old=True,
                                     target_map_id="old_test_map", target_global_pcd_sha256="0" * 64)
                    alignment["transform"]["translation_xyz"] = [10.0, 20.0, 30.0]
                    alignment["transform"]["matrix_4x4"] = [
                        [1., 0., 0., 10.], [0., 1., 0., 20.],
                        [0., 0., 1., 30.], [0., 0., 0., 1.]]
                    manifest["output"]["fixed_map_transform"].update(
                        enabled=True, translation_xyz=[10.0, 20.0, 30.0])
                    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
            fake = root / "ros2"
            fake.write_text("""#!/usr/bin/env bash
if [[ "$1 $2 $3" == 'run lightning_lm export_map_manifest' ]]; then
  exec "$MAP_MANIFEST_TOOL" "${@:4}"
elif [[ "$1 $2 $3" == 'run lightning_lm check_localization_inputs' ]]; then
  echo 'mock input delivery check'
  exit "${PROBE_EXIT}"
elif [[ "$1 $2 $3" == 'run lightning_lm run_loc_online' ]]; then
  echo 'mock algorithm failure' >&2
  exit "${ALGORITHM_EXIT}"
elif [[ "$1 $2" == 'topic list' ]]; then
  exit 0
fi
echo "unexpected ROS command: $*" >&2
exit 91
""")
            fake.chmod(0o755)
            (root / "tegrastats").write_text("#!/bin/sh\nexit 0\n")
            (root / "tegrastats").chmod(0o755)
            env = dict(os.environ)
            # Do not inherit field overrides from the shell running this test.
            for name in list(env):
                if name.startswith(("SANY_", "LIGHTNING_LM_")):
                    del env[name]
            env.update({
                "PATH": str(root) + os.pathsep + env["PATH"],
                "LIGHTNING_LM_REPO_DIR": str(REPO),
                "LIGHTNING_LM_ROS_SETUP": str(root / "setup.bash"),
                "LIGHTNING_LM_INSTALL_SETUP": str(root / "setup.bash"),
                "LIGHTNING_LM_OUT_ROOT": str(root / "runs"),
                "LIGHTNING_LM_RUN_MODE": "production",
                "SANY_VEHICLE_CONFIG": str(root / "vehicle.yaml"),
                "SANY_MAP_PATH": str(root / "map"),
                "PROBE_EXIT": str(probe_exit),
                "ALGORITHM_EXIT": str(algorithm_exit),
                "MAP_MANIFEST_TOOL": str(MANIFEST_TOOL),
            })
            if fixed_map_transform is not None:
                env["SANY_ENABLE_FIXED_MAP_TRANSFORM"] = fixed_map_transform
            # Allow cold Python/NumPy imports and metadata I/O on mounted workspaces.
            result = subprocess.run(
                ["bash", str(REPO / "scripts/run_sany_lidar_loc.sh"), "test"],
                env=env, text=True, capture_output=True, timeout=30)
            run = root / "runs/test"
            files = {str(path.relative_to(run)): path.read_text()
                     for path in run.rglob("*") if path.is_file()}
            return result, files

    def test_failed_input_check_is_visible_and_never_launches_algorithm(self):
        result, files = self.run_launcher(7, 0)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("required sensor deliveries are unavailable", result.stderr)
        self.assertIn("logs/input_check.log", result.stderr)
        self.assertIn("mock input delivery check", files["logs/input_check.log"])
        self.assertNotIn("logs/run_loc_online.stderr.log", files)

    def test_algorithm_exit_code_and_logs_are_preserved(self):
        result, files = self.run_launcher(0, 23)
        self.assertEqual(result.returncode, 23, result.stderr)
        self.assertIn("algorithm_exit_code=23", files["run_metadata.txt"])
        self.assertIn("mock algorithm failure", files["logs/run_loc_online.stderr.log"])

    def test_site_transform_is_enabled_only_in_run_snapshot(self):
        result, files = self.run_launcher(0, 0, "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        snapshot = yaml.safe_load(files["config.yaml"])
        self.assertTrue(snapshot["output"]["fixed_map_transform"]["enabled"])
        self.assertEqual(snapshot["output"]["fixed_map_transform"]["translation_xyz"], [10., 20., 30.])
        self.assertIn("map_manifest", snapshot["output"]["fixed_map_transform"]["provenance"])
        shared = yaml.safe_load((REPO / "config/sany/localization.yaml").read_text())
        self.assertFalse(shared["output"]["fixed_map_transform"]["enabled"])

    def test_missing_map_manifest_fails_before_sensor_or_algorithm_start(self):
        result, files = self.run_launcher(0, 0, missing_manifest=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("map manifest validation failed", result.stderr)
        self.assertIn("missing", files["logs/map_manifest.log"])
        self.assertNotIn("logs/input_check.log", files)
        self.assertNotIn("logs/run_loc_online.stderr.log", files)

    def test_native_map_disables_shared_site_transform_in_snapshot(self):
        result, files = self.run_launcher(0, 0)
        self.assertEqual(result.returncode, 0, result.stderr)
        snapshot = yaml.safe_load(files["config.yaml"])
        self.assertFalse(snapshot["output"]["fixed_map_transform"]["enabled"])
        self.assertEqual(snapshot["output"]["fixed_map_transform"]["translation_xyz"], [0., 0., 0.])
        self.assertIn("manifest_sha256", snapshot["map_package"])

    def test_legacy_switch_cannot_override_selected_map(self):
        result, files = self.run_launcher(0, 0, "0", aligned_map=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        snapshot = yaml.safe_load(files["config.yaml"])
        self.assertTrue(snapshot["output"]["fixed_map_transform"]["enabled"])
        self.assertIn("deprecated and ignored", result.stderr)

    def test_illegal_site_transform_fails_before_algorithm_start(self):
        result, files = self.run_launcher(0, 0, "yes")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("SANY_ENABLE_FIXED_MAP_TRANSFORM must be 0 or 1", result.stderr)
        self.assertNotIn("logs/run_loc_online.stderr.log", files)


if __name__ == "__main__":
    unittest.main()
