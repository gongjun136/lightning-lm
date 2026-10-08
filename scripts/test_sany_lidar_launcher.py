#!/usr/bin/env python3
"""Exercise launcher failure paths without starting ROS nodes or hardware."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]


class LauncherFailureTest(unittest.TestCase):
    def run_launcher(self, probe_exit, algorithm_exit):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "setup.bash").write_text(":\n")
            fake = root / "ros2"
            fake.write_text("""#!/usr/bin/env bash
if [[ "$1 $2 $3" == 'run lightning_lm check_localization_inputs' ]]; then
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
                "PROBE_EXIT": str(probe_exit),
                "ALGORITHM_EXIT": str(algorithm_exit),
            })
            result = subprocess.run(
                ["bash", str(REPO / "scripts/run_sany_lidar_loc.sh"), "test"],
                env=env, text=True, capture_output=True, timeout=15)
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


if __name__ == "__main__":
    unittest.main()
