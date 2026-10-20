#!/usr/bin/env python3
"""Exercise the dependency gate with real temporary Git repositories."""
import contextlib
import io
from pathlib import Path
import subprocess
import tempfile
import unittest

import yaml

from check_workspace_repos import check


class DependencyGateTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.repo = self.root / "consumer"
        self.repo.mkdir()
        self.dep = self.root / "src" / "common_msgs"
        self.dep.mkdir(parents=True)
        self.git("init")
        self.git("config", "user.name", "Test")
        self.git("config", "user.email", "test@example.invalid")
        (self.dep / "message").write_text("old interface\n")
        self.git("add", ".")
        self.git("commit", "-m", "old interface")
        self.old = self.git("rev-parse", "HEAD")
        self.manifest(self.old)

    def git(self, *args):
        return subprocess.check_output(
            ["git", "-C", str(self.dep), *args], stderr=subprocess.STDOUT, text=True
        ).strip()

    def manifest(self, version, filename=None):
        data = {"repositories": {"common_msgs": {
            "type": "git", "url": "example.invalid/common_msgs.git", "version": version}}}
        for name in ([filename] if filename else ["workspace.repos", "common_msgs.repos"]):
            (self.repo / name).write_text(yaml.safe_dump(data))

    def result(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return check(self.repo, self.root)

    def test_matching_dependency_passes(self):
        self.assertEqual(self.result(), 0)

    def test_newer_local_interface_cannot_hide_stale_pin(self):
        (self.dep / "message").write_text("new active_fault field\n")
        self.git("commit", "-am", "new interface")
        self.assertEqual(self.result(), 1)
        self.manifest(self.git("rev-parse", "HEAD"))
        self.assertEqual(self.result(), 0)

    def test_conflicting_manifests_fail(self):
        self.manifest("0" * 40, "common_msgs.repos")
        self.assertEqual(self.result(), 1)

    def test_dirty_dependency_fails(self):
        (self.dep / "message").write_text("uncommitted interface\n")
        self.assertEqual(self.result(), 1)

    def test_missing_dependency_fails(self):
        self.dep.rename(self.dep.with_name("missing"))
        self.assertEqual(self.result(), 1)

    def test_branch_instead_of_pin_fails(self):
        self.manifest("main")
        self.assertEqual(self.result(), 1)


if __name__ == "__main__":
    unittest.main()
