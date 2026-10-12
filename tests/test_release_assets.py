#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Exercise the filenames actually preserved by upload-artifact."""

import hashlib
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("recovery", ROOT / "tools/ci_verify_release_recovery.py")
recovery = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(recovery)


class Recovery(unittest.TestCase):
    def setUp(self):
        self.tag = "v2.4.1"
        self.commit = "a" * 40
        self.run = {"event": "push", "head_branch": self.tag, "head_sha": self.commit,
                    "path": recovery.WORKFLOW, "status": "completed"}
        self.jobs = [{"name": name, "conclusion": "success"} for name in recovery.REQUIRED_JOBS]
        self.jobs.append({"name": "publish release", "conclusion": "failure"})

    def test_only_original_validated_tag_build_is_accepted(self):
        recovery.validate_run(self.run, self.jobs, self.tag, self.commit)
        for field in self.run:
            with self.subTest(field=field), self.assertRaises(ValueError):
                recovery.validate_run({**self.run, field: "wrong"}, self.jobs, self.tag, self.commit)
        with self.assertRaises(ValueError):
            recovery.validate_run(self.run, self.jobs[1:], self.tag, self.commit)
        with self.assertRaises(ValueError):
            recovery.validate_run(self.run, self.jobs + [{"name": "tests", "conclusion": "failure"}],
                                  self.tag, self.commit)

    def test_sbom_must_equal_verified_signed_predicate(self):
        with tempfile.TemporaryDirectory() as directory:
            dist = Path(directory)
            for name in ("linux.tar.gz", "macos.tar.gz", "windows.zip"):
                (dist / name).touch()
                (dist / (name + ".spdx.json")).write_text('{}')
            verified = [{"verificationResult": {"statement": {"predicate": {}}}}]
            with patch.object(recovery, "gh_json", return_value=verified) as gh:
                recovery.verify_assets(dist, self.tag, self.commit)
                self.assertEqual(gh.call_count, 6)
                self.assertIn(self.commit, gh.call_args.args)
                self.assertIn("refs/tags/" + self.tag, gh.call_args.args)
                (dist / "linux.tar.gz.spdx.json").write_text('{"changed": true}')
                with self.assertRaisesRegex(ValueError, "differs from the signed predicate"):
                    recovery.verify_assets(dist, self.tag, self.commit)


@unittest.skipUnless(sys.platform != "win32" and shutil.which("bash"), "Release publishing uses Bash")
class ReleaseAssets(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.dist = self.root / "dist"
        self.dist.mkdir()
        self.names = [f"mbelib-neo-2.4.1-{suffix}" for suffix in (
            "linux-x64.tar.gz", "macos-arm64.tar.gz", "windows-x64.zip")]
        for name in self.names:
            (self.dist / name).write_bytes(name.encode())
            (self.dist / (name + ".spdx.json")).write_text("{}\n")

    def prepare(self):
        return subprocess.run(["bash", str(ROOT / "tools/ci_prepare_release_assets.sh")],
                              cwd=self.root, env={**os.environ, "ARCHIVE_PREFIX": "mbelib-neo-2.4.1"},
                              capture_output=True, text=True)

    def test_recovery_requires_valid_explicit_dispatch_inputs(self):
        for event, tag, run_id in (("push", "v2.4.1", "123"),
                                   ("workflow_dispatch", "invalid", "123"),
                                   ("workflow_dispatch", "v2.4.1", "")):
            result = subprocess.run(["bash", str(ROOT / "tools/ci_validate_release_tag.sh")],
                                    env={**os.environ, "GITHUB_EVENT_NAME": event,
                                         "RELEASE_RECOVERY_TAG": tag, "RELEASE_RECOVERY_RUN_ID": run_id},
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Release recovery requires", result.stderr)

    def test_package_named_sboms_and_checksums(self):
        result = self.prepare()
        self.assertEqual(result.returncode, 0, result.stderr)
        expected = "".join(f"{hashlib.sha256(name.encode()).hexdigest()}  {name}\n"
                           for name in sorted(self.names))
        self.assertEqual((self.dist / "SHA256SUMS.txt").read_text(), expected)
        for name in self.names:
            self.assertEqual((self.dist / (name + ".spdx.json")).read_text(), "{}\n")

    def test_missing_package_or_sbom_fails(self):
        for name in (self.names[0], self.names[0] + ".spdx.json"):
            path = self.dist / name
            content = path.read_bytes()
            path.unlink()
            self.assertNotEqual(self.prepare().returncode, 0)
            self.assertFalse((self.dist / "SHA256SUMS.txt").exists())
            path.write_bytes(content)

    def test_wrong_sbom_name_fails(self):
        (self.dist / (self.names[0] + ".spdx.json")).rename(self.dist / "unrelated.spdx.json")
        self.assertNotEqual(self.prepare().returncode, 0)
        self.assertFalse((self.dist / "SHA256SUMS.txt").exists())

    def test_unsupported_architecture_fails(self):
        for suffix in ("", ".spdx.json"):
            (self.dist / (self.names[0] + suffix)).rename(
                self.dist / (self.names[0].replace("x64", "unknown") + suffix))
        self.assertNotEqual(self.prepare().returncode, 0)


if __name__ == "__main__":
    unittest.main()
