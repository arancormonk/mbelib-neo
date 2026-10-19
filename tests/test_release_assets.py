#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Exercise the filenames actually preserved by upload-artifact."""

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
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
            for name in ("mbelib-neo-2.4.1-linux-x64.tar.gz", "mbelib-neo-2.4.1-macos-arm64.tar.gz", "mbelib-neo-2.4.1-windows-x64.zip"):
                (dist / name).touch()
                (dist / (name + ".spdx.json")).write_text('{}')
            verified = [{"verificationResult": {"statement": {"predicate": {}}}}]
            with patch.object(recovery, "gh_json", return_value=verified) as gh:
                recovery.verify_assets(dist, self.tag, self.commit)
                self.assertEqual(gh.call_count, 6)
                self.assertIn(self.commit, gh.call_args.args)
                self.assertIn("refs/tags/" + self.tag, gh.call_args.args)
                (dist / "mbelib-neo-2.4.1-linux-x64.tar.gz.spdx.json").write_text('{"changed": true}')
                with self.assertRaisesRegex(ValueError, "differs from the signed predicate"):
                    recovery.verify_assets(dist, self.tag, self.commit)

    def test_rebuild_requires_every_original_tag_check(self):
        with self.assertRaisesRegex(ValueError, "from 2.6.0"):
            recovery.validate_run(self.run, self.jobs, self.tag, self.commit, rebuild=True)
        self.tag = "v2.6.0"
        self.run["head_branch"] = self.tag
        jobs = [{"name": name, "conclusion": "success"} for name in recovery.TAG_CHECKS]
        jobs.append({"name": "package release (${{ matrix.os }})", "conclusion": "skipped"})
        recovery.validate_run(self.run, jobs, self.tag, self.commit, rebuild=True)
        for name in recovery.TAG_CHECKS:
            with self.subTest(missing=name), self.assertRaises(ValueError):
                recovery.validate_run(self.run, [j for j in jobs if j["name"] != name],
                                      self.tag, self.commit, rebuild=True)
        for conclusion in ("failure", "skipped", "cancelled"):
            invalid = [{**j, "conclusion": conclusion} if j["name"] == "clang-tidy" else j for j in jobs]
            with self.subTest(conclusion=conclusion), self.assertRaises(ValueError):
                recovery.validate_run(self.run, invalid, self.tag, self.commit, rebuild=True)
        with self.assertRaises(ValueError):
            recovery.validate_run({**self.run, "head_sha": "b" * 40}, jobs, self.tag, self.commit, rebuild=True)

    def test_recovery_requires_the_dispatch_check_from_the_release_that_has_it(self):
        dispatch = "linux • Debug • SIMD (SSE2 baseline + AVX2 dispatch)"
        self.assertNotIn(dispatch, recovery.tag_checks("2.6.0"))
        self.assertIn(dispatch, recovery.tag_checks("2.6.1"))
        self.assertIn(dispatch, recovery.tag_checks("3.0.0"))
        with self.assertRaises(ValueError):
            recovery.tag_checks("2.6")
        for tag in ("v2.6.1", "v2.7.0"):
            run = {**self.run, "head_branch": tag}
            jobs = [{"name": name, "conclusion": "success"} for name in recovery.tag_checks(tag[1:])]
            for options in ({"rebuild": True}, {"aur_only": True}):
                with self.subTest(tag=tag, options=options):
                    recovery.validate_run(run, jobs, tag, self.commit, **options)
                    without = [j for j in jobs if j["name"] != dispatch]
                    with self.assertRaisesRegex(ValueError, "lacks successful"):
                        recovery.validate_run(run, without, tag, self.commit, **options)
                    for conclusion in ("failure", "skipped", "cancelled"):
                        invalid = [{**j, "conclusion": conclusion} if j["name"] == dispatch else j for j in jobs]
                        with self.subTest(conclusion=conclusion), self.assertRaises(ValueError):
                            recovery.validate_run(run, invalid, tag, self.commit, **options)
        # Releases from before the check stay recoverable without it.
        run = {**self.run, "head_branch": "v2.6.0"}
        jobs = [{"name": name, "conclusion": "success"} for name in recovery.TAG_CHECKS]
        recovery.validate_run(run, jobs, "v2.6.0", self.commit, rebuild=True)
        recovery.validate_run(run, jobs, "v2.6.0", self.commit, aur_only=True)

    def test_rebuild_identity_requires_reviewed_main_workflow(self):
        run = {"event": "workflow_dispatch", "head_branch": "main", "path": recovery.WORKFLOW, "head_sha": "b" * 40}
        with patch.object(recovery, "gh_json", side_effect=[run, {"status": "ahead"}]):
            self.assertEqual(recovery.rebuild_identity("123"), ("refs/heads/main", "b" * 40))
        for field in run:
            with self.subTest(field=field), patch.object(recovery, "gh_json", return_value={**run, field: "wrong"}), self.assertRaises(ValueError):
                recovery.rebuild_identity("123")
        with patch.object(recovery, "gh_json", side_effect=[run, {"status": "diverged"}]), self.assertRaises(ValueError):
            recovery.rebuild_identity("123")
        with self.assertRaises(ValueError):
            recovery.rebuild_identity("not-a-run")

    def test_aur_recovery_accepts_failed_aur_but_requires_every_tag_check(self):
        tag = "v2.6.0"
        run = {**self.run, "head_branch": tag}
        jobs = [{"name": name, "conclusion": "success"} for name in recovery.TAG_CHECKS]
        jobs.append({"name": "Update AUR package (release)", "conclusion": "failure"})
        recovery.validate_run(run, jobs, tag, self.commit, aur_only=True)
        for name in recovery.TAG_CHECKS:
            with self.subTest(missing=name), self.assertRaises(ValueError):
                recovery.validate_run(run, [j for j in jobs if j["name"] != name], tag, self.commit, aur_only=True)
        for conclusion in ("failure", "skipped", "cancelled"):
            invalid = [{**j, "conclusion": conclusion} if j["name"] == "clang-tidy" else j for j in jobs]
            with self.subTest(conclusion=conclusion), self.assertRaises(ValueError):
                recovery.validate_run(run, invalid, tag, self.commit, aur_only=True)
        for options in ({"rebuild": True}, {"rebuild": True, "aur_only": True}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                recovery.validate_run(run, jobs, tag, self.commit, **options)

    def test_rebuild_signed_sbom_binds_original_tag_source(self):
        with tempfile.TemporaryDirectory() as directory:
            dist = Path(directory)
            sbom = {"packages": [{"SPDXID": "SPDXRef-Package-mbelib-neo", "versionInfo": "2.6.0",
                                  "sourceInfo": f"Release commit {self.commit}; see legal/status.md"}]}
            for suffix in ("linux-x64.tar.gz", "macos-arm64.tar.gz", "windows-x64.zip", "source.tar.gz"):
                archive = dist / ("mbelib-neo-2.6.0-" + suffix)
                archive.touch()
                Path(str(archive) + ".spdx.json").write_text(json.dumps(sbom))
            verified = [{"verificationResult": {"statement": {"predicate": sbom}}}]
            builder = ("refs/heads/main", "b" * 40)
            with patch.object(recovery, "gh_json", return_value=verified) as gh:
                recovery.verify_assets(dist, "v2.6.0", self.commit, builder)
                self.assertEqual(gh.call_count, 12)
                self.assertIn(builder[0], gh.call_args.args)
                self.assertIn(builder[1], gh.call_args.args)
                with self.assertRaisesRegex(ValueError, "validated release source"):
                    recovery.verify_assets(dist, "v2.6.0", "c" * 40, builder)


@unittest.skipUnless(sys.platform.startswith("linux") and shutil.which("bash"),
                     "Release publisher uses Linux Bash and sha256sum")
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

    def test_rebuild_requires_explicit_main_dispatch(self):
        for ref in ("refs/tags/v2.6.0", "refs/heads/unreviewed"):
            result = subprocess.run(["bash", str(ROOT / "tools/ci_validate_release_tag.sh")],
                                    env={**os.environ, "GITHUB_EVENT_NAME": "workflow_dispatch", "GITHUB_REF": ref,
                                         "RELEASE_RECOVERY_TAG": "v2.6.0", "RELEASE_RECOVERY_RUN_ID": "123",
                                         "RELEASE_REBUILD_PACKAGES": "true"}, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("reviewed workflow on main", result.stderr)

    def test_package_named_sboms_and_checksums(self):
        result = self.prepare()
        self.assertEqual(result.returncode, 0, result.stderr)
        names = sorted(self.names + [name + ".spdx.json" for name in self.names])
        expected = "".join(f"{hashlib.sha256((self.dist / name).read_bytes()).hexdigest()}  {name}\n"
                           for name in names)
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


    def test_new_releases_require_source_and_all_sboms(self):
        from release_assets import collect, prepare
        for old in list(self.dist.iterdir()):
            old.rename(old.with_name(old.name.replace("2.4.1", "2.6.0")))
        with self.assertRaises(ValueError):
            collect(self.dist, "2.6.0")
        source = self.dist / "mbelib-neo-2.6.0-source.tar.gz"
        source.write_bytes(b"source")
        with self.assertRaises(ValueError):
            collect(self.dist, "2.6.0")
        Path(str(source) + ".spdx.json").write_text("{}")
        packages, sboms = collect(self.dist, "2.6.0")
        self.assertEqual((len(packages), len(sboms)), (4, 4))
        self.assertEqual(len(prepare(self.dist, "2.6.0").splitlines()), 8)
        with patch.object(recovery, "gh_json", return_value=[{"verificationResult": {"statement": {"predicate": {}}}}]) as gh:
            recovery.verify_assets(self.dist, "v2.6.0", "a" * 40)
            self.assertEqual(gh.call_count, 12)
        with patch.object(recovery, "gh_json", side_effect=ValueError("tampered artifact")):
            with self.assertRaisesRegex(ValueError, "tampered"):
                recovery.verify_assets(self.dist, "v2.6.0", "a" * 40)

    def test_unsupported_architecture_fails(self):
        for suffix in ("", ".spdx.json"):
            (self.dist / (self.names[0] + suffix)).rename(
                self.dist / (self.names[0].replace("x64", "unknown") + suffix))
        self.assertNotEqual(self.prepare().returncode, 0)


if __name__ == "__main__":
    unittest.main()
