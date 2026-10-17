#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Exercise licensing evidence, exact source delivery and conservative data handling."""

import copy
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools/quality"))
import legal_inventory as inventory
import release_assets
import release_sbom
import source_bundle
import reference_provenance as provenance


class Evidence(unittest.TestCase):
    def test_inventory_notices_exceptions_and_coverage(self):
        data = inventory.load()
        read = lambda name: (ROOT / name).read_bytes()
        inventory.validate(data, read, ["src/core/mbelib.c", "src/tetra/tetra_acelp_codebooks.c"])
        changed = copy.deepcopy(data)
        next(c for c in changed["components"] if c["id"] == "etsi-data").pop("exception")
        with self.assertRaisesRegex(ValueError, "exception"):
            inventory.validate(changed, read, [])
        changed = copy.deepcopy(data)
        changed["components"] = [c for c in changed["components"] if c["id"] != "isc-origins"]
        with self.assertRaisesRegex(ValueError, "incomplete"):
            inventory.validate(changed, read, [])
        with self.assertRaises((ValueError, KeyError, FileNotFoundError)):
            inventory.validate(data, lambda n: b"" if n == "LICENSE" else {}[n], [])
        with self.assertRaisesRegex(ValueError, "Unrecorded source license"):
            inventory.validate(data, lambda n: b"// SPDX-License-Identifier: GPL-3.0-only\n" if n == "src/new.c" else read(n), ["src/new.c"])
        with self.assertRaisesRegex(ValueError, "Unrecorded source license"):
            inventory.validate(data, lambda n: b"// SPDX-License-Identifier: GPL-2.0-or-later AND LicenseRef-Unknown\n" if n == "src/new.c" else read(n), ["src/new.c"])

    def test_declared_and_concluded_licenses_and_relationships(self):
        for kind in ("source", "binary"):
            doc = {"packages": [], "relationships": [], "hasExtractedLicensingInfos": [], "files": []}
            release_sbom.add_components(doc, inventory.load(), source_bundle.read_pins(), kind,
                                        "a" * 40, lambda n: (ROOT / n).read_bytes())
            indexed = {p["SPDXID"]: p for p in doc["packages"]}
            self.assertEqual(indexed["SPDXRef-Component-etsi-data"]["licenseConcluded"], "NOASSERTION")
            self.assertEqual(indexed["SPDXRef-Component-pffft"]["licenseDeclared"], "LicenseRef-PFFFT")
            self.assertEqual("SPDXRef-Component-op25-tooling" in indexed, kind == "source")
            self.assertNotIn("SPDXRef-Component-dvsi-reference", indexed)
            extracted = {e["licenseId"]: e["extractedText"] for e in doc["hasExtractedLicensingInfos"]}
            self.assertIn("Copyright", extracted["LicenseRef-PFFFT"])
            self.assertIn("unresolved", extracted["LicenseRef-ETSI-EN-300-395-2-data"])
            relation = next(r for r in doc["relationships"] if r["relatedSpdxElement"] == "SPDXRef-Component-pffft")
            self.assertEqual(relation["relationshipType"], "STATIC_LINK" if kind == "binary" else "CONTAINS")

    def test_binary_notices_and_unexpected_data(self):
        data = inventory.load()
        files = {"share/doc/mbelib-neo/" + n: (ROOT / n).read_bytes()
                 for c in data["components"] if "binary" in c["distribution"] for n in c["notices"]}
        retained = (ROOT / "LICENSES/LicenseRef-PFFFT.txt").read_bytes()
        files["share/doc/mbelib-neo/pffft/PFFFT-NOTICES.txt"] = retained[retained.index(b"PFFFT: https://"):]
        release_sbom.validate_package(files, "binary")
        # Git and installation tools may use CRLF on Windows; notice text is retained.
        release_sbom.validate_package(
            {n: content.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n") for n, content in files.items()}, "binary")
        for n in ("LICENSES/ISC.txt", "LICENSES/LicenseRef-PFFFT.txt", "LICENSES/LicenseRef-ETSI-EN-300-395-2-data.txt"):
            changed = dict(files)
            changed.pop("share/doc/mbelib-neo/" + n)
            with self.assertRaisesRegex(ValueError, "notice"):
                release_sbom.validate_package(changed, "binary")
        with self.assertRaisesRegex(ValueError, "Unexpected"):
            release_sbom.validate_package({**files, "share/audio.wav": b"private"}, "binary")
        files["share/doc/mbelib-neo/pffft/PFFFT-NOTICES.txt"] = b"incomplete"
        with self.assertRaisesRegex(ValueError, "PFFFT notices"):
            release_sbom.validate_package(files, "binary")

    def test_versioned_release_notes(self):
        self.assertNotIn("source.tar.gz", release_assets.release_notes("2.5.1"))
        current = release_assets.release_notes("2.6.0")
        self.assertIn("2.6.0-source.tar.gz", current)
        self.assertIn("v2.6.0/legal/status.md", current)
        self.assertIn("v2.6.0/legal/patents.md", current)


@unittest.skipUnless(shutil.which("git"), "source builder tests require Git")
class SourceBundle(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.repo = self.root / "repo"
        self.repo.mkdir()
        self.version = inventory.load()["release"]
        self.tag = "v" + self.version
        required = json.loads((ROOT / "legal/source-inventory.json").read_text())["required"]
        required += ["legal/source-inventory.json", "docs/dependencies.md", "docs/tetra-acelp.md",
                     "tools/quality/build_op25_encoder.sh", "tools/quality/fetch_speech_corpus.py"]
        for name in set(required):
            target = self.repo / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / name, target)
        self.dependency = self.root / "pffft.tar.gz"
        # A synthetic dependency fixture; no external reference data or network.
        self.dependency.write_bytes(b"synthetic pinned archive")
        self.pins = {"PFFFT_SHA": "1" * 40, "PFFFT_ARCHIVE_SHA256": hashlib.sha256(self.dependency.read_bytes()).hexdigest()}
        (self.repo / "tools/ci-dependency-pins.env").write_text("\n".join(f"{k}={v}" for k, v in self.pins.items()))
        self.command("init", "-q")
        self.command("config", "core.hooksPath", str(self.root / "no-hooks"))
        self.command("config", "user.name", "Source test")
        self.command("config", "user.email", "source-test@example.invalid")
        self.command("add", ".")
        self.command("-c", "commit.gpgsign=false", "commit", "-qm", "test fixture")
        self.command("-c", "tag.gpgsign=false", "tag", self.tag)
        self.commit = self.command("rev-parse", "HEAD").strip()

    def command(self, *args):
        return subprocess.check_output(["git", "-C", str(self.repo), *args], text=True)

    def bundle(self, output="dist"):
        return source_bundle.create(self.repo, self.tag, self.commit, self.root / output, self.dependency)

    def test_exact_reproducible_source_and_identity(self):
        original = source_bundle.git(self.repo, "cat-file", "blob", f"{self.commit}:README.md")
        (self.repo / "README.md").write_text("dirty working copy")
        (self.repo / "private-audit.txt").write_text("must stay local")
        first = self.bundle()
        self.command("config", "tar.umask", "0077")
        self.command("config", "core.autocrlf", "true")
        self.command("config", "core.eol", "crlf")
        second = self.bundle("other")
        self.assertEqual(first.read_bytes(), second.read_bytes())
        files = release_sbom.package_files(first)
        self.assertEqual(files["README.md"], original)
        self.assertNotIn("private-audit.txt", files)
        _, _, identity, _ = release_sbom.validate_package(files, "source")
        self.assertEqual(identity["commit"], self.commit)
        sbom = release_sbom.create(first, "source", self.version, self.commit, identity["source_date_epoch"])
        self.assertEqual(sbom["packages"][0]["licenseDeclared"], "NOASSERTION")
        self.assertEqual(len(sbom["files"]), len(files))
        files["README.md"] = b"tampered"
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            release_sbom.validate_package(files, "source")

    def test_wrong_dependency_hash_and_commit_fail(self):
        self.dependency.write_bytes(b"wrong")
        with self.assertRaisesRegex(ValueError, "SHA256"):
            self.bundle()
        with self.assertRaisesRegex(ValueError, "validated commit"):
            source_bundle.create(self.repo, self.tag, "0" * 40, self.root, self.dependency)

    def test_archive_attributes_cannot_change_tag_contents(self):
        (self.repo / "README.md").write_text("Release commit: $Format:%H$\n")
        (self.repo / ".gitattributes").write_text("README.md export-subst\n")
        self.command("add", ".")
        self.command("-c", "commit.gpgsign=false", "commit", "-qm", "archive substitution")
        with self.assertRaisesRegex(ValueError, "changed tracked source"):
            source_bundle.tracked_source(self.repo, self.command("rev-parse", "HEAD").strip())
        (self.repo / ".gitattributes").write_text("docs/dependencies.md export-ignore\n")
        self.command("add", ".")
        self.command("-c", "commit.gpgsign=false", "commit", "-qm", "archive omission")
        with self.assertRaisesRegex(ValueError, "every tracked source file"):
            source_bundle.tracked_source(self.repo, self.command("rev-parse", "HEAD").strip())

    def test_missing_notice_and_forbidden_tracked_data_fail(self):
        (self.repo / "LICENSES/ISC.txt").unlink()
        self.command("add", "-u")
        self.command("-c", "commit.gpgsign=false", "commit", "-qm", "missing notice")
        with self.assertRaisesRegex(ValueError, "Missing required"):
            source_bundle.tracked_source(self.repo, self.command("rev-parse", "HEAD").strip())
        shutil.copyfile(ROOT / "LICENSES/ISC.txt", self.repo / "LICENSES/ISC.txt")
        (self.repo / "tests/private.wav").parent.mkdir(exist_ok=True)
        (self.repo / "tests/private.wav").write_bytes(b"restricted recording")
        self.command("add", ".")
        self.command("-c", "commit.gpgsign=false", "commit", "-qm", "bad data")
        with self.assertRaisesRegex(ValueError, "excludes"):
            source_bundle.tracked_source(self.repo, self.command("rev-parse", "HEAD").strip())

    def test_unsafe_paths_rejected(self):
        staging = json.loads((ROOT / "legal/source-inventory.json").read_text())
        for name in (".git/config", "build/audio.raw", "docs/audit.wav", "../LICENSE", "/LICENSE", "private-audit/foo"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                source_bundle.check_path(name, staging)


class ReferenceData(unittest.TestCase):
    def test_unknown_rights_stay_local_and_notices_survive(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            notice = root / "COPYING"
            notice.write_text("Synthetic fixture permission notice")
            sha = provenance.digest(notice)
            evidence = provenance.acquisition("https://example.invalid/audio", "https://example.invalid/terms",
                                              notices=[{"sha256": sha}])
            manifest = {"references": [{"name": "synthetic", "source": {"acquisition": evidence}}],
                        "corpus": {"notices": [{"path": "COPYING", "sha256": sha}]}}
            presentation = root / "pack"
            provenance.carry_notices(root, manifest, presentation)
            self.assertEqual((presentation / "notices/0-COPYING").read_bytes(), notice.read_bytes())
            self.assertIsNone(evidence["retrieved_at"])
            with self.assertRaisesRegex(ValueError, "permission"):
                provenance.check_public(presentation)
            evidence.update(redistribution="allowed", permission_source="synthetic test author grant")
            provenance.carry_notices(root, manifest, presentation)
            provenance.check_public(presentation)
            (presentation / "notices/0-COPYING").write_text("altered")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                provenance.check_public(presentation)

    def test_no_manifest_or_restricted_data_is_not_public(self):
        self.assertEqual(provenance.distribution({"references": [{"name": "unknown"}]})["redistribution"], "local-only")
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                provenance.safe_path(Path(directory), "../secret")

    def test_older_local_run_retains_notices(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "COPYING").write_text("Synthetic permission notice")
            notices = [{"file": "COPYING", "sha256": provenance.digest(root / "COPYING")}]
            (root / "corpus-manifest.json").write_text(json.dumps({"notices": notices}))
            manifest = {"references": [{"name": "old-run"}], "corpus": {"manifest": "corpus-manifest.json"}}
            provenance.carry_notices(root, manifest, root / "pack")
            self.assertEqual((root / "pack/notices/0-COPYING").read_bytes(), (root / "COPYING").read_bytes())
            with self.assertRaisesRegex(ValueError, "permission"):
                provenance.check_public(root / "pack")

    def test_dataset_report_keeps_evidence_and_checks_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "dataset.json"
            self.assertEqual(provenance.retain_dataset(path, root / "out")["redistribution"], "unknown")
            (root / "NOTICE").write_text("Synthetic reference notice")
            evidence = provenance.acquisition("https://example.invalid/source", "https://example.invalid/terms",
                                              status="restricted", notices=[{"file": "NOTICE", "sha256": provenance.digest(root / "NOTICE")}])
            path.write_text(json.dumps({"acquisition": evidence}))
            record = provenance.retain_dataset(path, root / "out")
            self.assertEqual(record["redistribution"], "restricted")
            self.assertEqual((root / "out/notices/0-NOTICE").read_bytes(), (root / "NOTICE").read_bytes())
            (root / "NOTICE").write_text("tampered")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                provenance.retain_dataset(path, root / "out")


if __name__ == "__main__":
    unittest.main()
