#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Recover publication only from a completed, validated tag build and its attestations."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess

from release_assets import collect, needs_source

REPO = "arancormonk/mbelib-neo"
WORKFLOW = ".github/workflows/ci.yml"
SPDX = "https://spdx.dev/Document"
REQUIRED_JOBS = {"Validate release tag", *(f"package release ({platform}-latest)"
                                         for platform in ("ubuntu", "macos", "windows"))}
TAG_CHECKS = {
    "Validate release tag", "clang-format", "clang-tidy", "cppcheck", "lizard", "scan-build", "semgrep",
    *(f"{platform}-latest • preset={preset}" for platform in ("ubuntu", "macos", "windows")
      for preset in ("dev-debug", "dev-release", "notones-debug")),
    "linux • Clang • Release", "linux • Debug • sanitizers",
    *(f"linux • ARM32 Release • FPU={fpu} • QEMU tests" for fpu in ("default", "neon")),
    "linux • GCC LTO link (x86-64-v4)", "linux • Release • ctest (x86-64-v3, fast-math)",
    "linux • Release • ctest (x86-64-v3, fast-math, clang)",
    "linux • Release • ctest (consumer defaults, GCC)",
    "install + consume (CMake, Windows)",
    *(f"install + consume (pkg-config/CMake) ({platform}-latest)" for platform in ("ubuntu", "macos")),
}


def gh_json(*arguments):
    return json.loads(subprocess.check_output(["gh", *arguments], text=True))


def validate_run(run, jobs, tag, commit, rebuild=False):
    if rebuild and not needs_source(tag[1:]):
        raise ValueError("Package rebuild requires the component-aware release format from 2.6.0")
    expected = {"event": "push", "head_branch": tag, "head_sha": commit,
                "path": WORKFLOW, "status": "completed"}
    if any(run.get(key) != value for key, value in expected.items()):
        raise ValueError("Source run must be a completed CI push build of the signed tag commit")
    successful = {job["name"] for job in jobs if job.get("conclusion") == "success"}
    if rebuild and not TAG_CHECKS.issubset(successful):
        raise ValueError("Source tag run lacks successful analysis, native, optimized, sanitizer or consumer checks")
    if not rebuild and not REQUIRED_JOBS.issubset(successful):
        raise ValueError("Source run lacks successful tag validation or all three attested packages")
    if any(job.get("conclusion") not in ("success", "skipped")
           for job in jobs if job["name"] != "publish release"
           and not (rebuild and job["name"].startswith("package release ("))):
        raise ValueError("Only the final publish job may have failed in the source run")


def rebuild_identity(run_id):
    if not re.fullmatch(r"[0-9]+", run_id):
        raise ValueError("Expected a numeric rebuild run ID")
    run = gh_json("api", f"repos/{REPO}/actions/runs/{run_id}")
    if (run.get("event") != "workflow_dispatch" or run.get("head_branch") != "main"
            or run.get("path") != WORKFLOW or not re.fullmatch(r"[0-9a-f]{40}", run.get("head_sha", ""))):
        raise ValueError("Package rebuild must use the reviewed CI workflow on main")
    comparison = gh_json("api", f"repos/{REPO}/compare/{run['head_sha']}...main")
    if comparison.get("status") not in ("identical", "ahead"):
        raise ValueError("Rebuild workflow commit is not contained in main")
    return "refs/heads/main", run["head_sha"]


def verify_assets(directory, tag, commit, builder=None):
    packages, _ = collect(directory, tag[1:])
    build_ref, build_commit = builder or (f"refs/tags/{tag}", commit)
    for package in packages:
        arguments = ("attestation", "verify", str(package), "--repo", REPO,
                     "--signer-workflow", f"{REPO}/{WORKFLOW}", "--source-ref", build_ref,
                     "--source-digest", build_commit, "--deny-self-hosted-runners", "--format=json")
        gh_json(*arguments)
        verified = gh_json(*arguments, "--predicate-type", SPDX)
        sbom = json.loads(Path(str(package) + ".spdx.json").read_text(encoding="utf-8"))
        if not any(item["verificationResult"]["statement"]["predicate"] == sbom for item in verified):
            raise ValueError(f"Downloaded SBOM differs from the signed predicate for {package.name}")
        if builder and not any(p.get("SPDXID") == "SPDXRef-Package-mbelib-neo"
                               and p.get("versionInfo") == tag[1:]
                               and p.get("sourceInfo") == f"Release commit {commit}; see legal/status.md"
                               for p in sbom.get("packages", [])):
            raise ValueError("Signed rebuild SBOM does not identify the validated release source")
        if needs_source(tag[1:]):
            # New releases also attest the SBOM file as a provenance subject.
            gh_json("attestation", "verify", str(package) + ".spdx.json", "--repo", REPO,
                    "--signer-workflow", f"{REPO}/{WORKFLOW}", "--source-ref", build_ref,
                    "--source-digest", build_commit, "--deny-self-hosted-runners", "--format=json")
        print(f"Verified original tag provenance and SBOM: {package.name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets", type=Path)
    parser.add_argument("--rebuild", action="store_true", help="Require original tag checks before rebuilding packages")
    parser.add_argument("--builder-run-id", help="Verify a reviewed main-workflow rebuild and its signed source identity")
    args = parser.parse_args()
    tag = os.environ.get("RELEASE_TAG", "")
    commit = os.environ.get("RELEASE_COMMIT", "")
    run_id = os.environ.get("RELEASE_RECOVERY_RUN_ID", "")
    if not re.fullmatch(r"v\d+\.\d+\.\d+", tag) or not re.fullmatch(r"[0-9a-f]{40}", commit):
        parser.error("Expected a validated release tag and full commit SHA")
    if args.assets:
        verify_assets(args.assets, tag, commit, rebuild_identity(args.builder_run_id) if args.builder_run_id else None)
    else:
        if not re.fullmatch(r"[0-9]+", run_id):
            parser.error("Expected a numeric source run ID")
        endpoint = f"repos/{REPO}/actions/runs/{run_id}"
        run = gh_json("api", endpoint)
        pages = gh_json("api", f"{endpoint}/jobs?per_page=100", "--paginate", "--slurp")
        validate_run(run, [job for page in pages for job in page["jobs"]], tag, commit, args.rebuild)
        print(f"Validated original release build {run_id} for {tag}")


if __name__ == "__main__":
    main()
