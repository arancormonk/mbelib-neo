#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Recover publication only from a completed, validated tag build and its attestations."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess

REPO = "arancormonk/mbelib-neo"
WORKFLOW = ".github/workflows/ci.yml"
SPDX = "https://spdx.dev/Document"
REQUIRED_JOBS = {"Validate release tag", *(f"package release ({platform}-latest)"
                                         for platform in ("ubuntu", "macos", "windows"))}


def gh_json(*arguments):
    return json.loads(subprocess.check_output(["gh", *arguments], text=True))


def validate_run(run, jobs, tag, commit):
    expected = {"event": "push", "head_branch": tag, "head_sha": commit,
                "path": WORKFLOW, "status": "completed"}
    if any(run.get(key) != value for key, value in expected.items()):
        raise ValueError("Source run must be a completed CI push build of the signed tag commit")
    successful = {job["name"] for job in jobs if job.get("conclusion") == "success"}
    if not REQUIRED_JOBS.issubset(successful):
        raise ValueError("Source run lacks successful tag validation or all three attested packages")
    if any(job.get("conclusion") not in ("success", "skipped")
           for job in jobs if job["name"] != "publish release"):
        raise ValueError("Only the final publish job may have failed in the source run")


def verify_assets(directory, tag, commit):
    packages = sorted([*directory.glob("*.tar.gz"), *directory.glob("*.zip")])
    if len(packages) != 3:
        raise ValueError("Expected three release packages")
    for package in packages:
        arguments = ("attestation", "verify", str(package), "--repo", REPO,
                     "--signer-workflow", f"{REPO}/{WORKFLOW}", "--source-ref", f"refs/tags/{tag}",
                     "--source-digest", commit, "--deny-self-hosted-runners", "--format=json")
        gh_json(*arguments)
        verified = gh_json(*arguments, "--predicate-type", SPDX)
        sbom = json.loads(Path(str(package) + ".spdx.json").read_text(encoding="utf-8"))
        if not any(item["verificationResult"]["statement"]["predicate"] == sbom for item in verified):
            raise ValueError(f"Downloaded SBOM differs from the signed predicate for {package.name}")
        print(f"Verified original tag provenance and SBOM: {package.name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets", type=Path)
    args = parser.parse_args()
    tag = os.environ.get("RELEASE_TAG", "")
    commit = os.environ.get("RELEASE_COMMIT", "")
    run_id = os.environ.get("RELEASE_RECOVERY_RUN_ID", "")
    if not re.fullmatch(r"v\d+\.\d+\.\d+", tag) or not re.fullmatch(r"[0-9a-f]{40}", commit):
        parser.error("Expected a validated release tag and full commit SHA")
    if args.assets:
        verify_assets(args.assets, tag, commit)
    else:
        if not re.fullmatch(r"[0-9]+", run_id):
            parser.error("Expected a numeric source run ID")
        endpoint = f"repos/{REPO}/actions/runs/{run_id}"
        run = gh_json("api", endpoint)
        pages = gh_json("api", f"{endpoint}/jobs?per_page=100", "--paginate", "--slurp")
        validate_run(run, [job for page in pages for job in page["jobs"]], tag, commit)
        print(f"Validated original release build {run_id} for {tag}")


if __name__ == "__main__":
    main()
