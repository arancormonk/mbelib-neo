# Dependency Management

mbelib-neo keeps runtime dependencies small and reviewable.

## Compiled Dependencies

The project compiles `pffft.c` and `pffft_common.c` from
<https://github.com/arancormonk/pffft> into both mbelib libraries. The fork's
separate FFTPACK reference implementation is not needed. `cmake/Pffft.cmake`
populates a commit archive with CMake FetchContent, SHA-256 verification and
TLS verification, without configuring the dependency's C++ project.
`PFFFT_SHA` and `PFFFT_ARCHIVE_SHA256` in `tools/ci-dependency-pins.env` are the
source of truth. Builds never follow a branch or tag.

Normal clones need no submodule commands. A first configure needs HTTPS access
to codeload.github.com. For offline packaging, provide an extracted archive or
checkout with `-DFETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT=/absolute/path/to/pffft`.
This explicit override bypasses download/hash checking; the packager verifies
its provenance and must use the compatible fork. It needs no Git metadata.
For an already populated build, `-DFETCHCONTENT_FULLY_DISCONNECTED=ON` avoids
network access. Disable it before changing the pin so CMake can fetch the update.

The fork retains upstream license headers. Installation includes LICENSE.txt
and the compiled-source notices under the project's documentation `pffft/`
directory. No public PFFFT headers or external link requirements are installed.

### Manual update and rollback

1. Land focused fixes and tests in the fork; do not patch downloaded sources.
2. Run `python3 tools/update_pffft.py FULL_40_CHARACTER_COMMIT_SHA`. It downloads
   the commit archive, checks the production layout, computes SHA-256 and
   replaces both pins together. Review the fork diff and both pin changes.
3. Configure with downloads enabled and no local override. Run Debug scalar
   and SIMD, optimized GCC/Clang, sanitizer and NOTONES suites, allocation
   failure tests, quality preflight and OSV. Compare DVSI scores against a
   baseline before accepting numerical changes or refreshing golden hashes.
4. To roll back, restore **both** pin lines from the previous revision, or run
   the same command with the previous full SHA; reconfigure and repeat checks.

ClusterFuzzLite checks out the same SHA with `tools/fetch-pinned-git.sh` during
image preparation. `tools/pffft_metadata.py` generates OSV's explicit commit
input from the pin and enriches release SPDX SBOMs with the fork URL, revision
and archive checksum before upload and attestation. Pin and integration edits,
including deletions, trigger PR fuzzing and dependency scanning.

System dependencies used by consumers are:

- C99 compiler
- CMake 3.20 or newer
- platform math library where required by the toolchain

## Tooling Dependencies

CI and local quality tools are tracked through:

- `.github/requirements/*.in`
- `.github/requirements/*.txt`
- `.github/dependabot.yml`
- `.github/workflows/*.yml`
- `tools/*.sh`

Hashed Python requirements are used where Python tooling is installed in CI.
GitHub Actions are pinned to immutable commit SHAs by policy.

Semgrep 1.179.0 pins `mcp==1.29.0`, which is affected by
[GHSA-84m7-p3x7-pcfv](https://github.com/modelcontextprotocol/python-sdk/security/advisories/GHSA-84m7-p3x7-pcfv).
The Semgrep requirements explicitly override that pin with `mcp==1.30.0`.
Regenerate the complete hashed lockfile with the command in its header, including
`--override .github/requirements/semgrep.in`. Install it with
`python -m pip install --require-hashes --no-deps -r .github/requirements/semgrep.txt`
in a dedicated Python 3.14 environment. CI checks that `pip check` reports only
the exact Semgrep/MCP metadata conflict, then runs the Semgrep guardrails.
Remove the override, `--no-deps`, and the expected-conflict check together when
a Semgrep release accepts the patched SDK; regenerate the lockfile and require
an ordinary `pip check` to pass. This affects CI tooling only.

## Monitoring

The project monitors dependencies through:

- Dependabot for GitHub Actions and Python requirements
- GitHub dependency review on pull requests
- OSV-Scanner through `tools/osv_scan.sh`
- scheduled guardrail CI
- OpenSSF Scorecard

If a vulnerability is reported in a dependency:

1. Determine whether the vulnerable code is present and reachable.
2. Update or patch the dependency when exploitable.
3. If not exploitable, document the reason in the narrowest available
   suppression or advisory note.
4. Add tests or checks when the issue could regress.

## Update Policy

Dependency updates should:

- preserve license notices
- avoid unrelated refactors
- update documentation when dependency behavior or requirements change
- run focused codec tests and `tools/osv_scan.sh`
- receive human review when they affect compiled code, workflows, packaging, or
  release behavior

### AUR packaging

PFFFT remains embedded in mbelib. Neither `mbelib-neo` nor its consumer
`dsd-neo` needs a separate PFFFT `depends` or `makedepends` entry. The mbelib
release and `-git` AUR workflows run `tools/update_aur_pffft.py` to append the
pinned archive to `source`, its checksum to `sha256sums`, and the local source
override plus disconnected mode to the CMake invocation. makepkg fetches and
verifies the archive before `build()`. The recipe rejects a source pin mismatch
(for example, a newer `-git` checkout with an older AUR recipe).

To prepare an existing recipe manually:

```sh
python3 tools/update_aur_pffft.py /path/to/aur/PKGBUILD
# Run in the AUR package directory:
makepkg --printsrcinfo > .SRCINFO
```

Review and publish the AUR change with the corresponding mbelib update. The
helper preserves the original source checksum and package dependencies.

### PFFFT migration measurements

The Debug scalar baseline at mbelib commit `c2d0b801238cd6b23e0d533e1dd8453362c61e88`
was compared with the pinned fork across all 237 locally available DVSI vectors
(development, validation, final and tones). D-STAR and AMBE+2 encoder rows were
identical. IMBE changed one bit in one of 3,073 frames of `irstia`.
Decoder PCM changed by at most one int16 unit: 96 D-STAR, 197 IMBE and 96
AMBE+2 samples, out of 13,678,400 samples per codec. The maximum per-vector
absolute change in decoder spectral error was 7.16e-7 dB; encoder spectral
error changed by at most 2.63e-6 dB and encoder level by 6.46e-6 dB.

The AMBE+2 and D-STAR encoder tone gates, call-progress timing and tone-edge
gates passed. Historical pitch/shelf/low-band/encoder gates require an
improvement or bit-identical PCM and do not pass this FFT replacement; their
thresholds were not changed. The legacy acoustic tone gate also fails its
absolute frequency/agreement thresholds on the baseline. Scalar and SIMD
AMBE float golden hashes were refreshed after these comparisons and the FFT
correctness tests; every int16 golden hash stayed unchanged. Vectors and
per-frame work files remain local under `build/pffft-migration/`.

Acquisition regression checks are opt-in because they access the network:
`python3 tests/test_pffft_acquisition.py [--cmake /path/to/cmake]` covers fresh
fetches, cached reconfiguration, offline overrides with spaces in their paths,
source distributions without Git metadata, bad paths/layouts, and bad hashes.

Local migration validation passed with CMake 3.20.5 and 4.4.4, GCC 16 and
Clang 23: Debug scalar/SIMD (including x86-64-v3), optimized Release with
fast-math/LTO, ASan/UBSan, NOTONES, and 32-bit x87. Clang LTO used LLD.
The full quality preflight (including scan-build, OSV and workflow checks),
release hardening, separate installed static/shared consumers, and all three
fuzz harnesses (500 sanitizer runs each) passed. The shared library's 108
exports are unchanged. A generated AUR recipe built with the local override
and disconnected mode, and both release/VCS recipes generated the expected
`.SRCINFO` with makepkg. Windows, macOS and ARM execution remains covered by
the existing CI matrix and was not run on this Linux x86 host.
