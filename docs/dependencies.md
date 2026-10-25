<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
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

Normal clones need no submodule commands. A first configure from a Git clone
needs HTTPS access to codeload.github.com. Release source bundles instead use
the included original archive in `third_party/`, verified against the same hash. For offline packaging, provide an extracted archive or
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
   Refresh `LICENSES/LicenseRef-PFFFT.txt` from the dependency's LICENSE.txt
   and CMake-generated PFFFT-NOTICES.txt, preserving their full text. Update
   its notice hash and review record in `legal/components.json`; regenerate
   the component table with `python3 tools/legal_inventory.py --write-table`.
3. Configure with downloads enabled and no local override. Run Debug scalar
   and SIMD, optimized GCC/Clang, sanitizer and NOTONES suites, allocation
   failure tests, quality preflight and OSV. Compare DVSI scores against a
   baseline before accepting numerical changes or refreshing golden hashes.
4. To roll back, restore **both** pin lines from the previous revision, or run
   the same command with the previous full SHA; reconfigure and repeat checks.

ClusterFuzzLite checks out the same SHA with `tools/fetch-pinned-git.sh` during
image preparation. `tools/pffft_metadata.py` generates OSV's explicit commit
input from the pin. `tools/release_sbom.py` uses `legal/components.json` to
record embedded PFFFT with its revision, archive checksum and exact retained
`LicenseRef-PFFFT` text before upload and attestation. Source-only OP25 tooling
and unresolved ETSI license conclusions are recorded separately. Pin and integration edits,
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

Semgrep 1.182.1 pins two vulnerable dependency ranges:

- `mcp==1.29.0`, which is affected by
  [GHSA-84m7-p3x7-pcfv](https://github.com/modelcontextprotocol/python-sdk/security/advisories/GHSA-84m7-p3x7-pcfv).
  The Semgrep requirements override it with `mcp==1.30.0`.
- `opentelemetry-*~=1.37.0` and `~=0.58b0`, which hold `opentelemetry-api` below
  the fix for [GHSA-xcq7-p3pf-p88h](https://github.com/advisories/GHSA-xcq7-p3pf-p88h).
  The Semgrep requirements override the OpenTelemetry packages Semgrep names
  with 1.45.1 and 0.66b1, so the family stays in step.

Regenerate the complete hashed lockfile with the command in its header, including
`--override .github/requirements/semgrep.in`. Install it with
`python -m pip install --require-hashes --no-deps -r .github/requirements/semgrep.txt`
in a dedicated Python 3.14 environment. CI checks that `pip check` reports only
the exact Semgrep/MCP and Semgrep/OpenTelemetry metadata conflicts, then runs
the Semgrep guardrails. When a Semgrep release accepts a patched version, drop
that override and its expected conflict lines. When none remain, also remove
`--no-deps` and the expected-conflict check; regenerate the lockfile and require
an ordinary `pip check` to pass. Bumping the Semgrep pin changes the version in
the expected conflict lines too. This affects CI tooling only.

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

### PFFFT update to 548918d

The update from `77a7ecd` to `548918d` adds the fork's 8-lane AVX2 float
backend. On x86-64, `pffft_new_setup()` uses it only on CPUs it was measured
faster on (AMD Zen 3 and Zen 4) and only at the sizes listed for that CPU in
`src/pffft_wide_sizes.h`, after a CPUID/XGETBV check; everything else keeps
the 4-lane code. mbelib's only transform, the real 256-point FFT, is in none
of those tables, so its FFT code path is unchanged on every CPU. GCC and
Clang build the backend inside `pffft.c` with target attributes (MinGW GCC
only in optimized builds). MSVC builds it as a translation unit of its own,
`src/pffft_wide_tu.c`, which `CMakeLists.txt` now adds for MSVC x86-64 with
`/arch:AVX2 /GL-` and `PFFFT_WIDE_TU`. The fork's 4-lane complex transforms
also fuse their deinterleave/interleave into the first and last passes.

The new compiled headers (`pffft_passes_impl.h`, `pffft_internal.h`,
`pffft_setup.h`, `pffft_cpu_x86.h`, `pffft_wide_impl.h`,
`pffft_wide_select.h`, `pffft_wide_sizes.h`, `pffft_wide_api.h`) and
`pffft_wide_tu.c` were added to the compiled-source notices; the retained
PFFFT text and its inventory hash were refreshed. Both pins were regenerated
from the immutable archive.

Every golden hash is unchanged. Local suites passed for Release GCC with
warnings as errors, Release Clang, Debug SIMD (SSE2 and AVX2 golden PCM),
Debug scalar, ASan/UBSan and NOTONES, and the acquisition regressions. MSVC
is covered by the existing Windows CI.

### PFFFT update to 77a7ecd

The update from `0ac3d39` to `77a7ecd` changes one compiled header,
`src/pffft_real_impl.h`: the ordered real transform's copy of the high float
pairs of two SIMD vectors now addresses their bytes (`unsigned char`
pointers and a byte offset) instead of float pointers into the vector
objects, which CodeQL reported. The copied bytes are the same, so the output
is unchanged. The other changes are the fork's own CMake options (host-tuned
defaults for its top-level builds, which never reach this library: it compiles
PFFFT from its own source list and flags), CI, README and tests. Both pins
were regenerated from the immutable archive; the compiled-source notices are
unchanged, so only the pinned revision in the retained PFFFT text and its
inventory hash changed.

Every Debug golden hash is unchanged: `gen_golden` for scalar, SSE2 and AVX2
(the dispatched copy and a `-march=x86-64-v3` build) matches the values in
`tests/test_golden_pcm.c`.

### PFFFT update to 0ac3d39

The update from `94d3c83` to `0ac3d39` includes upstream FFT accuracy and
allocation fixes, more accurate twiddle factors, and fused ordered-transform
layout passes. Both pins were regenerated from the immutable archive. The
compiled-source notices now include the new complex and real implementation
headers, and the retained PFFFT text and inventory hash were refreshed.

The twiddle changes alter the six AMBE float golden hashes in Debug scalar,
SSE2 and AVX2 builds. All int16 golden hashes and the single-frame float hashes
stay unchanged. The replacement hashes agree between GCC 16 on the development
host and GCC 13 with Ubuntu 24.04's glibc 2.39.

Baseline and candidate Debug scalar/SSE2 builds were compared on 75 DVSI
development/validation vectors each. Every alignment
was unchanged. Out of 8,847,360 decoded samples per codec and configuration,
scalar changed 81 D-STAR, 111 IMBE and 100 AMBE+2 samples; SSE2 changed 85, 107
and 91 respectively. Every changed sample differs by one int16 unit. The
maximum absolute per-vector log-spectral-distance change was 0.000015 dB
against DVSI and 0.000005 dB against the input. Encoder comparisons and the
final holdout partition were not rerun for this dependency update. Vectors,
audio and detailed reports remain local under `build/pin-quality/`.

Local suites passed for Debug scalar/SSE2/AVX2, optimized GCC/Clang with
fast-math and LTO (Clang uses LLD), ASan/UBSan and NOTONES. The acquisition
regressions, full quality preflight, license inventory, OSV and secret scans
also passed.
