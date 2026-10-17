<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# Building a corresponding-source archive offline

Starting with 2.6.0, download `mbelib-neo-<version>-source.tar.gz` and its SBOM
from the same release as the binary packages. Verify checksums and attestations
as described in [release-verification.md](release-verification.md), then extract
the archive. `SOURCE-IDENTITY.json` identifies the validated tag, commit,
dependency pin and file hashes. Git metadata is unnecessary.

With a C99 compiler, CMake >=3.20, a build tool and Python 3 available locally:

```sh
cmake --preset dev-debug
cmake --build --preset dev-debug --parallel
ctest --preset dev-debug --output-on-failure
cmake --install build/dev-debug --prefix /your/install/prefix
```

The original pinned PFFFT archive is in `third_party/`; CMake verifies its
SHA-256 before extracting it locally. No download or populated build cache is
needed. An explicit `FETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT` override takes
precedence and makes its provenance the caller's responsibility. Ordinary Git
clones without the bundle retain the pinned network-download fallback.

Optional speech measurements, corpus acquisition, analyzers requiring their
own tools, and external conformance vectors are not prerequisites for this
offline build. They are not bundled. Public tests use synthetic fixtures.
Read [licensing status](../legal/status.md) before redistributing the result;
complete source delivery does not resolve the disclosed ETSI or patent issues.
