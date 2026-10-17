<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# Reference data provenance and distribution

Reviewed 2026-10-17. DVSI materials and ETSI reference archives, vectors,
reference builds and generated audio stay outside release artifacts. Local
measurement is an intended use, not a finding of permission. Public tests and
examples use synthetic or specifically licensed fixtures. Restricted-reference
comparisons remain optional and local.

Acquisition manifests record source and terms URLs, retrieval date (null when
not known), observation date, content hashes, notice hashes and redistribution
status. An old cache without retrieval evidence stays unknown; recording it
today does not invent a historical retrieval date. Use
`tools/quality/reference_provenance.py` to record user-acquired ETSI archives
and vectors without opening reference program code. Supply the actual retrieval
date only when known. Neither a hash nor a downloaded notice establishes rights.

DVSI vectors are marked `restricted` for project distribution. ETSI reference
materials and legacy Codec2 speech samples remain `unknown` unless their own
applicable terms establish otherwise. Codec2's software license is not presumed
to license every contributed recording. CMU ARCTIC acquisition retains each
speaker's COPYING file and its hash; distribution still requires review of the
specific recording and notice. An absent record always means unknown.

The quality runner copies source manifests and notices into its immutable run.
Reports carry distribution metadata, and listening presentations include a
`DISTRIBUTION.json` and the original notices. Local analysis and listening
remain available for unknown permissions. The metadata marks those packs as
local-only. `reference_provenance.py --check-public PACK` rejects a pack unless
every reference has an explicit allowed status, permission source and intact
notices. This is an evidence check, not a legal opinion or a publishing command.

Release staging uses an explicit inventory and rejects audio, downloaded
archives, build trees and private audit records. The only generated dependency
archive permitted in corresponding source is the original hash-pinned PFFFT
archive. The extracted TETRA codebook data is a separately disclosed retained
component; see [ETSI-RETAIN-1](status.md#etsi-retain-1).
