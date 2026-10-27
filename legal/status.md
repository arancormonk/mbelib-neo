<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# Licensing status

Planned release **2.6.4**, retaining `SOVERSION 2`. The original licensing review
was recorded on 2026-10-17 for release **2.6.0**.
The project licenses its contributions under GPL-2.0-or-later. This objective
does not establish permission for every component in the combined library.
[components.json](components.json) is the maintained inventory used by release
SBOMs and checks. Declared licenses, copyright permission, GPL compatibility,
and patent status are separate fields. Aggregate conclusions remain
`NOASSERTION`; checks enforce records and packaging, not legal clearance.
Notice hashes use UTF-8 text with LF line endings, including on Windows;
archive and SBOM checksums always cover the actual delivered bytes.

Original mbelib contributions retain their ISC notice. PFFFT's exact retained
license texts use `LicenseRef-PFFFT`, covering its compiled backends. Both sets
of notices accompany installations. The OP25 patch is GPL-3.0-or-later and
retains Max H. Parke's attribution and the complete GPLv3 text. It is separate
developer tooling in the source archive, never part of the installed library.

## ETSI-RETAIN-1

TETRA remains enabled. Copyright permission and GPL compatibility for
`src/tetra/tetra_acelp_codebooks.c` remain **unresolved**. Its custom identifier
and attribution are a notice, not an established redistribution license.

The [ETSI IPR Policy](https://www.etsi.org/images/files/IPR/etsi-ipr-policy.pdf)
retrieved 2026-10-17 is dated December 2022. Clauses 9.2.1–9.2.3 describe
limited software grants and a contributor alternative; 9.2.4 excludes implied
patent grants. Outstanding questions are whether the policy applies to the
1995 table contributions, who held the necessary rights, whether alternative
terms were elected, and whether any applicable restrictions permit GPLv2
redistribution and modification. The policy's software definition alone does
not answer these questions. Conformance results and tables in other projects
do not establish authorization.

The explicit release exception is to retain the tables with this disclosure
and `NOASSERTION` for their assessed license. It is not a finding of clearance.
Closing it requires evidence of applicable rights and terms, with compatibility
review, or independently authorized replacement data. No rights-holder contact
or legal opinion has been obtained as part of this engineering change.

## Patents and release review

AMBE+2 remains enabled. [The patent review](patents.md) records implementation
observations and open questions, not a non-infringement opinion. Educational,
research or local use is not asserted as an exemption. The GPL grant does not
confer third-party patent rights.

Update the inventory, provenance record and patent review when affected codec
implementations change and during every release preparation. Review evidence
and actual status; passing a predicted expiration date does not close an issue.
Record dates and sources, and preserve unresolved conclusions until supporting
evidence changes them. External reference data follows
[reference-data.md](reference-data.md).

Releases before 2.6.0 are not retroactively cleared by these changes. A local
historical exposure and source-backfill record identifies affected tags and
their original dependency revisions. Existing tags, release assets and Git
history remain unchanged. Backfill must use each release's own source and
dependency, and retain the disclosures applicable to that release.
