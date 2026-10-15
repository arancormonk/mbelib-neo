<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# TETRA ACELP speech decoder

`src/tetra/` decodes TETRA full-rate speech: the ACELP codec of ETSI EN 300 395-2 clause 4,
137 type-1 bits per 30 ms frame in, 240 samples at 8 kHz out, with the bad-frame concealment of
clause 4.2.3.2. The encoder, channel coding (clauses 5-6), the AMR option (clauses 9-12) and
homing-sequence detection are out of scope. The API is described in the README
([TETRA decoding](../README.md#tetra-decoding)) and `include/mbelib-neo/mbelib.h`.

## Status

EN 300 395-2 clause 8 makes ETSI's fixed-point C code the normative definition of the codec, and
EN 300 395-4 tests a decoder bit for bit: TE1_46.CHD decoded must equal TE2_128.PCM. **This
decoder passes that test: its output equals TE2_128.PCM byte for byte.** It is written from the
text of the standard and from black-box measurements of ETSI's reference decoder. Outside the
test sequence two differences from the reference decoder are known: bad frames received before
the first good frame, and one synthetic stream of random bits whose output clips at full scale
(see [Known differences](#known-differences)).

The decoder uses integer arithmetic only, through saturating 16/32-bit operators with portable C
semantics (`src/internal/mbe_fixed_point.h`), so its own output is identical on every compiler,
architecture and build configuration (`test_tetra_decoder` asserts stream hashes everywhere).

## Clean-room process

- Sources: the text of EN 300 395-2 (V1.3.1, 2005-01, and V1.3.3, 2025-02), general fixed-point
  speech-coding practice (the semantics of the usual basic operators), and black-box observation
  of ETSI's reference decoder: frames in, PCM out.
- Never opened: ETSI's .C/.H files, the C-WORD documents, and derivatives of ETSI's code
  (osmo-tetra's codec patches, tetra-kit's codec, install-tetra-codec, dsd-neo's
  `tetra_codec_etsi.c` beyond its ABI header).
- The only ETSI files read are C-CODE/CLSP_334.TAB and C-CODE/ENER_QUA.TAB, by
  `tools/tetra/extract_acelp_tables.py`, which copies their numbers and nothing else.
- The reference decoder was built from ETSI's archive with dsd-neo's `tools/tetra/acelp_adapter`
  (`-DBUILD_REFERENCE_SDECODER=ON`, 32-bit, no source edits), and reproduces TE2_128.PCM from
  TE1_46.CHD byte for byte. It was only run on input frames.
- ETSI archives, the conformance files, the reference build and its output stay outside the
  repository. Test hashes in this repository are hashes of this decoder's output on
  self-generated input.

## Decoder

Clause references are to EN 300 395-2. Values in Qn have n fractional bits.

- **Bits** (4.2.2.7): table 3 order; table 3 misprints subframe 3's pulse fields
  ("B95-B99", "B100-B102"), and table 4 gives pulse 2 = B95-B97 and pulse 1 = B98-B102.
- **LSPs** (4.2.2.3, 4.2.3.1.1): the three split-VQ sub-vectors (Q15 cosines) are looked up, and
  the pairs on either side of a split boundary closer than 916 (q3/q4) or 1244 (q6/q7) are pushed
  apart, each moving by half the shortfall rounded down (so an odd shortfall leaves the gap one
  short of 916 or 1244). If any adjacent pair is then still not strictly decreasing
  (q(i) <= q(i+1)), the frame's LSPs are discarded and the previous frame's are used, for
  interpolation and as the next frame's previous vector.
- **Interpolation** (equation 22): weights 3/4, 1/2, 1/4 of the previous vector in Q15, the
  result rounded down.
- **LSPs to A(z)** (4.2.2.2): the sum and difference polynomials in Q24 with a 32 x 16-bit
  double-precision multiply that keeps the product of the low half's top bit exactly: (hi, lo) q
  = 2 hi q + b q + 2 floor(l q / 2^15), with b bit 14 of lo and l its lower 14 bits (for an even q
  this equals the usual 2 hi q + 2 floor(lo q / 2^15)). A(z) coefficients are rounded to Q12 (16
  bits; a coefficient beyond +-8 wraps).
- **Pitch** (4.2.2.4): subframe 1 has thirds from 19 1/3 to 84 2/3 (indices 0-196) and integers
  85-143; subframes 2-4 have 32 thirds from tmin - 2/3 around subframe 1's integer lag T1, with
  tmin = T1 - 5 limited so that [tmin, tmin + 9] lies within 20..143.
- **Adaptive codebook** (4.2.3.1.2): fractional lags use two 32-tap filters, the Hamming-windowed
  sinc of 4.2.2.4 (window 0.54 + 0.46 cos(pi t / 16) over +-16 samples, i.e. +-48 thirds) sampled
  at +-1/3 and rounded to Q15; a lag below 60 repeats the subframe's own new samples.
- **Algebraic codebook** (4.2.2.5): pulse positions 2i, 2 + 8j, 4 + 8j and 6 + 8j, plus the
  shift bit; a position of 60 or more is outside the subframe and contributes nothing (pulse 1
  indices 30 and 31, and pulses 3 and 4 at index 7). Amplitudes +1.4142 (Q14 23168), -1, +1, -1,
  all negated by the sign bit. The pulses are filtered by F(z) = A(z/0.75)/A(z/0.85) (gamma powers
  as recursive rounded Q15 products) and, for an integer lag T below 60, sharpened by
  f(n) += 0.8 f(n - T); each sample is rounded down.
- **Gains** (4.2.2.6, equations 31-39): energies in Q8 log2, both 0 after a reset. Prediction
  0.5 E(p) + 0.25 E(c) - 3 (and symmetrically for the code), computed with the quarter rounded
  down and floored at zero, plus the codebook's prediction error; the results are limited to 27
  and 25. Each gain is 2^((E - Ea) / 2) (equations 38 and 39), Ea being the log2 of an energy
  product (equations 31 and 33) formed from 16-bit mantissas: pg, the energy of the synthesis
  filter's response to an impulse of 1024, and the adaptive-vector energy 1 + sum v^2 are
  normalized and their high halves taken; the code-vector energy sum c^2 contributes its high half
  at a fixed scale, without normalization. Each product of two mantissas goes through one log2
  (33-entry table round(32767 log2(1 + i/32)) with linear interpolation), the normalization
  shifts are subtracted (both for the pitch gain, pg's for the code gain), the result is truncated
  to Q8, and one offset (4434 in Q8) accounts for the Q formats of both gains. E - Ea, read in Q9,
  is converted by a 33-entry pow2 table; the gains are Q12 and the pitch gain is limited to 1.2
  afterwards (the stored energies are not changed by the limit).
- **Excitation and synthesis** (equations 40-41): u = gp v + gc c rounded to 16 bits (it is also
  the adaptive codebook's memory); 1/A(z) with Q12 coefficients and a saturating 32-bit
  accumulator; the output doubled with saturation (4.2.1).
- **Concealment** (4.2.3.2): previous LSPs; the previous good frame's subframe-4 integer lag in
  all four subframes, fraction 0; both energies lowered by 0.5 per subframe and not below zero;
  the previous good frame's innovation indices. A bad frame's bits are not read, and any nonzero
  indicator marks a frame bad.

### Spec-silent details resolved against the reference

The boundary gaps 916 and 1244, the rejection of LSP vectors left crossed (gap 0 rejected, gap 1
kept), the pitch index mapping, the interpolation filter taps, the pulse-1 indices 30 and 31,
the order of the sign and the rounding in the codevector, the prediction's floor at zero, the
gain computation, the multiply of the LSP polynomials and the concealment energy floor were each
settled by comparing candidates with the reference decoder's output on TE1_46.CHD and on frames
generated to isolate one stage. For the gains, the candidates followed the operator counts in the
standard's informative complexity tables (annex C table C.9 and annex D tables D.25 to D.31: two
log2 and two pow2 calls per subframe, two normalizations and two 16 x 16-bit products in the
energy decoding), which rule out separate logs of the three energies. For the polynomial
multiply, frames with a single code pulse at the end of every subframe and a lag of 120 make the
reference's output decay freely through 1/A(z), so its Q12 coefficients can be read back from
the PCM. They equal those of the usual double-precision product whenever every LSP is even, and
LSP transitions whose interpolated vector has a single odd value isolate the one to four
multiplications that can differ; their measured deviations (over 1 300 single-multiplication
cases) fit the multiply above without exception.

### Departure from the reference

After a reset, bad frames that arrive before any good frame decode to silence in both decoders.
The reference then starts the first good frame late and sometimes with a burst (in one test the
first good frame peaked at 10 730 where a fresh decode peaks at 216). This decoder instead keeps
its reset state until the first good frame, which then decodes as it would from a reset. The text
does not cover this case.

## Known differences

- **Bad frames before the first good frame.** This decoder keeps its reset state (see the
  departure above); the reference does not, so the first good frames after such a start differ.
- **Clipping in random data.** In one synthetic stream of random bits with 5% of the frames lost
  (one of `test_tetra_decoder`'s streams), the output differs from frame 1283 for four frames, in
  a passage that clips at full scale. Every other stream tested decodes identically.

## Measured results

`tools/quality/tetra_conformance.py compare` (see [testing](testing.md#tetra-conformance))
compares the output with the reference decoder's and exits 0 only for identical files.

- TE1_46.CHD decodes to TE2_128.PCM exactly (5 340 of 5 340 frames).
- TE3_46.SPE (EN 300 395-4's encoding of TE2_128.PCM, a second stream of the same speech) decodes
  exactly as the reference decoder decodes it.
- TE1_46.CHD with bad frames set by `tetra_conformance.py bfi` decodes exactly as the reference
  decodes it for `--rate 1`, `--rate 3` and `--rate 10` (seeds 1 and 2 each) and for
  `--burst 100:20 --burst 1000:50 --burst 3000:300`. With `--first 3` the two differ from frame 3
  (the bad-start departure): 5 197 of 5 340 frames exact, segmental SNR 97.2 dB, frame 10
  (-31.1 dBFS) at 13.7 dB.
- Concealment on its own was compared on 34 twelve-frame excerpts of TE1, each with nine bad-frame
  patterns (a single bad frame at three places, runs of 2, 3 and 6, two bad frames around a good
  one, and one or two at the start; 306 cases): the first concealed frame is exact in all 306
  cases and the whole excerpt in all 238 cases without a bad start.
- Of the seven synthetic streams behind `test_tetra_decoder`'s regression hashes, five decode
  exactly as the reference decodes them; the bad-start stream and the random stream with losses
  differ as described above.

For runs that are expected to differ, `compare --gate` takes its exit status from tolerance
thresholds instead: a segmental SNR of at least 30 dB over frames above -60 dBFS, a level
difference within 0.1 dB, no frame above -50 dBFS below 14 dB SNR after the first 10 frames, and
reference silence reproduced within 2 LSB.

## Tables

The decoder needs two tables that exist only as data in ETSI's archive: the LSP split-VQ
codebooks (CLSP_334.TAB) and the gain codebook (ENER_QUA.TAB). `tools/tetra/extract_acelp_tables.py
ARCHIVE` checks the archive against the SHA-256 of en_30039502v010301p0.zip (V1.3.1) and
en_30039502v010303p0.zip (V1.3.3), whose C-CODE members are identical, reads those two members
only, and writes `src/tetra/tetra_acelp_codebooks.c` with the archive and member hashes;
`--check` compares instead of writing. The file carries `LicenseRef-ETSI-EN-300-395-2-data`: the
values are (C) ETSI and they are not covered by the project's GPL grant
(`LICENSES/LicenseRef-ETSI-EN-300-395-2-data.txt`). They are carried for an implementation of the
standard that passes its decoder test bit for bit: ETSI's IPR Policy (clause 9.2.2, with SOFTWARE
defined in clause 15, definition 14, to include codec code and data) grants implementers a
royalty-free licence to software contributed to a standard "for the purpose of using the SOFTWARE
in any compliant implementation". Whether that grant reaches these 1995 contributions has not been
established, a contributor could have chosen fair, reasonable and non-discriminatory terms instead
(clause 9.2.2(ii)), and the grant is limited to that purpose, so it does not make the data
available under the GPL; no licence that does is known. Other open-source decoders of standardized
speech codecs, for example FFmpeg's AMR and G.729 decoders, carry their standards' tables in the
same way.

Everything else is computed from formulas by `tools/tetra/gen_acelp_constants.py` into
`src/internal/tetra_acelp_constants.h` (log2/pow2 tables, interpolation filters, gamma powers);
`test_tetra_tables` regenerates it and compares.

## IPR search

ETSI's IPR database exports (ISLD-export.zip and GD-export.zip) were searched on 2026-10-15 for
declarations whose standard names "300 395" in any ETS/EN form, and for every declaration on the
TETRA project without a standard, with titles screened for speech and codec terms. Two bear on
the codec:

- ISLD-200101-005, Motorola, against ETS 300 395-2: US 5,097,508 / EP 0 415 163 (priority
  1989-08-31) and US 5,359,696 / EP 0 450 064 / WO 91/03790 (priority 1988-06-28).
- ISLD-190001-073, Thomson-CSF, TETRA project: the algebraic-codebook family US 5,444,816,
  US 5,699,482, US 5,701,392, US 5,754,976, EP 0 516 621 and CA 2 010 830 (priority 1990-02-23),
  and EP 0 275 751 (priority 1986).

Google Patents lists their anticipated expirations (it disclaims legal conclusions):
US 5,097,508 2009-08-31; EP 0 450 064 2010-06-25; US 5,359,696 2011-10-25; US 5,444,816
2012-08-22; US 5,701,392 and US 5,754,976 2012-09-10; US 5,699,482 2014-12-16; EP 0 275 751
2007-12-15; the other European filings date from 1990 and are past 20 years. GD-export holds
general declarations (29 rows marked applicable to all standards) that name no patents; none is
TETRA-specific. The exports also hold national family members that were not checked one by one,
and ETSI does not guarantee that the database is complete. Other TETRA project declarations
(Alcatel, Motorola, Nokia, Ericsson, OTE, Securicor) have air-interface or system titles.
