<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# Testing Policy

mbelib-neo uses automated tests, static analysis, sanitizers, and fuzzing to
reduce regression and security risk.

## Test Suites

The CTest suite includes:

- API/version/result helper checks
- ECC tests for hard and soft Golay/Hamming paths
- AMBE 2400 encoder round trips, exact `log2Ml` and frame-mark parity (on frames with and without voiced bands), pitch endpoints, independent contexts and reset replay in one thread, plus behaviour on synthetic speech: noise of any colour and level stays unvoiced, steady vowels from 70 to 310 Hz are voiced and on pitch, harmonics below 2 kHz with noise above voice only the lower bits, and 240 Hz, missing-fundamental, strong-second-harmonic, glide and onset cases show no octave errors; the decoded level follows the input down 40 dB and silent or -70 dBFS input is coded as voice frames that decode near silence; the 16 DTMF digits, the call-progress tones and single tones at several phases are sent as tone frames in DVSI's layout with the index the decoder plays at their frequencies, decode within 0.5 dB of the input level and leave the prediction history unchanged, KNOX tones stay voice, voice-tone-voice prediction stays with the decoder's, a malformed history is rejected on a tone frame without touching state, and call-progress tones are timed as below (`test_ambe2400_encoder`)
- AMBE+2 2450 encoder: frame FEC round trips of all 49 bits and single-error correction; decoded frames requantized to the same model (gain, prediction, block DCTs, PRBA and higher-order quantizers, packing); the TIA-102.BABA-1 eq 4 voicing quantizer against an independent evaluation; DTMF, KNOX, call-progress and single tones at several phases sent as tone frames the decoder plays at the input level, the DTMF frequency and twist limits, voice-tone-voice prediction, and a malformed history rejected on a tone frame without touching state; call-progress timing, checked for the D-STAR encoder too: a tone is sent once it has filled three spans and for two frames after, 45 ms bursts at any alignment stay voice, a frame rejected during that confirmation leaves the stream unchanged, a 20 ms dropout is bridged and a DTMF digit that follows is sent at once; tone edges, also checked for the D-STAR encoder: 60 ms DTMF and single tones are sent for three or four frames wherever they start, and a DTMF digit that changes directly to another leaves no voice frame between them wherever the change falls, also at a level just above the detector's floor; plus the checks shared with IMBE below (`test_ambe2450_encoder`)
- IMBE 7200x4400 encoder: the (15,11) Hamming encoder over every data word, frame FEC round trips of all 88 bits and single-error correction, the alternating synchronization bit, and decoded frames of every harmonic count requantized to identical bits (`test_imbe7200_encoder`)
- Shared by both (`tests/encoder_test_support.h`): `cur_mp` bit-identical to the decoder chain, its frame mark included, with `prev_mp` untouched, rejected samples leaving the stream unchanged, prediction boundary taps, 16-bit/float parity, independent contexts and reset replay, noise staying unvoiced, steady vowels from 70 to 310 Hz voiced and on pitch, no octave errors (240 Hz, missing fundamental, strong 2nd harmonic, glide, and 133, 200 and 390 Hz harmonic series), the decoded level following the input down 40 dB, quiet input coded as voice, and a one-frame flush. Also for the D-STAR encoder: a stale frame mark handed in with `cur_mp` replaced by the decoder's, and frames with no voiced band sent at the codec's fixed pitch code and harmonic count: noise decoding at the input level, vowels alternating with noise, a frame rejected at the fixed code (shown to take the fixed code with a valid history) leaving the stream and `cur_mp` untouched, and noise band-passed to 1-2 kHz and 2.4-3.2 kHz decoding with the harmonics inside the band at least 20 dB above those outside it on the fixed code's grid
- AMBE joint vector search, for D-STAR's and AMBE+2's tables: its linear model reproduces the decoders' own mean-removed residual for random rows at every L from 9 to 56, the refinement never raises the amplitude-weighted error measured from the decoders' reconstruction and sends only rows the decoder receives (D-STAR's b8 even, 0 for blocks of at most two harmonics), and a frame the decoder reproduces exactly keeps the nearest-row search's rows (`test_ambe_vector_search`)
- AMBE block lengths (`test_ambe_block_lengths`): the four tables, the D-STAR and rate-33 decoders' and the two encoders' (`mbe_ambe2400_enc_tables()`, `mbe_ambe2450_enc_tables()`), each row of L 9 to 56 summing to L and equal to the test's copy of TIA-102.BABA-1 Annex C but row 20, 4, 5, 5, 6 where the annex prints 4, 4, 6, 6 (rows 0 to 8 zero). The decoders, in each codec, on pairs of frames built from chosen codes at 19, 20 and 21 harmonics (every pair of the three, the second predicted from the first; every code of 20 harmonics, D-STAR b0 54 to 57 and rate 33 b0 51 to 53; four D-STAR voicing fields and twelve rate-33 V/UV rows, the band-grid rows 18, 21, 24, 26, 29 and 31 among them, which the encoder never sends, so only the decoder is tested there), through the public parameter decoder: (1) each frame's `log2Ml` within 2e-4 of a transform formed in double (PRBA rows, 8-point inverse DCT, each block's coefficients and inverse DCT at the expected lengths, the prediction at the position the decoder rounds), and every division that moves one harmonic between two neighbouring blocks (and 5, 5, 5, 5 at 20) more than 0.05 from it (least margin measured 0.405 D-STAR, 0.732 rate 33). The 2e-4 is about 210 single-precision ulps at a value of 10 (more below): the decoder's sums each round in float and the double reference does not; it catches a wrong row by three orders of magnitude and no smaller error. (2) Each frame's `log2Ml` and gamma against an exact result formed in float by the decoder's arithmetic written out (`exact_frame()`: the codebooks from the codecs' headers and D-STAR's three constants, `cosf()` and `log2f()`, every operation in the decoder's order), compared by bit pattern where the library is IEEE, `FLT_EVAL_METHOD` is 0 and contraction is off (dev-debug, asan-ubsan-debug, notones-debug, the SIMD Debug build and the `-march=x86-64-v3` Debug build: measured difference 0) and within 1.6e-5 elsewhere (dev-release 5.25e-6, the two Clang fast-math LTO builds 5.72e-6, 32-bit Release 1.43e-6, 32-bit Debug 0, the largest differences measured). An error of an ulp or a few (the review's 1e-6 at L 20) fails it in the five exact configurations and passes in the five others. The golden PCM sequences hold no frame of 20 harmonics, and rate 33's encoder takes its model from the same public decoder, so neither they nor the parity check pin the arithmetic at L 20: the exact comparison does, where it is exact. The encoders: (3) decoder and encoder state parity by bit pattern (w0, L, gamma, the frame mark, and Vl, log2Ml and Ml over 1..L, the fields the encoders document as the decoder's) over a 1,000-frame signal at about 180 Hz that sends frames of 20 harmonics (147 D-STAR, 120 rate 33; D-STAR's existing parity test reaches 35 such frames too), an encoder error ending it before its bits are read and a tone status told from an error. Parity shows that each encoder predicts from the history the decoder will hold, not that either uses the right row: a codec's table header serves its decoder and encoder alike (a wrong row in it leaves no parity mismatch), and rate 33's encoder takes its model from the decoder; only a table copied into one translation unit of D-STAR's breaks it. (4) The joint refinement acting at L 20: on 300 random frames of 20 harmonics the shared quantizer `ambe_enc_quantize_residual()` changes the nearest rows on most (225 D-STAR, 206 rate 33 measured) and lowers the weighted error in all (by 4.6% and 4.9%), and rate 33's codec path (`mbe_ambe2450_quantize_amplitudes()` at b0 52) sends the shared quantizer's rows on every frame; D-STAR's own encode path is static and not reached. (5) Residuals that blocks of 4, 5, 5, 6 represent exactly (random rows decoded by the search's model of the decoder) coded through the quantizer within 1e-4 rms with the encoders' own tables and at least 0.05 rms off in the mean with row 20 put back to the annex's (measured 1.6e-7 and 1e-7 against 0.209 and 0.206). An error that is not finite fails every comparison before it is aggregated. Mutations, each failing it (Debug, the whole suite run for each; other tests failing too named): row 20 back to 4, 4, 6, 6 in D-STAR's table, in rate 33's, in both; 4, 5, 6, 5 and 5, 5, 5, 5 in both; row 19 4, 5, 4, 6 in both and in D-STAR's; row 21 4, 5, 5, 7 in both and in rate 33's (`test_golden_pcm`); a local copy of the table with the annex's row 20 in D-STAR's encoder or decoder (`test_ambe2400_encoder`, `test_ambe2400_frame_types`) or rate 33's (`test_ambe2450_encoder`, `test_ambe_vector_search`); the annex's row restored for rate-33 band-grid frames of 20 harmonics only; `log2Ml[1]` off by 1e-6 at L 20 in rate 33's decoder and in D-STAR's (the latter `test_ambe2400_encoder`, `test_ambe2400_frame_types`); the joint refinement skipped at L 20; a NaN out of the search's model of the decoder at L 20 (`test_ambe_vector_search`), and only for the own row (the same); an encoder failing on its 400th call (`test_ambe2400_encoder`, `test_ambe2450_encoder`); and, in the test itself, a NaN for the own-table round trip's error and a NaN in an alternative row's reference. At `30249f6` the restored grid row, rate 33's 1e-6 error, the skipped refinement and the two NaNs in the test passed it (checked).
- MBE speech analysis numerics: Kaiser windows, the window transform, harmonic fits and magnitudes of off-grid harmonics, noise statistics without magnitude ripple, the FFT autocorrelation against the direct sum, the voicing thresholds, and the measures taken again on another grid (`mbe_analysis_measure_at()`: five tones at unequal levels on harmonics of the IMBE fixed grid, which the analysis does not track, up to the top of the harmonics IMBE codes; every harmonic below Nyquist returned, each tone read in its own harmonic within 1 dB, its neighbours at least 10 dB and every other harmonic at least 25 dB below) (`test_speech_analysis`)
- Encoder context, FFT and autocorrelation plan allocation failures (each of the eleven, for the D-STAR, AMBE+2 and IMBE encoders) clean up fully; after successful allocation, encoding (tone frames included) and reset allocate nothing (`test_ambe2400_encoder_oom`, `test_ambe2450_encoder_oom`, `test_imbe7200_encoder_oom`, GNU link wrapping when LTO is disabled)
- noise determinism and frame-state determinism checks
- parameter and synthesis behavior checks, including an `L * w0 < pi` bound for every model the decoders emit and rate-33 V/UV rows 18-31 on their band grid (`w0 = 2.81 / L`, the unvoiced factor from it, the next frame's prediction unchanged), V/UV row 17 decoded by the public parameter decoder with no voiced band and every amplitude at the unvoiced scale, and the marks of rows 13, 15 and 17 on every row's decode, and P25's decoded w0 for all 208 pitch codes against an independent evaluation of DVSI's pitch (eq 46's nominal pitch in float times `1 + g` in double, rounded once: bit for bit where the library is built without fast-math, the 32-bit x87 build included, and within an ulp with it; the law wholly in double, unrounded, within 1.5 ulps), with L and K from the nominal pitch and ProVoice keeping the nominal pitch (`test_params`)
- the codec paths' unvoiced synthesis: DVSI's noise generator (its period, its first samples, a cold start and a skipped frame, and the samples and next state of a frame drawn from each of its 65,536 states against the recurrence), a subframe against the equations it implements, the window's squares summing to one, a level step inside a frame, the public `mbe_synthesizeSpeechf()` keeping the specification's rendering, a state handed from one rendering to the other, invalid arguments and models, and two streams replaying alike (`test_unvoiced_subframes`)
- P25 tones sent as voice frames: frames built with the IMBE quantizer at a tone's pitch code and voicing word pass the decoder's tone test (a DTMF tone, a single tone, the dominant harmonic's row where single tones share a pair) and fail it with another voicing word or a call-progress pair; the test's statistic is the rms of all the other decoded magnitudes (the square root of the mean of their squares) over that of the tone's, and with one other harmonic (the 10th of DTMF 1's frame) raised step by step a frame passes exactly while that statistic stays under 0.375, and a single tone, its other harmonics raised together, exactly while it stays under 0.0235; a tone frame's synthesis model (the 2^-19 pitch, the tone harmonics at one level over the synthesis droop, the others unvoiced at 0, the mark) with the decoded frame kept as the prediction history, its own model at mid-frame, advancing it as rendering it, a repeat playing the history's tone model when the frame the repeat criteria reject is another valid tone frame, and ProVoice untouched; three tone frames rendering two steady components at one model amplitude and the frame after fading the tone out; words that the decode's own test, run again on the model they sit in, does not give playing as no word, sample for sample as the same history without one (a code past the table or without rows, a row of another code, another code's row in this code's model, the real word in a model of another pitch code, a tone magnitude negative, zero, subnormal, infinite or not a number, code 1's word in a real code-2 tone frame's history, of the same L and K, and in a real code-1 frame of the same voicing whose other magnitudes fail the single-tone limit); every one of the 116 rows through a real decode playing the tone model, and three repeats after it; a tone mark with no tone model giving way, through the P25 and ProVoice paths, to the mark a voiced model carries (`MBE_FRAME_MARK_ONSET_COLUMNS`) and to no mark on a model with no voiced band, in all three states, on a frame without channel errors (the decoded mark the caller's state overwrote not recovered; where smoothing voices a harmonic, the model rendered carries the onset mark instead), the decoders' own marks kept; the public decoder clearing the word on a special frame without touching the history's word, and a status error (an invalid bit after a tone frame, a missing state) writing nothing; ProVoice keeping 0 in every `Vl[0]`, and a P25 repeat after a ProVoice frame playing the repeated model as speech; a tone across the repeat limit (three repeats play it, the fourth mutes, the reset clears the word and the mark); and the table pinned (its 116 rows and the index, by structure and an FNV-1a hash that every one-unit change of a field or an entry was checked to move) (`test_p25_tones`)
- AMBE+2 and D-STAR tones through the public frame paths: each at the level of the DVSI-calibrated laws and at the frequencies DVSI's decoders play (harmonics of a fundamental on the 2^-19 grid: for the dual tones Annex J's, except IDs 145 and 147, whose fundamentals are DVSI's own, with DTMF "1" at 693.909 + 1214.340 Hz; for the single tones their own rule, harmonic ID / 16 + 1, not Annex J's harmonic numbers, with ID 46 at 1437.469 Hz, harmonic 3 of 31402 units; D-STAR's tones as the AMBE+2 tone of the same nominal pair) (`test_tones`)
- the codec paths' voiced synthesis in subframes: a steady model adding up to its harmonics, an onset starting the pulse train, a late onset (in the rate-33 rendering no band voiced at mid-frame whatever the frame before, and after a voiced frame the train running on through the first subframe at that frame's pitch and through the second at this one's; in the other renderings acting after a frame with no voiced band only), a pitch step keeping the train, the train taken an octave down or up before a subframe step whose period ratio lies within 20% of 2 or of 1/2 (`test_octave_steps`: ratios 1.70 and 0.55 take the step and 1.55, 2.45, 0.62 and 0.39 glide, each case first shown to tell the two apart by 0.05 cycle, the stored pulse phase within 1e-5 and in [0, 1) as stored, read from its bits, and advancing leaving the same pulse state, `PSIl` and `PHIl` as rendering; `test_stored_phase_range`: a phase that reaches 0.9999999917 cycle at a frame's end stored as 0), the mid-frame model against an independent evaluation, phases following `log2Ml` rather than the rendered amplitudes (the mid-frame model's between frames of different pitch too) and untrusted `log2Ml` values, the mid-frame pitch of the subframe renderings (`test_ambe2450_mid_pitch`: rate 33's, between voiced frames of different pitch only, the geometric mean times the stored factor 1 - 3.0e-4 (the float, promoted) times the half-octave wave, against the formula computed in the test in double and rounded to float once, bit for bit where IEEE arithmetic gives the library and the test the same bits (no fast-math, no excess precision, no fused multiply-add the compiler may contract) and within an ulp elsewhere, for pairs of rate-33 pitch codes as the decoder stores them, among them the pair beside the wave's maximum over the codes and the pair beside its minimum, and off the code grid; D-STAR's and the IMBE rendering's the geometric mean in float for the same pairs, bit for bit in an IEEE build; equal pitches that pitch in every rendering, bit for bit in an IEEE build and within an ulp with fast-math, for pitches a decoder produces (equal pitches take the inherited `sqrtf(a * b)`, which overflows or underflows at the ends of the float range, as before the wave); one or no voiced frame, a pitch-carrying frame and a silence frame; pitches at the ends of the float range: the two largest floats, in either order, whose pitch passes `FLT_MAX` before rounding, giving no pitch and a mid-frame model with no voiced band (a fallback that returned `FLT_MAX` fails), the others a positive normal pitch, and finite samples; `test_ambe2450_mid_pitch_codes`: every ordered pair of different rate-33 pitch codes, the pitch as above and the mid-frame model's harmonic count the one the formula's pitch gives, including the pairs whose count the wave changes, 0 -> 94 among them, so that a count taken from the pitch without the wave fails), advancing and rendering leaving the same pulse train in the rate-33 rendering, the waveform test's independent rendering (which regenerates its phases from `log2Ml` on its own, computes rate 33's mid-frame pitch with the wave itself, from the stored factor, and the mid-frame model's harmonic count from its own mid-frame pitch, matching the library's bit for bit where IEEE arithmetic allows, over shaped spectra rendered with a shelf, has its own copy of the octave rule, through steps an octave down and up and a glide, its train running on from a frame marked `MBE_FRAME_MARK_MID_COLUMN_0` in the rate-33 rendering alone (a row-17 frame followed by a row-13 frame), and computes the mid-frame voicing itself from the rules as the parity document states them, `reference_mid_voicing()`, in D-STAR's and P25's renderings the 500 Hz columns' rule, from its own column samples in double, every pitch there played 0.1% above the one listed so that no position the rule rounds falls on a whole number) in all three codec renderings (D-STAR's, P25's and rate 33's), linear crossfades into and out of silence, a muted frame moving the train and the noise generator on in all three renderings, advancing leaving the train, the phases and the noise generator where rendering leaves them through every mark that shapes the mid-frame model (`test_advance_follows_the_marks`, in all three renderings: P25 onsets read in columns before a voiced run and the mark changing nothing inside one, late onsets after unvoiced and voiced frames, rate 33's rows 13 and 17), without a transform plan and with a real one (where a rendered frame leaves its own subframe's second half in `previousUw` and an advanced one none), stored pulse state that is not a frame's own being dropped, invalid models (harmonic counts, pitches and magnitudes of any value) leaving the neighbouring state alone (a harmonic count outside 1..56 in either frame, in all three renderings, leaving the current state byte for byte as it was, rendered and advanced: P25's rendering marks the model it renders only on a frame it takes), untrusted unvoiced overlaps (`test_untrusted_unvoiced_overlap`: not a number of either sign, plus and minus infinity, plus and minus `FLT_MAX` and 2^127 in the subframe rendering's `previousUw`, directly and through the P25 data path, and in `mbe_synthesizeSpeechf()`'s `previousUw` and noise overlap: the frame finite, and playing as with 0 there where a `previousUw` sample is not a finite number or a noise-overlap sample lies outside the 0 to 53124 its generator writes, 53124 itself playing as stored, and a huge finite `previousUw` sample playing as stored: its own value from the subframe rendering alone, the output limit of its sign through the P25 data path and `mbe_synthesizeSpeechf()`; a guard against a defect inherited from main, found in review), the noise-overlap helper itself (`test_noise_overlap_guard`: in-range values, 0.5 and a subnormal among them, copied bit for bit, every value outside the range read as +0, and a buffer that overlaps the overlap filled as `memmove()` fills it), and the frame marks: each stored as its value in float and read only from that exact value, 6 and 7 included (`test_frame_marks`), every P25 and D-STAR pitch code of an unvoiced frame, the previous voicing lasting to mid-frame through the P25 data path (`test_p25_voicing_persists`), held there in the bands whose column, or the column 250 Hz higher, the frame before covers, D-STAR's 500 k Hz and P25's from 505 Hz (`test_voicing_held_by_columns_without_a_pitch`), the 500 Hz columns a P25 onset starts with (`test_p25_onset_columns`: a frame voiced in its first band only late at pitch code 56 and early at 55, harmonics 1 and 2 sounding at mid-frame and the third not; the word `1000001100` early at 86, harmonics 19 to 23 sounding, and late at 92; `1011` then unvoiced bands early at 79, harmonics 7 to 11, and late at 92, 99 and 102; a frame voiced throughout at code 0 silent at mid-frame at its 9th harmonic, 3.65 kHz; through the P25 data path the pulse train starting with the frame or at mid-frame, ProVoice alike), the onset-columns mark changing nothing after a voiced frame and a frame played as a P25 tone model keeping its own model (`test_p25_onset_columns_only_at_onsets`), the mid-frame voicing of D-STAR and P25 by 500 Hz columns against the rule stated in the test on the decoders' own inputs (`test_dstar_mid_columns_by_words`: all 256 pairs of D-STAR voicing words at twelve pairs of pitch codes, the same code, b0 0 (418 Hz, two or three harmonics a band) among them, and steps between codes, from the words' band bits, onsets of all fifteen voiced words and each followed by a frame lacking a pitch included, and in the rate-33 rendering this frame's voicing, or the hold through a frame lacking a pitch that rendering keeps; `test_dstar_mid_columns_table`: 1101 -> 1111 and 1111 -> 1101 starting late and letting go column 5 under column 6 voiced in both, 1000 -> 0100 voicing column 1, an onset starting the top column late, a frame lacking a pitch letting the top column go, a frame at an unvoiced code holding nothing; `test_p25_mid_columns_by_words`, from the voicing words: covered and uncovered columns and the 250 Hz reach at onsets, the isolated last band of a frame of 9 harmonics (word `001`) covering column 5, the exception at 37 harmonics and not at 36 nor with band 11 voiced, frames lacking a pitch after one voiced from 1005 Hz (held) and after ones voiced in isolated bands only (nothing held), in-run changes both ways and pitch steps, and the rate-33 rendering as for D-STAR; `test_p25_mid_columns_cases`: the exception's harmonics, a frame after a tone model taking the rule as after the same model unmarked; `test_p25_nearest_harmonic_across_a_step`: b0 134, word `000000000010`, then b0 161, word `111111111101`, the mid-frame harmonic 36 at 33.4999986 of the first frame's harmonics read at 33 (voiced) and kept voiced; `test_p25_mid_columns_across_steps`: 3,797 pairs of different pitch codes (steps of up to 205 codes up and 206 down), the pair 134 and 161 both ways among them, four voicing words each, 60,752 frame pairs, against the rule stated on the words, each frame's nearest harmonic in double with ties upward; `test_dstar_exact_column_edge`: a 250 Hz model whose harmonic 2 lies exactly on the 500 Hz edge, in the upper column; `test_rate33_hold_past_l`: rate 33's kept 505 Hz hold against a reference in the library's float arithmetic, a harmonic past L read at L, at 100 Hz with L 1, at the pitch whose bits are 0x3e0761b7 (the float quotient exactly 3, where a double one reads harmonic 4; not compared in builds that evaluate float wider than float, where 32-bit x87 Release reads harmonic 4 too) and at `FLT_MIN`; `test_p25_forged_tone_in_state`: `MBE_FRAME_MARK_TONE` forged in any of the three states of a stream leaving the samples and all three states exactly as without it, P25 and ProVoice through `mbe_processImbe4400Dataf()` and ProVoice through its own entry point, `mbe_processImbe7100x4400Framef()`, on frames packed for it; the columns' rule reads no mark of the previous model, while `mid_pitch()` and the train's loading still read its `PITCH_CARRIES`, `PITCH_ABSENT` and `PREVIOUS_HELD`); D-STAR's column positions at the ordinary pitch codes 0 to 123 the decoder's at its own f0, at least 6.0e-5 of a column from a whole number and 1.5e-5 from a half (6.58e-5 and 2.00e-5 in IEEE builds, 1.81e-5 with fast-math and x87; `test_dstar_column_edges`; code 124's silence model, at 250 Hz, lies exactly on positions 0.5 and 1.0, a margin of zero, and decodes with no voiced band); the pulse train through a frame lacking a pitch after a frame voiced above its first band running on when a column holds and starting afresh when none does, rendered and advanced alike, D-STAR and P25 (`test_train_through_held_columns`); a band with no state in either frame playing nothing at mid-frame where the columns decide, and pitches and harmonic counts no decoder gives leaving no voiced band (`test_mid_columns_bands_with_no_state`); every mark value and a stored value that is none, in the frame before, this frame or both, in each subframe rendering, the mid-frame voicing the reference's and rendering and advancing alike (`test_marks_with_columns`); every mark value in each of a stream's three states through each codec's data path playing finite samples, exactly as the stream with no mark there in every slot (`test_marks_in_codec_states`), `MBE_FRAME_MARK_MID_COLUMN_0` set in each of a stream's three states after voiced frames and a frame with no voiced band through whose first half the train ran on (rate 33 V/UV row 17, D-STAR b0 34 and P25 and ProVoice b0 25, which lack a pitch), before a voiced frame, directly or after a lost frame, through D-STAR's, P25's and rate 33's data paths and ProVoice's own entry point (`mbe_processImbe7100x4400Framef()`, frames packed for it), against the same state with no mark in that slot: on the last frame played in the rate-33 path, before the voiced frame, the train that frame stored loaded, the samples and the three returned states those of the stream with `MBE_FRAME_MARK_PITCH_ABSENT` there; everywhere else the samples bit for bit those with no mark, the returned states after the voiced frame exactly those, and after a lost frame alike in everything but the mark word, the forged 6 in `prev_mp` carried into all three returned states by D-STAR's, P25's and ProVoice's repeat and kept in rate 33's history, every other mark word the stream's; and in D-STAR's, P25's and ProVoice's paths the train the last frame played stored shown to change the samples when it runs on, so that 6 put over its own `MBE_FRAME_MARK_PITCH_ABSENT` plays as no mark does (`test_row17_mark_in_codec_states`; placed after `test_p25_forged_tone_in_state`), every column edge and every edge 250 Hz lower at least 5.97e-5 of a harmonic from a whole number at each of the 208 pitch codes at the pitch P25 plays (b0 152 at 2005 Hz), and 1.55e-4 at ProVoice's nominal one (b0 36 at 2755 Hz) (`test_p25_onset_edges`), an onset frame muted by the error rate leaving the pulse train and the noise generator where the same frame heard leaves them but no unvoiced overlap, which on the test's input (an onset with unvoiced bands, and the same error history after it in both streams) makes the frame after differ in its first 80 samples only (`test_p25_muted_onset_advances`), a frame that decodes with no voiced band after frames with none, whose adaptive smoothing voices harmonics, rendered as an onset read in the columns, decoded, repeated and muted alike (`test_p25_smoothed_onset`: one harmonic voiced covers no column and starts the train at mid-frame, harmonics 5 to 12 cover one and start it with the frame; the history keeps the decoded model and mark, the model rendered carries `MBE_FRAME_MARK_ONSET_COLUMNS`), each mark of the mid-frame voicing (5, 6, 7 and the late onset) acting in its own codec's rendering only, while the other marks act as marks in every rendering (`test_marks_act_in_their_own_rendering`: in D-STAR's rendering 5, 6 and 7 play as no mark and a late onset acts after a frame with no voiced band only; in the rate-33 rendering 5 plays as none; on the previous model, a frame with no voiced band through whose first half the train ran on, 6 renders the next frame as `MBE_FRAME_MARK_PITCH_ABSENT` and `MBE_FRAME_MARK_PREVIOUS_HELD` do in the rate-33 rendering and as no mark in D-STAR's and P25's, the two told apart in each), P25's rendering marking the model it renders (`test_p25_rendering_marks_the_model_it_renders`: on a voiced model every mark but a tone model's `MBE_FRAME_MARK_TONE`, the late onset and the held previous model included, rendering as the onset mark and replaced by it; a model with no voiced band keeping any other mark and losing the onset mark), rate 33's marks in a caller's P25 history or last frame played playing as none before a rejected frame and the frame after it (`test_p25_foreign_marks_in_state`), and the marks that shape the mid-frame voicing on a caller's state, with any pitch bit pattern in either frame (none, subnormal, the largest float, just under and past Nyquist), harmonic counts 1 and 56 and voicing values no decoder gives, rendering finite samples and state, rendered or advanced, from silence or not, under the sanitizer build too (`test_untrusted_mid_voicing_marks`), and the previous model held to mid-frame by a rate-33 silence frame; a frame that starts from silence: the independent rendering of the waveform test through such frames, a steady voiced pair rising as n / 80 of the ordinary rendering from an empty train, stale pulses, phases and overlap neither sounding nor moving the new train, a silence frame rising to its own model, and advancing leaving the pulse train and the noise generator where rendering leaves them (`test_subframe_synthesis`)
- the pulsed state of D-STAR and P25 (`test_subframe_synthesis`): the waveform test's independent
  rendering has the rule from its statement: the 56 patterns typed from the published table, each
  code found from the pitch the codec's public decoder gives it, compared by bit pattern; the bands
  by their 500 Hz channel and that of 250 Hz more, the frequency l w0 4000 / pi in double as the rule
  states it (no channel edge avoided); the alternation of the pulses; the pulse 24 samples after its
  model's point through the two models' pulsed bands within one interval; the voiced response with
  the unvoiced factor taken out, its window and the low-band roll-off's gain and phase in this test's
  own double arithmetic, a harmonic rendered only below Nyquist as the library compares it (the
  product l w0 in float); clause 5's mid-frame pitch; and P25's validity counter and read flag
  (`reference_read_codes()`, `reference_next_count()`). Its pieces from the library: in the waveform
  tests the three models as `mbe_subframe_own_model()` and `mbe_subframe_mid_model()` build them,
  the two own models whole (pitch, harmonic count, voicing and magnitudes) and the mid-frame model's
  magnitudes, the mid-frame model's pitch and harmonic count checked against the test's own and
  replaced by them, and its voicing too but on a frame marked `MBE_FRAME_MARK_PREVIOUS_HELD` (where
  the library's copy of the frame before's own model stands), the phases the test's own
  (`reference_own_phases()`, `reference_mid_phases()`); in the data-path tests the models' pitch,
  harmonic count, magnitudes and regenerated phases as the library builds them (the mid-frame
  model's pitch and harmonic count checked against the test's), the unvoiced synthesis of the bands
  left unpulsed (`mbe_synthesizeUnvoicedSubframes()`, DVSI's generator), and, where a test says so,
  the spectral enhancement, the repeat synthesis and the tone renderer. It plays a stream through the state
  (`test_waveform_through_pulsed_codes`: frames with no voiced band without a code and with codes
  whose own, mid-frame or both patterns pulse, in every order of two and three (a de Bruijn
  sequence), default-code frames, voiced frames before and after each kind, tone models before and
  after frames with codes and a caller's tone mark on a frame with a code, frames from silence, and
  P25's counter walked to its floor and its ceiling), sample for sample in D-STAR's, P25's, rate
  33's and ProVoice's renderings, the last two playing no pulse at the same pitches, every frame's
  stored state compared with the test's own constants (slots 34, 35 and 36, `1.0f` by its bits);
  D-STAR's pulsed bands at a channel's edge (the 3rd harmonic at 749.999964 Hz, which a float
  product would put on the edge) and, in IEEE builds, at Nyquist (`test_pulsed_channel_edges`);
  every one of the 56 codes of each codec through its data path
  (`test_pulsed_codes_through_data_paths`: a stream of the default code, the reserved code three
  times, the default code, the code and the default code, every frame checked (D-STAR's stream first
  plays one unchecked frame, after which its first call's replacement of the initial state is
  behind it), amplitude fields random, P25's built with the IMBE quantizer, each frame against the
  rule's pulses plus the library's unvoiced synthesis (DVSI's generator) of the bands the reference
  leaves unpulsed and the state each frame stores; ProVoice, through its own entry point, fed the
  bits the P25 stream played, frame for frame, and playing them with no pulsed band); P25's counter
  through its data path (`test_p25_validity_counter_walk`: a fresh stream reading its first code,
  frames at b0 127, the code this library's P25 encoder sends such frames, taking it to -10 and
  holding it there, the next ten code frames not read and the eleventh read, the ceiling of 3,
  ten frames at the default code b0 25 bringing it from -10 to 0, a voiced frame leaving it, and
  D-STAR reading every code after twelve frames outside its table); the counter through P25's
  repeats, its repeat limit's mute and the error rate's (`test_p25_counter_through_a_mute`: every
  one of fifteen lost frames checked for its repeat count, mute flag, model, counter and read flag,
  losses 1 to 4 replaying the decoded model, the 4th muted, loss 5 installing the reset model with a
  repeat count of 0 and the mute flag clear, mutes again at losses 9 and 14; after three code-48
  frames the counter stands at 3 through the run and the next code frame is read, after twelve
  frames at b0 127 it stands at -10 and the next code frame is not read and raises it to -9; three
  frames muted by the error rate at b0 127, advanced unvoiced, lowering it from 1 to -2 and the
  next two code frames not read); the counter's steps from every value it holds for a voiced model,
  a tone model, the default code, a reserved code, b0 127's pitch, a pitch an ulp off a reserved
  code, the two floats beside the default code's pitch and the reset model's pitch, rendered and
  advanced alike (`test_pulsed_count_steps`); the frame
  before's read flag deciding its own model's pulse and clause 5 (`test_pulsed_read_flag`); the
  three slots read by their bits (`test_untrusted_pulsed_state`: 59 bit patterns, not a number,
  infinities, subnormals, the zeros, the largest floats, 0.5, 1.0 and its neighbours, every whole
  number the counter holds and both its neighbours (the test checks that the list holds them) and
  whole numbers past its range, in each slot, a
  frame after each playing and storing what it does after the value the slot reads as; in D-STAR's
  rendering the read flag and the counter not read; a state without the subframe mark read as a
  fresh stream's); no pulsed state in rate 33's, ProVoice's or the TIA rendering's subframes
  (`test_pulsed_renderings`); every pitch code, P25 0..207 and D-STAR 0..123, decoded by the public
  decoders with no voiced band (the test requires all 208 and all 124 to decode so), each ordinary
  one playing no pattern and each reserved one its own, and P25's counter stepped by each as the
  rule says, -1 for an ordinary code and +1 for b0 25 and the 56, D-STAR's left at 0
  (`test_ordinary_codes_carry_no_pattern`); what has a code
  (`test_pulsed_state_is_the_models`: pitches one ulp from each code's having none, a voiced band
  taking it away, a band with no state left unpulsed, every mark of either frame leaving the pulses
  alone but the two that put another model at mid-frame, `MBE_FRAME_MARK_PREVIOUS_HELD` and
  `MBE_FRAME_MARK_TONE`); a caller's `MBE_FRAME_MARK_PITCH_ABSENT` holding voiced bands at
  mid-frame on a frame with a code, whose mid-frame model then takes no pattern
  (`test_pulsed_held_by_a_forged_mark`); a muted frame storing the emission as rendering it does
  (`test_pulsed_advance_stores_the_emission`); adaptive smoothing that voices a band taking the code,
  the frame storing no emission, no read code and the counter it found, rendered and, muted by the
  error rate, advanced (`test_pulsed_smoothing_takes_the_code`); D-STAR's lost frames, its repeat limit's mute and the
  voice after it (the resume from the prediction history at 0.3 of its amplitudes rebuilt by the
  test, with the library's spectral enhancement and its unvoiced overlap from the library's repeat
  synthesis, and compared wholly; the two frames after wholly), its tone's first frame (the repeat of the last frame fading out under the tone over the
  first 80 samples) and the voice frame after the tone (fading in under the tone's tail), the tone
  and its tail from the library's tone renderer, and P25's repeats, each frame against the rule
  within one stream (`test_pulsed_repeats_mutes_and_tones`); and ProVoice after a caller's model at
  a P25 reserved pitch, with each of the 59 bit patterns in each of the three slots and each of 13
  mark values (`test_provoice_has_no_pulsed_state`). Mutations checked to fail them, each pulsed-state test run alone (in the whole suite
  every one fails `test_subframe_synthesis`), named by the tests that fail (the waveform stream, the
  channel edges, the data paths, the counter's walk, the counter's steps, the mute test, the read
  flag, the slots, the renderings, the ordinary codes, the smoothing test, the repeat and tone test,
  and the rest by their names): the pulse at 23 or 25 samples (waveform D-STAR and P25, edges, data
  paths, walk, repeats); the alternation dropped (waveform, data paths, walk, slots, steps,
  repeats) or inverted (all but the waveform's rate-33 and ProVoice runs, the mute test, D-STAR's
  slots, the smoothing, the forged-mark and the ProVoice tests); a voiced model counting as having
  emitted (waveform, all four renderings); the pulse ringing on past the next model's point
  (waveform, data paths, repeats); the voiced train through pulsed bands, or the pulse through voiced
  bands (waveform D-STAR and P25); pulsed bands also playing noise (data paths, walk, repeats); the
  reach of 0 or 500 Hz (waveform, edges, data paths, repeats); the two patterns swapped (waveform,
  edges, data paths, read flag, `test_pulsed_state_is_the_models`, repeats); the table shifted by
  one code (waveform, edges, data paths, walk, repeats); D-STAR taking P25's mid-frame pitch (its
  waveform, edges, data paths, read flag, repeats) or P25 D-STAR's (its waveform, data paths, walk,
  read flag, repeats); `mbe_advance_subframes()` not storing the emission (steps,
  `test_pulsed_advance_stores_the_emission`); ProVoice taking P25's pulsed state (its waveform,
  renderings, the ProVoice test); the counter ignored, or read at -1 (P25's waveform, walk, the mute
  test, slots, steps, read flag); its ceiling 2 or 4 (P25's waveform, data paths, walk, the mute
  test, slots, steps, repeats); its floor -9 or -11 (walk, the mute test, slots, steps); b0 25 not
  counted (P25's waveform, data paths, walk, steps, ordinary codes,
  `test_pulsed_advance_stores_the_emission`, repeats); a voiced frame lowering it (P25's waveform,
  walk, steps, smoothing) or raising it (walk, steps, smoothing); a model at no code's pitch
  lowering it (P25's waveform, the mute test, steps); the counter applied to D-STAR (D-STAR's
  waveform, edges, data paths, walk, slots, steps, read flag, ordinary codes,
  `test_pulsed_advance_stores_the_emission`, repeats); the read flag of the frame before ignored
  (P25's waveform, walk, slots, read flag); an ordinary code (b0 12) leaving the counter (ordinary
  codes); the float just above b0 25's pitch counting as b0 25 (steps); a forged counter word
  (`0xC0E00001`) read as -7 (slots); a reserved frame the smoothing voices still read and raising the
  counter (smoothing); the first 80 samples of D-STAR's voice after the repeat limit's mute erased
  after a b0 69 frame (repeats); the outgoing repeated speech erased when a tone follows a b0 69 frame
  (repeats); a harmonic's frequency formed in float (edges); P25 b0 1 or b0 12 taking pattern 48
  (ordinary codes); "emitted" stored as `3.0f` (waveform, edges, data paths, walk, slots, advance,
  repeats); the TIA rendering falling through to D-STAR's codes (renderings); the Nyquist comparison
  in double (edges, IEEE builds); the low-band roll-off's corner at 101 Hz (waveform, all four
  renderings, edges, data paths, walk, repeats, and six other test programs). Five of them (b0 12 leaving the
  counter, the float above b0 25's pitch, the forged word, the smoothed frame, the erased resume)
  passed the tests as they stood at `6a97d51`. A scratch harness outside the
  repository, built against the baseline's and the head's libraries in each of the ten local
  configurations, compares rate 33's and ProVoice's data paths, `mbe_synthesizeSpeechf()` and the
  three encoders byte for byte over streams whose caller states hold reserved pitches, every mark
  and values in the slots ([DVSI parity](dvsi-parity.md#november-2026-the-pulsed-state)).
- rate 33's pulsed rows (`test_subframe_synthesis`; [DVSI parity](dvsi-parity.md#rate-33-the-pulsed-rows)):
  the waveform test's independent rendering has the rule from its statement: the fourteen masks of
  V/UV rows 18 to 31 typed as the document prints them, the word a frame carries in `PHIl[0]` and
  the one it stores in slot 37 formed in the test (`test_row_word()`), the frame's model checked
  against the band-grid model of the word's code by bit pattern, the bands by their 500 Hz channel
  and that of 250 Hz more in double, the pulse times computed in thousandths of a sample and rounded
  to float once (`reference_row_time()`), and each pulse rendered through its own model's pulsed
  bands only, under that model's crossfade weight, rising in the interval that ends at its point and
  falling in the one that starts there, from its time or the interval's start
  (`reference_row_render()`), with the voiced response of the pulsed state's reference (its window,
  the unvoiced factor taken out in the library's arithmetic, `reference_pulsed_magnitude()`: the
  factor in float rounded to float, the quotient in double rounded to float and capped at `1e30f`,
  which a fast-math build of the library may evaluate otherwise; the low-band roll-off in the test's
  own arithmetic). What it takes
  from the library for the rows: the band-grid model of each code (its w0 and L) as the public
  decoder, `mbe_decodeAmbe2450Parms()`, gives it, and the three models as the pulsed state's
  reference takes them (`mbe_subframe_own_model()` and `mbe_subframe_mid_model()`, the mid-frame
  model's pitch, harmonic count and voicing checked against the test's own). It plays every row at
  every code 0..119 (all 1,680 pairs, so that each pair's playback, not only its word, is checked;
  rows 24 and 25 at every pair of grid indices that occurs) between frames with no voiced band, so
  that pulses fall before their model's window, in its rising half, on either side of its point (b0
  79 puts a single-row pulse 0.626 sample before it, b0 80 0.390 after it; no legal code puts a pulse
  exactly on a model's point), in its falling half and (own models) in the next frame, L 9 and 56
  among them; plain, own-pulsed, mid-pulsed and
  both-pulsed frames in every order of two and three (a de Bruijn sequence); voiced frames before and
  after each kind; silence frames after pulsed ones; repeats (held words) of each kind, twice for
  own-only and both-pulsed rows, and the frame after; frames from silence of each kind and after
  pulsed frames; and words on models with a voiced band
  (`test_waveform_through_pulsed_rows`), sample for sample in rate 33's rendering and in D-STAR's,
  P25's and ProVoice's, which play no pulsed row at the same frames, every frame's stored slot
  37 compared with the test's own word. Through the data path and the decoder: every row 0..31 at
  every code 0..119 decoded by the public decoder and through `mbe_processAmbe2450Dataf()` from a
  fresh state, a word in all three returned states exactly for rows 18..31, none for silence
  frames, erasures, a tone frame or a rejected repeat from a fresh state, and a rejected repeat
  after a word frame holding the repeated frame's held word (`test_r33_rows_decoder_word`); the
  repeat of rows 19, 24 and 27, twice, and the frame after, against the same stream with the word
  cleared (`test_r33_rows_repeat`); the repeat under a tone's first frame, rows 20 and 27, its
  first 80 samples differing from the stream without the word by 1 - n / 80 of the difference
  between the two repeats (each rendered with the library's repeat loader and its repeated-speech
  synthesis, unlimited) (`test_r33_rows_under_tone`); a frame advanced unrendered storing what rendering stores
  (`test_r33_rows_advance`); adaptive smoothing that voices a band making the frame no pulsed-row
  frame, rendered and, muted by the error rate, advanced (`test_r33_rows_smoothing`); fully pulsed
  frames playing no noise (`test_r33_pulsed_bands_play_no_noise`, with and without the unvoiced
  synthesis); the 2,640 values that are words (the decoder's 1,680, rows 18..31 not held, and the
  repeat's 960, rows 18..25 held; a held form of rows 26..31 is none) and both float neighbours of
  each, 78 bit patterns of every class (`k_slot_values` and `k_row_values`, 77 distinct, three of
  them words; 4088.0f, held row 25 with the code 120, and 4096.0f, held row 26 code 0, among the
  others), words of every row, code and held form in `PHIl[0]` on a grid model and in slot 37 of
  the last frame's state, a word on a model with a voiced band (harmonic 1 only, its mid-frame model
  unvoiced), on a pitch an ulp off the grid, on the grid pitch with another harmonic count, of a
  code of another L, a held word of row 28, the stored word with the last model an ulp off its grid
  pitch, and a state without the subframe mark, each playing and storing as the state with no word,
  or, taken, playing otherwise than it and storing in slot 37 its word for rows 18..25 and 0 for
  rows 26..31, whose own mask is empty (`test_untrusted_row_words`); the
  repeat loader leaving a caller's word on a model that is not its code's grid model (a word of row
  20, one of row 28, a held word of row 20, and row 20's word of code 40 on code 0's grid model)
  through `mbe_processAmbe2450Dataf()`: an erasure returning the branch point's `PHIl[0]` in all
  three states (the value in `cur_mp` and `prev_mp_enhanced`, 0 in `prev_mp`), and the samples and
  the states but for that word of the same stream without the value
  (`test_r33_repeat_leaves_foreign_words`); the scope of that guarantee: on the grid model of code 40
  with harmonic 1 voiced, a word of row 20 or of row 28, then three repeats and the repeat limit's
  mute, the loader, which tests the grid and not voicing, returning the held form 3368.0f in
  `cur_mp` and `prev_mp_enhanced` (row 20) or clearing the word (row 28) where the branch point
  returned the value, while the samples, slot 37 and every other field of the three states are
  those of the stream without the word, its mute noise seeded alike
  (`test_r33_repeat_on_a_voiced_grid_model`); what is not
  promised, as the limit of the representation: a caller's word on the grid model of its code plays
  as the history of a frame of that row, in `PHIl[0]` before an erasure and in slot 37 before a frame
  of row 16 (`test_r33_forged_history_plays`); a pulsed band at the cap, every `Ml` at 1e30, against
  the independent rendering at the capped magnitude (`test_r33_pulsed_magnitude_cap`); and every
  pair's pulse times by their bits, a float replica of the library's pulse rendering fed the times
  in double rounded once reproducing the library's frame and the frame after it sample for sample
  (`test_r33_pulse_times_by_bits`): the times in float differ at 132 of rows 24 and 25's 480
  evaluations (and the single law's at 212 of 360). It compares bit for bit where the test can
  reproduce the library's arithmetic (`PULSE_REPLICA_EXACT`: the library built without fast-math,
  float expressions in float or, unoptimized, stored at each statement as the library stores them,
  and no contracted multiply-adds where the target fuses them, which CMake states for the private
  tests with `MBELIB_TEST_FP_CONTRACT_OFF` in GCC's and Clang's Debug builds), and elsewhere within
  a thousandth of the two frames' peak, a tolerance that cannot tell the times formed in float from
  those formed in double. Of the ten local configurations, dev-debug, asan-ubsan-debug,
  notones-debug, the Debug build with SIMD, the Debug `-march=x86-64-v3` build and the 32-bit Debug
  build compare bit for bit, and the float-times mutation fails in each; dev-release, the two Clang
  fast-math LTO builds and the 32-bit Release build take the tolerance, and the mutation passes in
  each: it is caught only in the former six. P25's slots 34 to 36 not read as a row by rate 33's
  rendering, and slot 37 not read by D-STAR's, P25's or ProVoice's rendering, each writing 0 there
  (`test_row_slot_between_renderings`); and words in each of the three states' `PHIl[0]` and in slot
  37 of D-STAR's, P25's and ProVoice's streams through their own entry points, then a lost frame and
  a voiced frame, playing as the same state with no mark there, the returned states alike but for
  the mark word (`test_row_words_in_codec_states`). A word of another code with the model's L is
  taken, with that code's times: the band grid's pitch depends on L alone. Mutations checked to
  fail them, in throwaway builds outside the repository (`test_subframe_synthesis`, and
  `test_golden_pcm` where marked G; the test named is the one that fails first in the suite's
  order): the slope 1.0 and 1.03, every single-row pulse one sample early and late, rows 24 and 25
  with `mod` and `div` swapped, with a grid step of 11 and of 12 samples and a grid of 12 positions,
  on the single law, row 23 on the grid law, the own and mid-frame masks swapped (G), the table
  shifted by one row (G), one channel of row 21 flipped, alternation (G), a pulse played through
  the neighbouring model too, a pulse cut at its model's point (G), a pulse ringing past its model's
  window, the carried pulse dropped at the frame boundary (G), the voiced train through pulsed bands
  (G), the reach 0 and 500 Hz (G), the word honoured in D-STAR's, P25's and ProVoice's renderings,
  the repeat literal, the repeat's mid-frame pulse missing, the last frame's pulse played from
  silence, a smoothed-voiced frame still pulsed (the waveform stream: its stored word, with words on
  voiced models) (`test_waveform_through_pulsed_rows`); `mbe_advance_subframes()` not storing the
  slot (`test_r33_rows_advance`); pulsed bands also playing noise (G;
  `test_r33_pulsed_bands_play_no_noise` when the rows' tests run first, the pulsed state's data
  paths in the suite's order); a word accepted on a model of another L, on an off-grid pitch, on a
  model with a voiced band, and a neighbour of the stored word accepted (`test_untrusted_row_words`);
  the repeat dropping the state (`test_r33_rows_decoder_word`); a repeat of rows 26..31 pulsed
  (`test_r33_rows_repeat`); and, found surviving by the review and caught since, the harmonic-count
  test dropped (`test_untrusted_row_words`, the grid pitch with another L), row 21 at code 41
  ignored (`test_waveform_through_pulsed_rows`, every pair), the last model's pitch test dropped
  (`test_untrusted_row_words`, the stored word an ulp off), the grid times computed in float
  (`test_r33_pulse_times_by_bits`) and the pulsed magnitude's cap removed
  (`test_r33_pulsed_magnitude_cap`); and the fix round's own two, the repeat loader acting on a word
  on any model (`test_r33_repeat_leaves_foreign_words`) and the reader taking a held form of rows
  26..31 (`test_untrusted_row_words`), the latter as the reader was before the fix (the range
  bound back at 4855, the last value of the old form, and the row test of a held word removed): the
  bound at 4087 and the row test each reject those values alone, so restoring the old bound alone,
  or removing the row test alone, leaves the reader as it is and survives, an equivalent mutation
  (removing the bound altogether was not run: it would accept 5120.0f as row 18 code 0 with a held
  field of 2 and give out-of-range floats to the integer conversion); and the second fix round's
  own, the loader testing voicing too (`test_r33_repeat_on_a_voiced_grid_model`), which the loader
  acting on any model leaves passing (on the grid model it acts as the library does) and the
  foreign-word test catches. The word accepted on a model with a voiced band passed the tests until
  the case of a model voiced at harmonic 1 alone was added (`pulse_bands()` refuses a model with a
  voiced band itself, so only a mid-frame model left unvoiced shows it). The shared pulsed state's
  mutations (the voiced train through pulsed bands, the reach, pulsed bands playing noise) fail the
  D-STAR and P25 tests first in the suite's order and the rows' tests when those run first. A scratch
  harness outside the repository, built against the baseline's and the head's libraries in each of
  the ten local configurations, compares D-STAR's, P25's, rate 33's and ProVoice's paths,
  `mbe_synthesizeSpeechf()` and the three encoders byte for byte over streams whose caller states
  hold reserved pitches, every mark, pulsed-row words and values in the slots, rate 33's as well as
  the other paths' (rate 33's own stream holding no frame of rows 18..31; its first form gave rate
  33's states no word and did not say so, which hid the repeat loader's defect) ([DVSI
  parity](dvsi-parity.md#november-2026-rate-33s-pulsed-rows)).
- DVSI's low band (`test_subframe_synthesis`, `test_unvoiced_subframes`, `test_params`, `test_p25_tones`;
  [DVSI parity](dvsi-parity.md#dvsis-low-band)): the fixed third-order high-pass at 75 Hz against its
  factored form in double, gain and phase within 2e-7 at every bin of the unvoiced transform and on a
  sweep of harmonic frequencies, 1/sqrt(2) and 3 pi / 4 at 75 Hz, and every class of float by its bits
  (`test_low_band_dvsi`): 0 for both at +0, -0, subnormals, negative values, -infinity and NaNs of
  either sign, 1 and 0 at +infinity, a gain of 0 and the phase's limit 3 pi / 2 (to 1e-6) at FLT_MIN,
  and far above the corner, from 75 kHz to FLT_MAX across the switch to the uncancelled form at 76.8
  kHz, the phase against its series in 75 / f to 2e-7 of its value and the gain exactly 1; the
  second-order filter's own tests unchanged; the unvoiced subframes under DVSI's filter with each
  subframe's corner (160 and 190 Hz) against the equations, the bins under the corner zero after the
  band's scale is taken from its whole power, a corner exactly at a bin's centre (62.5 and 187.5 Hz)
  keeping that bin, which the float above it removes (each by more than four times the comparison's
  tolerance), and bands above both corners playing the same bits cut and uncut
  (`test_dvsi_subframes_match_equations`); the tracker's rule (`test_low_cut_tracker`, in D-STAR's,
  P25's and rate 33's renderings, none in ProVoice's): 110 Hz on a fresh state and as the first
  frame's corner, then fifteen kinds of model (fully voiced; a top, bottom or middle band unvoiced; a
  band with no state; no voiced band; a tone mark; a silence mark on a fully voiced model; the AMBE
  silence model fully voiced; the silence pitch with the other codec's harmonic count; a reserved
  code with and without every band voiced; every band voiced at a subnormal and at an infinite pitch,
  no pitch by the library's class), each moving it or not as the rule says in that rendering
  and each advanced (`mbe_advance_subframes()`) storing what rendering stores, the stored corner the
  tracker before capped at 200 Hz, the ceiling while a 252 Hz stream takes the tracker past 240 Hz and
  the plateau coming down to 180 Hz, a caller's pitch keeping the tracker at 500 and 50 Hz; the range
  and the silence models (`test_low_cut_range`: every pitch code of the three codecs, the pulsed rows'
  grid and every D-STAR voicing included, decoding to 50 to 500 Hz; D-STAR's b0 124 at the silence
  pitch with 14 harmonics, rate 33's b0 124 and 125 with 15 and the silence mark, and no voice code of
  either at the silence pitch whatever its harmonic count); the knee (`test_low_cut_knee`); the cut's level (`test_low_cut_level`); the unvoiced
  cut through the rendering (`test_low_cut_unvoiced_frame`); the slots as caller state
  (`test_untrusted_low_cut`: 24 bit patterns in each slot, each playing and storing as the value it
  reads as, and every value read, the bounds 50 and 500 (tracker), 50 and 200 (corner) among them,
  playing otherwise than the value a rejected one falls back to; a state without the subframe mark and
  one the TIA rendering wrote reading as a fresh one; ProVoice reading neither and storing 0 in both)
  and moved between renderings (`test_low_cut_between_renderings`); through the public data paths,
  one test per row of the table of frame kinds in [DVSI parity](dvsi-parity.md#dvsis-low-band), each
  frame's stored tracker against the rule computed in the test on the model the frame renders (or the
  tracker before, by bits, for one that holds it): D-STAR's voice, lost and repeat-limit frames, its
  error-level mute, the resume, tone entry after an audible and an inaudible frame, tone tails
  (`test_low_cut_dstar_path`), b0 124 and a reserved code voiced throughout by adaptive smoothing and
  a caller's silence mark on the history (`test_low_cut_dstar_smoothed`); P25's voice, mute, repeats,
  reset model at five error rates, a reserved code voiced throughout and a caller's silence mark
  (`test_low_cut_p25_path`) and its tone model (`test_tone_model`); rate 33's voice, silence, repeat,
  mute and tone frames (`test_low_cut_r33_path`), one to three valid frames the error rate mutes, a
  muted invalid frame after them and after a heard frame, tone entry at an ordinary and at a muting
  error rate (`test_low_cut_r33_mutes`), a silence frame voiced throughout heard and muted, its repeat
  and the fade after it, a pulsed row voiced throughout and a caller's silence mark
  (`test_low_cut_r33_smoothed`, the review's reproductions among them); the table's conditional rows
  on the histories that hold the tracker (`test_low_cut_conditional_rows`: D-STAR's tone entry and
  single-tone tails after a frame with no voiced band, a caller's silence mark through D-STAR's fade
  and muted repeat, rate 33's muted row-16 frames, a fade that drops a caller's mark, the lost frame
  after a muted tone with the error rate still high, the repeat limit under a held tone); and caller
  states of every float class of pitch at harmonic counts 1 and 56 playing finite samples
  (`test_low_cut_caller_pitches`). The independent rendering of the waveform tests (`run_waveform()`)
  renders the low band in its own arithmetic in all four subframe renderings (the filters from their
  factors at the one float frequency the library forms, the cut at each model's corner, the knee at
  that frequency) and carries the tracker as the library does, a float rounded once a frame; each
  frame's stored tracker is checked as one step of the rule from the tracker the frame before left,
  and its stored corner by bits. The tracker comparisons (`track_is()`) are bit for bit in the builds
  where the test forms the step as the library does (`MID_PITCH_INEXACT` 0: dev-debug,
  asan-ubsan-debug, notones-debug and the Debug build with SIMD), and elsewhere (dev-release and both
  Clang LTO builds with fast-math, the Debug `-march=x86-64-v3` build, which fuses multiply-adds, and
  both 32-bit builds, whose evaluation method is 2) the stored value must be the exact step rounded to
  nearest, within half its spacing and 1e-7 of it, which a value one ulp off fails unless the
  exact step lies within 1e-7 of a spacing from a tie (`TRACK_TIE_BAND`), where either float passes:
  every compared step is checked for that, and none of the file's 741 steps lies there
  (`test_low_cut_ties`, which also takes a case inside the band: T 110 Hz, pitch bits 0x3d8504c0, 3.9e-8
  of a spacing from the tie); what that bound cannot catch is a step formed
  otherwise that rounds to the same float on every frame of the streams (none found: the step formed
  in float rather than double and the tracker an ulp off fail both in the Debug build and in the Debug
  `-march=x86-64-v3` build, where the bound applies, and so does the pitch's class weakened to
  "greater than zero"). Mutations, each failing a test (`test_subframe_synthesis` unless named; the table, run
  in the Debug build, is in the change's report): the power 1 and 3; a Butterworth magnitude for the
  knee; the ceiling removed or on the tracker; a first value of 100 and of 120; the rate halved and
  doubled; the tracker moved by a frame with one band unvoiced, by one with no voiced band, by a tone
  model, by a silence mark, by the AMBE silence model; the silence model told by its pitch alone, or
  given to P25's rendering too; not moved by `mbe_advance_subframes()`; stepped twice; the corner read
  after the frame's own step; the cut on the mid-frame model only, applied to ProVoice, to the TIA
  rendering (`test_params`) or to the stored model; the third-order filter in ProVoice; the old filter
  in rate 33's, P25's or D-STAR's rendering, or in P25's unvoiced bins alone; the last frame's own model
  at this frame's corner and the mid-frame model at the last frame's; the stored tracker read without
  its range or without the subframe mark; stored corners of exactly 50 and of exactly 200 Hz
  rejected; the pulsed bands not cut; the tracker not kept in range, an ulp off or its step formed in float;
  the pitch's class weakened to "greater than zero"; the unvoiced bins
  scaled instead of zeroed, zeroed under 0.9 of the corner, the bin at the corner removed, the old
  table in DVSI's bins, the new one in ProVoice's (`test_unvoiced_subframes`); the unvoiced cut on the
  mid-frame subframe only; the filter's gain at +infinity, subnormals and NaNs taken as frequencies,
  the phase formed directly far above the corner (`test_params`); the tones' start phase on the old
  filter (`test_ambe2400_frame_types`, `test_ambe2450_frame_types`, `test_frame_paths`). No encoder
  synthesizes a model, so a cut applied to the stored model changes no encoder's bits; it is caught by
  the state the rendering returns.
- Four rules of DVSI's P25 decoder, and the two of its pulse train (`test_params`,
  `test_subframe_synthesis`, `test_frame_paths`, `test_p25_tones`, `test_p25_rendering_copies`; [DVSI
  parity](dvsi-parity.md#november-2026-four-rules-of-dvsis-p25-decoder)): the pitch of every P25
  code, all 208, by bit pattern against `2 pi floor(2^21 / (2 b0 + 79)) / 2^19` formed in the test in
  double and rounded once, at or below eq 46 by less than a grid step, L pinned against a table of
  the 208 harmonic counts and K from it, through the public decoder and the data path, ProVoice's
  nominal pitch unchanged (`test_imbe_pitch_law`); the pulsed state's code tests on the grid
  (`test_p25_code_tests_on_the_grid`); the onset columns' edges (`test_p25_onset_edges`: at least
  5.97e-5 of a harmonic from a whole number at the grid pitch, 1.55e-4 at ProVoice's nominal one).
  P25's mid-frame pitch (`test_p25_mid_pitch_pairs`) against the floored quadratic with the fitted
  excess formed in the test from its statement (the mantissa by halving and doubling, not
  `frexp()`) over every ordered pair of the 208 codes, each code against itself included, bit for bit
  in every build, on the grid; on every pair D-STAR's and ProVoice's renderings the geometric mean
  in float and rate 33's its own (the wave against the test's reference rounded once, by its bits
  where the build's double arithmetic is the test's, within an ulp elsewhere; the geometric mean at
  equal codes); its range: pitches that are no positive normal float (zero, negative zero, the
  smallest and largest subnormals, both infinities, both NaNs, negative pitches) in either model
  giving no pitch and no voiced band, the largest pair of equal floats whose fundamentals' product is
  at most 2^64 giving the quadratic's pitch and the next none, two models at w0 2^70 none (the
  review's pair, which the first form gave a pitch), pairs at `FLT_MAX` none, 2^70 with its inverse
  the quadratic's pitch, a product of 2^-37 a pitch and one of 2^-40 none. The pulse train
  (`test_p25_train_by_steps`, run first): from a stored phase of exactly 0, the phase after 12
  frames at each of the 37 codes DVSI's pulse-train vectors hold long enough to read, in whole steps
  of 2^-19 cycle, equal to `12 x 160 floor((n + M) / 2)` modulo a cycle and `floor((n + M) / 2) - n`
  equal to DVSI's (pinned), rendered and advanced unrendered alike, any code that misses named on
  stderr; codes 80 and 81, 53 and 54, 21 and 22 alternating, odd and even sums both met; an interval
  off the grid (code 21's pitch an ulp up) and one across an octave (code 21 then 81, the second
  marked as a tone model so that its own model stands at mid-frame) against the test's linear sweep
  (`reference_advance()`) within 1e-7 cycle (the stored phase is a float) and not the floored
  count; D-STAR's, rate 33's and ProVoice's renderings at 80 then 81 and at 80 then 81 marked as a
  tone model (every end on the grid, the first sum odd: a floor would move the phase by 7.6e-5)
  keeping their linear sweep. The phase source (`test_p25_phase_source`): a P25 own model's phases
  against the library's kernel run in the test on `log2Ml` plus `-log2(MBE_UNVOICED_MAGNITUDE_FACTOR /
  sqrt(w0))` on the bands held as unvoiced, a voiced band and a band with no state as decoded;
  compared by bits where the library is built without fast-math (dev-debug, asan-ubsan-debug,
  notones-debug, the SIMD and `x86-64-v3` Debug builds and both 32-bit builds) and within 1e-6 rad
  where it is built with it (dev-release and the two Clang LTO builds), where the tolerance would not
  catch a step formed in float instead of double; the order of the step and the floor: an unvoiced
  band at -20.1 (b0 140) that the step lifts across the floor, against the same band floored first
  (it differs), one at -40 that stays under, and a NaN, a negative NaN and both infinities, in an
  unvoiced and in a voiced band, entering at the floor exactly; the tail following an unvoiced band
  L; the other renderings' own models the plain kernel's; the step's response at each distance from
  the voicing boundary; the mid-frame envelope; a tone model as the decoder plays one (DTMF 1, b0 53,
  harmonics 4 and 7) read as decoded, and the same structure on b0 140 with code 53's word and with
  the review's word `0x01ff0704`, the mark without the word, other harmonics, another flag and
  another band sounding read as speech. The tone model's test through the public paths
  (`test_p25_tone_model_forgeries`): a stream plays a real DTMF 1 frame, quantized by the IMBE
  encoder's quantizer (the played levels about 201.638 and 201.304), and then a voiced b0 80 frame; each of 23 forgeries of its last frame (the review's, a code with no rows,
  a row of another code, harmonics no row has, the pitch an ulp off, another code's L, another code's
  K, L one less with K kept (21 -> 20), another magnitude at 1e-3, another at -0, a tone harmonic
  unvoiced, voiced at 0 and at a NaN, both tone harmonics at 0 and at -0 (in the droops' ratio, so
  that only the levels' positivity refuses them), the two levels 1 and 100, their ratio off the
  player's by 1e-3 and by 2e-6 (outside the 2^-20 allowed), a `log2Ml` a NaN, +infinity and
  -infinity, the mark alone, the word alone) plays exactly as the same history without mark and word,
  samples and all three returned states, through the P25 data path, ProVoice's frame entry point and
  the P25 data path with ProVoice's flag, and D-STAR's and rate 33's data paths (their own words
  kept); the real state, the ratio off by 2e-7 (inside 2^-20) and a `log2Ml` at the smallest
  subnormal (finite) act in P25's alone, and differ from the history without them. The tone's levels
  decided where it is taken (`test_tone_levels_decided_once`, `test_p25_tones`): a repeat of a played
  DTMF 1 history whose tone magnitudes a caller set to `FLT_MIN`, the others 0, would be played at
  subnormal levels; it plays as the same history with the word cleared, samples and all three states
  alike but for the word, unmarked; over all 116 rows and levels from `FLT_MIN` to 1e30 (three floats
  at each of 47 rungs: `FLT_MIN` times 31^k for k 0 to 45, and 1e30) every model the player marks passes `mbe_imbe4400_tone_model()` and every
  repeat not taken is unmarked (16,008 marked, 348 refused). The own models by columns
  (`test_p25_own_columns`): at every code, with the voicing ending after, starting at and lying only
  at each harmonic, against the reading formed in the test, its own models built in the test
  (`reference_own_model()`, not the library's) and its magnitudes moved across the factor rounded to
  float, as the library forms it, compared by bits without fast-math and within 4e-7 with it (up to
  2.05e-7 seen in dev-release; there a factor formed in double would pass); the forged structure on b0 60 read by the columns. P25's
  copies of the shared code (`test_p25_rendering_copies`, which compiles `src/core/mbe_subframe.c`'s
  static functions into the test): over every ordered pair of the 208 codes, each with one of 64
  pairs of eight voicing classes, the copies give what the shared code gives where P25's rules do not
  apply, band states exactly and floats by their bits (in the 32-bit Release build, where the two
  bodies round apart under x87 excess precision, within 1e-5 with both values finite by their bits,
  1.34e-6 the largest seen) (43,264 pairs, about 1.36 million mid-frame bands by the columns, 676
  fills under the rule); the test compiles the file into its own translation unit with the test's
  flags, so in a library built with fast-math it compares the two bodies as those flags compile them,
  not the library's object, and the static analyzers see it whole (the models' harmonic counts
  bounded before the copies are called). The same test runs the floored train's grid test
  (`p25_grid_steps()`) on every point of the grid, 1 to 2^18 steps with the Nyquist point, each with
  its canonical float giving its count, and on both floats beside each, which fail; and the phase a
  floored interval stores at a whole cycle, tested as stored (6553 steps from 48 2^-19 - 2^-56 reach
  1 - 2^-56 cycle, which leaves exactly 0, and the octave step to 3276 after it launches its pulse at
  sample 80; with the phase tested unrounded the 32-bit Release build stores 1.0 and fails). Its
  comparator fails an infinity or a NaN in either mode, its own bits included, and +0 against -0 in
  the exact mode (`check_comparator()`). The independent rendering of the waveform tests (`run_waveform()`) takes P25's
  mid-frame pitch with the excess, the train's floored means, its phase inputs and its own models'
  voicing from its own arithmetic (`reference_p25_mid_pitch()`, `reference_advance()`,
  `reference_phase_input()`, `reference_own_model()`, `reference_own_voicing()`), over a stream that
  now holds runs of frames at code pitches (80, 81, 53, 54; at the level of the frames before them:
  the tolerance is absolute, and at amplitude 16 a fully voiced frame of 27 harmonics differs by up
  to 0.055 of a sample's 100 or more with or without the floors). The train's ownership of its
  pulses (`test_p25_train_ownership`): the review's stream (codes 172, 65 and 203, every band voiced,
  through `mbe_processImbe4400Dataf()` from a fresh state: the third frame's stored times hold no
  duplicate) and its second reproduction (code 7, the stored phase `0x1.afp-5f`, then code 52);
  chains of floored frames over neighbouring codes from nine stored phases (0; 0.25 and 2^-18, on the
  grid; 2^-19 + 2^-30, 2^-11 of a step past a grid point; the smallest subnormal; 2^-100; 1 - 2^-24,
  31/32 of a step past one; two floats 0.516 and 0.047 of a step past one), the test's reference
  restarting each frame from the phase the library stored, so that it checks each interval's
  ownership for the phase presented to it, not a chain's exact count, each frame's stored pulses equal to the test's own count of the whole numbers crossed
  (`count_crossings()`, one by one in integers), no two stored times nearer than half the shortest
  period, and rendering and `mbe_advance_subframes()` ending alike (count, times and phase by bits);
  the grid's last point (262143 then 262144 steps from phase 0 leaves 524208 steps); the octave
  brackets on the whole numbers: 8736 then 5460 (exactly 8/5) from phase 0.25 floored, the nearest
  pairs on both sides of all four edges, and each edge exactly at 40 pairs, among them pairs whose
  floats' ratio a float test would put inside the bracket.
  Mutations (the third fix round's list, 61, each applied to the head in a copy, the six tests
  `test_params`, `test_subframe_synthesis`, `test_frame_paths`, `test_p25_tones`, `test_golden_pcm`
  and `test_p25_rendering_copies` built in Debug and run, but the phase stored unrounded, built as
  the 32-bit Release build, the only kind that shows it; on `e42a5a1`, the last commit that changes
  the library or a test, all 61 fail a test, the first-failing assertion of each in
  `build/mut/mutations4-final.log` outside the repository; the second round's 59 failed on
  `ecb78d0` and `5369c52` too): the tone taken on its decoded magnitudes alone, without the levels'
  normality (`test_p25_tones`); the phase stored unrounded (`test_p25_rendering_copies`); the grid 2^-18 and rounding in place of the floor; the mid-frame value rounded,
  unfloored, in float, with no excess (fails at b0 51, 53, 75, 80, 94, 95, 164), with its constant
  only (53, 75, 94), with twice the excess, with nodes 1, 1.4 and 2, only where the pitches differ, in
  ProVoice's, D-STAR's (also `test_golden_pcm`) or rate 33's rendering (also `test_params`), with no
  range check and with a range of 2^128; the train's mean unfloored, rounded up, floored only at equal
  codes, left out, in ProVoice's, D-STAR's or rate 33's rendering, across an octave and for pitches off
  the grid; the pulses owned by the rounded time (the first form); the grid's bound tested before n
  (refusing 262144) and after it; the octave decided on the two floats' ratio, and each bracket's edge
  taken inside (down and up); the phase step with its sign reversed, on voiced bands, on bands whose
  log2Ml is not finite, left out of the mid-frame envelope, and leaking into ProVoice's, D-STAR's
  (also `test_golden_pcm`) or rate 33's own phases with the column dispatcher untouched (phase only);
  the last frame's own model read by the bands; the reach 0 and 500 Hz; the first edge 500 Hz; a
  turned band keeping its magnitude; the unvoiced factor formed in double; P25's copy of the column
  test with a reach a hertz longer and of the filler with one magnitude moved
  (`test_p25_rendering_copies`); the tone model's test by structure alone, without the row, without
  the pitch, without L, without K, without the zero magnitudes, without the positive tone magnitudes,
  without the levels' ratio, with a ratio tolerance of 1e-5, with the ratio taken without the droops,
  without the finite `log2Ml`, with a finite test that lets infinities through, without the mark and
  without the flag. In the first fix round's list (46), five survived the tests as first written (the
  train in ProVoice's, D-STAR's and rate 33's renderings, across an octave, and the tone magnitudes'
  positivity; `build/mut/mutations2-first.log`); in this round's first run two more did (an octave
  bracket's upper edge taken inside, which an exact 3/5 with an even sum cannot show, and the tone
  magnitudes' positivity again, which the levels' ratio had come to cover) and one edit no longer
  applied to the code (the unfloored mean, then re-aimed; `build/mut/mutations3-edb277f.log`), and the
  tests were strengthened before the run above.
  The differential harness (outside the repository, `build/harness/pd_diff2.c`, built against the
  baseline's library and the head's in the ten local configurations) enumerates rather than draws:
  the pitch every D-STAR code (0 to 127) and every rate-33 code at every V/UV row (0 to 31) decodes
  to; caller states forged on every pitch of a list of 1,004 (D-STAR's 124 decoded pitches, P25's 208
  in three forms, the old law's, the grid's and eq 46's nominal, rate 33's 120 at row 16 and 120 on
  the band grid at row 18, and 16 extreme float classes: zeros, subnormals, `FLT_MIN`, `FLT_MAX`,
  infinities, NaNs, a negative, pi, 2^31, 2^70, 2^-32, 1) crossed with each of 12 voicing classes
  (none, all, a low region, alternate bands, a lone harmonic, the top only, the isolated last band,
  states 2, -1 and 7 among voiced and unvoiced bands, a tone structure, the ends); on each such state
  one of 24 mark-or-word cases (every mark, values that are none, pulsed-row words, tone words real
  and forged with the tone structure, the review's word, other words in `Vl[0]`), one of 12
  log-magnitude classes in `log2Ml` (NaN, infinities, -20.1, -40, 1e30, `FLT_MAX`, subnormal, zeros),
  one of 7 magnitude classes in one `Ml` (none, NaN, +infinity, 1e30, `FLT_MAX`, subnormal, -0), one
  of 9 slot cases in slots 34, 35, 37 to 39 and 95 and one validity count in slot 36 (-11 to 3), each
  in its own field; a forged state on every frame of every stream, put into one of the three slots
  before one frame of each of 23 paths (the D-STAR, P25 and rate-33 data paths, float and short;
  D-STAR's, rate 33's and P25's frame paths, hard and soft, float and short; ProVoice's four frame
  paths and the P25 data path with its flag), the frames enumerating each path's codes (rate 33's
  every b0 at every row, ordinary voice frames of rows 18 to 31 included) with their other bits from
  a fixed sequence and the error counts and flags enumerated; the decode-only frame entry points,
  D-STAR DV packing, `mbe_classifyAmbe2450Frame()`, the tone, silence and comfort-noise synthesizers,
  `mbe_spectralAmpEnhance()`, `mbe_applyAdaptiveSmoothing()`, `mbe_useLastMbeParms()`,
  `mbe_synthesizeSpeechf()` and `mbe_synthesizeSpeech()`, the three encoders through their float,
  short and frame entry points over forged histories, and TETRA (`mbe_tetraDecodeFrame()` over 6,000
  frames with bad-frame patterns and a reset); every byte of all three returned states hashed after
  every frame, no slot zeroed. It prints its own coverage, by factor assignment: `COVERAGE forged
  states 12048 (every frame of every stream), pitches 1004, pairs reached: vm 288/288 vl 144/144 vs
  108/108 ml 288/288 ms 216/216 ls 108/108 mmc 168/168; pitch x voicing 12048/12048; 1504 states take
  a tone structure over their voicing class` (v voicing, m mark, l log class, s slot case, mc
  magnitude class); and as the forged states hold them, after mark cases 19 to 21 write the tone
  structure (class 10's bands, harmonics 4 and 7 at 120) over the voicing class: `COVERAGE as held
  (after the tone structure's overwrite): avm 255/288 avl 144/144 avs 108/108; distinct (pitch
  index, Vl[1..56]) 10176 of 12048 states; the 23 decoder streams forge every index on every frame;
  the TIA rendering (mbe_synthesizeSpeechf()) forges every second index, the encoders' histories
  every seventh of their 2,000 frames, the public helpers one index a frame`. So 33 mark-by-voicing
  pairs are not reached (the three tone-structure marks with the other 11 voicing classes), and
  10,176 distinct pairs of pitch and voicing are, not 12,048 (classes coincide at some harmonic
  counts). ProVoice entered from the states the D-STAR and rate-33 streams leave, with the low cut's
  slots changed in all three, plays the same samples for 200 frames, and every byte of the returned
  states is compared after every frame: `prev_mp.noiseOverlap[38]` and `[39]` differ after the first
  frame (they keep the values they entered with), and from the second frame on the states are
  equal; the baseline (`f725a66`) and the branch below (`1e531f8`) do the same, so
  it is the branch point's behaviour, not this change set's. Built against the branch below's head
  `1e531f8` and the merged head `ecb78d0`, again against `5369c52` (the tone predicate's structure
  loop moved into a helper), each time with every PROTECTED line (4,252, the two coverage lines and
  the ProVoice lines among them) the same in each of the ten local configurations; and, in its
  present form (three coverage lines), against `c18d306`, the last commit that changes the library
  or a test (the third review's fixes in), in the ten configurations: every PROTECTED line (4,253)
  the same in each, and only P25's lines change (215: the 208 decoded pitches, its six streams and
  the IMBE encoder's `cur_mp`), P25's lines in the dev-debug build those of the head before the
  merge. In the second fix round each library edit was first run against the commit before it, every
  line the same: the pulses' ownership (`e249b64` against `f24ed64`, the first form of the harness)
  in the 32-bit Release and GCC fast-math LTO builds, the analyzers' zeroed arrays (`edb277f`
  against `e249b64`) in those and the Clang LTO build. What it leaves out: MSVC
  and macOS builds, concurrent threads, caller states outside these classes, and combinations beyond
  the pairs crossed (no triple of factors is enumerated in full).
  Valgrind (memcheck, full leak check, the merged head's dev-debug build) on `test_params`,
  `test_p25_tones`, `test_frame_paths`, `test_subframe_synthesis`, `test_p25_rendering_copies`,
  `test_golden_pcm` and `test_ambe_block_lengths`: no invalid access, no uninitialised read, nothing
  definitely or indirectly lost; the five tests that synthesize through the library report two
  blocks "possibly lost" (1,095 and 4,167 bytes) and two still reachable (144 and 2,752 bytes), all
  the thread-local transform plan's (`mbe_fft_plan_alloc()`, which a thread keeps until it exits),
  which valgrind counts as two errors; a Debug build of the branch below (`1e531f8`) shows the same
  four blocks in `test_params` and `test_subframe_synthesis`. On `c18d306` (`test_p25_tones` again on
  `e42a5a1`) the tests this round
  changed (`test_p25_tones`, `test_p25_rendering_copies`, `test_subframe_synthesis`) show the same:
  no invalid access or uninitialised read, the same two possibly-lost blocks in the two that
  synthesize through the library, none in the copies test.
  Preflights: `tools/preflight_ci.sh` and `tools/quality_preflight.sh`, both with `--base 1e531f8`
  (the branch below as merged), pass on `5369c52` (their first run, on `0e08319`, failed Lizard's
  complexity limit on the tone predicate, 16, and IWYU on the copies test's include of `mbe_math.h`;
  `5369c52` fixes both); and again on `e42a5a1`, after the third review's fixes (their first run
  there, on `c18d306`, failed clang-tidy and scan-build on the tone test's floating-point loop
  counter, which `e42a5a1` replaces with an integer rung). The quality preflight's Gitleaks scan again covered no commit (`0 commits
  scanned`, `scanned ~0 bytes`: the container cannot resolve a linked worktree's Git directory, and
  the job reports success). Instead the pinned image (`tools/ci-dependency-pins.env`) with the
  repository's `.gitleaks.toml` scanned `f725a66..HEAD` in a clone of the branch outside the
  worktree (`build/scripts/secretscan.sh`, its log `build/secretscan.log`, both outside the
  repository): at `0e08319` the range held 45 commits by `git rev-list --count`, two of them merges
  (`ecb78d0`, and `b88c102` of the branch below), and the scanner reported `43 commits scanned`, no
  leaks; the commits after it are scanned the same way before the branch is handed over.
- AMBE 3600x2450 frame-type handling per TIA-102.BABA-1 and as DVSI's decoder plays it: frame classification, silence history freeze, silence frames on 15 bands at their level with the previous model held to mid-frame, repeats, erasures (tone ID 255 among them, on the hard, soft and parameter-bit paths, and the all-ones frame of DVSI's DTX vectors, a tone frame with the unusable index 180), tones (including unverified tone frames; the components a tone ID plays, in 2^-19 cycle per sample, against values written out for DTMF, KNOX, call-progress and single tones, within 0.1 Hz of the nominal single tones and 20 Hz of the nominal dual ones, Table 9's nominal lookup unchanged; a tone after a voice frame starting as a fresh pulse train against an independent rendering in double, its second frame and a tone changed to another going on from the phases), muting (the mute noise level, the last frame fading out against an independent rendering, the model kept, without a row's mark after erasures, and the noise generator's single frame through every kind of mute), V/UV rows 13, 15 and 17 (each row's own voicing and its mid-frame voicing after a voiced and after an unvoiced frame, row 10 as the control; row 15 after a voiced frame leaving the pulse phase at the phase before plus 80 samples at each frame's pitch, row 17 then row 15 starting the train at mid-frame, and row 17 then row 16 launching no pulse; the pulse train running on through a row-17 frame (`test_train_through_row17`, every mid-frame voicing and pitch stated by the test itself, not read from the library's mid-frame model: after a voiced frame the row-17 frame leaving the phase before plus 160 samples at that frame's pitch; row 13, row 5 at another pitch and another row-17 frame after it launching their first pulse (1 - p) / g samples in, p the phase the row-17 frame stored and g the mid-frame pitch, and leaving p + 80 g + 40 (g + f') cycles (80 g for the row-17 frame's unvoiced end in place of the last term), each case first shown to differ from a fresh start by more than 0.01 cycle, and the pulses the row-17 frame left sounding in them, the same stream with those pulses dropped and their phase kept playing other samples before its first new pulse; row 15 after it starting a new train at mid-frame, its first pulse at sample 80, the dropped pulses changing nothing; after a tone the row-17 frame from silence launching a new train at its first sample, which the frame after it takes on; under an error-rate mute, at one pitch code, the unheard row-17 frame storing the phase before plus 160 samples at that pitch whatever its end voicing, and the frame that ends the mute starting from silence, its first pulse at its first sample (row 15: at sample 80); advancing, started from fields no frame leaves, leaving the same pulse state, `PSIl`, `PHIl` and mark as rendering, compared as the file's fast-math rule asks); a lost frame after each of them repeating the model with its own voicing through and leaving no mark, and a tone frame right after each playing as after the same frame without its mark, the speech it fades out repeating the model's own voicing; each as the first frame after a mute and after a tone through `mbe_processAmbe2450Dataf()`, its first pulse at the frame's start for rows 13 and 17 and at sample 80 for row 15, which pins what the library does; P25's onset mark in a caller's state playing as none before an erasure or a tone), the frame after a mute or a tone rendered from silence from the kept model (plus a tone's tail), and lost frames repeating a tone bit for bit until the repeat limit fades it out; a held-tone word in the caller's state that names an ID with no tone playing, frame after frame and state after state, as the word without the held tone; and the played tone table pinned (every AMBE+2 ID's and D-STAR index's components, by an FNV-1a hash that every one-unit change of a dual-tone row's fundamental or harmonics was checked to move) (`test_ambe2450_frame_types`)
- D-STAR voicing: every b1 value at four pitches, each voiced band reaching 250 Hz into the band below (`test_ambe2400_frame_types`)
- D-STAR frame classes outside plain voice: b0 125, tone frames without a tone and b0 127 frames repeating like error repeats through the hard, soft and parameter-bit APIs and sharing the repeat limit, the null frame and b0 120..124 as voice (b0 124 with the AMBE+2 silence model, all 16 b1 values playing alike on every API and committed to the prediction history), zero-volume tones, C0 words with four errors (every pattern of weight 4 and none lighter, hard and soft) playing as lost frames whose decode the next frame predicts from, through a lost frame or a tone after them, runs of such frames playing sample for sample as lost frames (the repeated model rendered with its own `log2Ml`: the last synthesized frame holding the repeated model's `log2Ml` and the prediction history the held one's, and the frame's phases those of the first, not the held history's nor those of the `Ml` it renders), a lost frame or a tone right after such a frame playing sample for sample as after a lost frame in its place (the speech under the tone's first frame included), the same with one state as `cur_mp` and `prev_mp_enhanced` as with three (`test_one_state_held_history`: sequences of voice, lost and uncorrectable-C0 frames, single and dual tones and the repeat-limit mute, none with channel errors, each playing sample for sample as its copy with lost frames in place of the uncorrectable ones up to the first voice frame after an uncorrectable one, every repeat after a repeat-limit frame up to the next voice or tone frame at the mute noise's level, and the prediction history, `log2Ml`, L, w0 and the voicing, the same after every frame in both arrangements; `test_repeat_limit_mute_holds`: the repeats after the repeat limit's mute at the mute noise with one state as with three), corrections repeating no frame and tones with corrections playing, lost frames after a single tone (the tail fading out, then silence) and after a dual tone (the tone going on until the repeat limit), tones muted with the error level, the repeat-limit mute (fade-out, DVSI-level noise, prediction kept), the voice after it starting from the prediction history (against a reference that takes the history's `log2Ml`, the stream's silent last frame given other `log2Ml` first; that frame alone changing; a held or initial history fading in from silence; a b0 124 or still-muted frame after it) and the error-level mute and its noise level, and the all-ones DV frame as a lost frame for runs of 1, 3 and 20, the 80-sample crossfades into and out of a tone against a stream with a lost frame in its place, a tone's volume step ramping over 80 samples before the clip while a new index starts at once, and D-STAR tones at the frequencies of the AMBE+2 tone of the same nominal pair (the DTMF row tone first, call progress 144..147 as 160..163, single tones as AMBE+2's), started as a fresh pulse train after a frame that was not a tone frame, a held-tone word in the caller's state that names an index with no tone playing as no held tone (the next tone frame, voice frame or lost frame and the frame after it bit for bit as without it), and the other codecs' frame marks in a caller's history or last frame played (P25's onset columns, rate 33's rows 13 and 17) playing as none before a lost frame and a tone's first frame and in the voice frame after, the late onset acting only after a frame with no voiced band, as before those marks (`test_foreign_marks_in_history`) (`test_ambe2400_frame_types`)
- float-to-int16 conversion parity checks, and the D-STAR and rate-33 paths' round-to-nearest conversion
  in every compiled variant against `roundf()` (exact halves and the floats beside them, the limit,
  infinities, NaN, subnormals and a sweep over float bit patterns, in every lane and at unaligned
  addresses), with `mbe_floattoshort()` still truncating (`test_pcm_rounding`); the AMBE `short` entry
  points rounding their float twins' output, the IMBE ones truncating, and both AMBE paths' muted
  16-bit noise; the pitch the P25 frame wrappers store, hard and soft (within an ulp of DVSI's
  pitch and no further from eq 46's than the law's largest offset, ProVoice's within an ulp of
  eq 46's), and a repeat playing the repeated model's pitch, ProVoice's nominal one and P25's
  DVSI's (`test_imbe_repeat_keeps_its_pitch`); and D-STAR tone streams against the fresh
  start computed in double from its formula (`cos(2 pi F (t - 16) + p)`) past the first frame's
  crossfade, within 3e-4 of the tone's amplitude, and against the public
  `mbe_synthesizeTonefdstar()` continuing the phases the first frame set, bit for bit in an IEEE
  build and within 1e-5 * (1 + peak) under fast-math, peak the largest magnitude of the reference
  frame; in a build without tones, silence past the first frame's crossfade, by the samples' bits
  (`test_frame_paths`)
- golden PCM hash regression checks, by kernel ISA: `test_golden_pcm` runs the CPU probe's choice,
  `test_golden_pcm_sse2` and `test_golden_pcm_avx2` force either copy of a run-time AVX2 dispatch build
  and skip (exit 77) where the build or CPU cannot run it; `MBELIB_TEST_REQUIRE_AVX2=ON` makes the AVX2
  skip a failure
- SIMD architecture detection checks, including where the AVX2 dispatch applies
- run-time AVX2 dispatch (`test_cpu_dispatch`): the CPU probe against the compiler's CPU-feature
  runtime, the ISA selection rules, voiced oscillators that differ between SSE2 and AVX2 only within
  tolerance, tone fits within 1e-9, and soft ECC decoders and the TETRA adaptive codebook bit-identical
  across ISAs; `test_simd_numerics`, `test_tetra_simd` and `test_tetra_decoder` repeat their oracles and
  stream hashes under every ISA the build and CPU can select
- the AVX2 translation units define only `*_avx2` externals, so no AVX2 code can be linked into a
  baseline caller (`test_avx2_symbols`, Linux, non-LTO configurations)
- the host-tuning rules in CMake script mode, then configure-only checks of per-configuration
  settings (`test_build_configuration`; the configure checks need Ninja and skip when the forwarded
  compiler command selects a target ISA): every skip rule of `cmake/MbeNativeArch.cmake` (including
  a target ISA in the compiler command, quoted or tab-separated flags and MSVC's `CL`/`_CL_`) and the
  MinGW/Cygwin assembler requirement; Ninja Multi-Config, with or without `CMAKE_BUILD_TYPE` and
  from a source path with spaces, gets SIMD in both configurations (64-bit targets) and LTO only in
  Release; Release gets the native flag (when the host allows it) and Debug never does; a user
  `-march` or `MBELIB_ENABLE_NATIVE=OFF` turns it off; and a portable Release keeps its AVX2 copies
  when only Debug targets AVX2; one build directory reconfigured for `-m32` or (Clang) AArch64 probes
  the architecture again; a universal `CMAKE_OSX_ARCHITECTURES` off Apple keeps the dispatch; and the
  parent's `CMAKE_OSX_ARCHITECTURES` is forwarded. Exported build settings (`CFLAGS`, `CL`, `CMAKE_GENERATOR`, ...)
  are removed from the nested configures' environment
- quality-tool filesystem helper checks (`test_quality_fs`, built with `MBELIB_BUILD_TOOLS=ON`)
- quality-tool CLI alignment, accepted/rejected frame widths, and private output
  permissions (`test_quality_tools`, with tools enabled and Python 3 available)
- DVSI test-vector import: deinterleave bijection, TIA-102.BABA 7.5 code spacing,
  fixture round trips and preserved channel bits (`test_quality_dvsi`)
- TETRA decoder: bit unpacking (tables 3 and 4), every pitch index, pulse positions and absent
  pulses, sharpening, the 16-bit wrap of an overflowing codevector, LSP lookup, boundary widening
  and the rejection of crossed vectors,
  interpolation, LSP to A(z) against a double-precision expansion and against coefficients read
  back from ETSI's reference decoder, the synthesis filter, the
  adaptive codebook, log2/pow2, the gains against values ETSI's reference decoder produced on the
  test's own streams (each ruling out another plausible computation) and at the extremes
  (including a code gain past 16 bits), the floor of the code-vector energy, the 27/25 energy
  limits, concealment state (energy decay and floor, held LSPs, lag and indices, silence and the
  reference decoder's LSP state before the first good frame, ended by any good frame), the
  API contract (NULL arguments, rejected frames leaving state and output untouched, ignored
  bad-frame payloads, reset, interleaved contexts), and stream hashes asserted on every build
  (`test_tetra_decoder`); the saturating fixed-point operators at their edges
  (`test_fixed_point`); independent ordered arithmetic oracles for SIMD and scalar
  synthesis, fractional pitch, pulse shaping, excitation, energy and LSP interpolation
  (`test_tetra_simd`), including saturation, wraps, short-lag dependencies, unaligned
  buffers and tails; the TETRA entry points through the shared library
  (`test_tetra_shared_api`); TETRA argument and bit validation (`test_input_validation`); the
  table generators and the committed tables' shape and provenance (`test_tetra_tables`); the
  TETRA comparison tool, and with tools enabled `mbe_tetra_decode` (`test_tetra_conformance_tools`,
  `test_tetra_decode_tool`). `fuzz_tetra_decode` fuzzes the decoder through the public API.

PCM conversion is always compiled with IEEE semantics, including fast-math builds, so NaN/Inf handling is preserved.

Test assertions remain enabled in every configuration. MSVC force-includes
`tests/enable_assertions.h` to undefine `NDEBUG` before reading test sources,
avoiding the D9025 warning from combining `/DNDEBUG` and `/UNDEBUG`.
`test_assertions` checks that the assertion expression actually executes.

Run the default test suite with:

```sh
cmake --preset dev-debug
cmake --build --preset dev-debug -j
ctest --preset dev-debug --output-on-failure
```

## Continuous Integration

GitHub Actions runs tests and quality checks on pull requests and pushes to the
primary branch. Required checks include cross-platform builds, sanitizer tests,
static analysis, CodeQL (C/C++, workflows and Python), PR fuzzing, repository
security guardrails, workflow linting, dependency review, secret scanning, OSV
scanning, and install/consume checks.

Fast-math changes what the optimizer may assume, so the shipped configuration
gets jobs of its own. On Linux the whole suite, quality tools included, also
runs in `dev-release` builds (fast-math, LTO, SIMD) with x86-64-v3 code
generation, once with GCC and once with a current Clang from the Arch
toolchain. Clang folds bit tests on floats it computed (an exponent-mask
NaN/Inf test once passed overflowed encoder models on macOS this way) without
a diagnostic, so only a test that feeds the guard a NaN or infinity catches
that. With `-Werror`, its `-Wnan-infinity-disabled` does turn `isnan()`,
`isinf()` and NaN or infinity constants in fast-math code into a build error;
the job first probes that the compiler still enforces it. A third job builds and links everything for x86-64-v4 with
GCC and LTO and warnings as errors, since GCC reports some diagnostics only for
AVX-512 targets at LTO link time; the runner may lack AVX-512, so that build is
never run. To test for NaN or infinity in fast-math code, compare the
magnitude of the bit pattern against a bound, or compile the translation unit
with IEEE semantics as `mbe_pcm_convert.c` is.

Consumers do not build a preset. dsd-neo's CI, distribution packages and the
vcpkg port configure a plain `CMAKE_BUILD_TYPE=Release`, whose defaults are
SIMD on, fast-math off, warnings as errors and host tuning
(`MBELIB_ENABLE_NATIVE`) on; packagers that ship binaries pass
`-DMBELIB_ENABLE_NATIVE=OFF`. The IEEE-only SIMD batches
(`!defined(__FAST_MATH__)`, such as the AMBE codebook search) compile only
there, so a fourth Linux job builds and tests exactly that configuration with
the runner's GCC and fails if those defaults change. v2.5.0's codebook search
failed to build under it (GCC 13 `-Waggressive-loop-optimizations` on an `int`
row offset) while every preset job passed. The same job also builds and tests
the portable variant (`-DMBELIB_ENABLE_NATIVE=OFF`), which keeps baseline x86-64
code generation and dispatches the AVX2 kernels at run time.

It then configures Ninja Multi-Config without a build type, with LTO, and runs
the suite for Release and Debug, so a multi-config Release keeps SIMD, LTO and
the host tuning.

Release assets and the install-and-consume jobs build `dev-release` with
`-DMBELIB_ENABLE_NATIVE=OFF`, so the shipped libraries need only x86-64 v1
(SSE2) or the ARM64 baseline and use AVX2 through the run-time dispatch; the
install-and-consume jobs also run that build's suite (on Linux and Windows
with `MBELIB_TEST_REQUIRE_AVX2=ON`). The `dev-release` legs of the build matrix
keep the native default and also require the AVX2 golden test on Linux and
Windows. A Debug job builds a portable SIMD tree and checks the strict golden
hashes for both dispatch paths (`test_golden_pcm_sse2`, `test_golden_pcm_avx2`
with `MBELIB_TEST_REQUIRE_AVX2=ON`), then builds a Debug `-march=x86-64-v3`
tree and requires its `gen_golden` output to equal the dispatched
`gen_golden avx2`.

Every required check runs on every pull request. A required check that a path
filter keeps from starting never reports, so the pull request waits on it
forever. A job that only matters for some paths therefore runs on every pull
request and skips its work when `tools/ci_changed_files.sh` finds nothing for
it, as PR fuzzing does. A job that `if:` skips still reports success under its
name, so a required job that only makes sense on a pull request, such as
dependency review, lives in `guardrails-pr`, a workflow that only pull requests
start, with no event condition.

## Regression Test Requirement

At least 50% of bugs fixed in the last six months should include regression
tests. A pull request that fixes a bug should add a regression test unless:

- the behavior cannot be reproduced reliably in automation
- the fix is entirely documentation or packaging metadata
- a better guardrail exists, such as a static-analysis rule or workflow check

When no regression test is added for a bug fix, the pull request must explain
why.

## Major Functionality Test Requirement

Major new functionality must add or update automated tests. Major functionality
includes:

- public API changes
- codec or DSP behavior changes
- external input handling changes
- dependency changes that affect compiled code
- security-sensitive workflow or release changes
- installation and packaging behavior changes

## Coverage Target

The project target is at least 80% statement coverage for project-owned source
files when measured with a FLOSS C coverage tool. Vendored code under
`src/external/` is excluded from project-owned coverage accounting.

Coverage evidence should be generated from an instrumented Debug build and
recorded in a pull request, issue, or CI artifact before claiming coverage-based
badge criteria.

## Dynamic Analysis

For memory-safety-sensitive C changes, run sanitizer tests:

```sh
cmake --preset asan-ubsan-debug
cmake --build --preset asan-ubsan-debug -j
ctest --preset asan-ubsan-debug --output-on-failure
```

Frame-processing paths are also covered by ClusterFuzzLite PR fuzzing with
AddressSanitizer and UndefinedBehaviorSanitizer, including fixed-size hard/soft
frame decode paths and parameter synthesis paths. The fuzzers run on pull
requests that change the library, its headers, the fuzz harnesses or the fuzz
build (`fuzz_targets` in `tools/ci_changed_files.sh`).

## SIMD and replay benchmarks

Build the opt-in benchmarks with the library configuration being measured:

```sh
cmake --preset dev-release -DMBELIB_BUILD_BENCHMARKS=ON
cmake --build --preset dev-release -j
build/dev-release/bench_synth
build/dev-release/bench_synth steady mixed
build/dev-release/bench_synth gradual voiced
build/dev-release/bench_synth jump voiced
build/dev-release/bench_replay encode ambe2400 3 7 < speech.raw
build/dev-release/bench_replay hard ambe2450 3 7 < ambe2450.frames
build/dev-release/bench_replay soft imbe7200 3 7 mixed < imbe7200.frames
build/dev-release/bench_tetra 3 7 < /path/to/TE1_46.CHD
```

`bench_synth [legacy|steady|gradual|jump] [mixed|voiced]` defaults to legacy
mixed, the original workload, unchanged: w0 jumps between 0.09 and 0.11 at a
fixed L = 40 (harmonics above Nyquist are skipped by the synthesizer), with no
warm-up and state carried across runs. The quality A/B runner times
`bench_synth` without arguments against older builds, so that default must stay
the same work. In the other workloads every model keeps its harmonics below
Nyquist: steady and gradual pitch exercise interpolated low harmonics, jumps
exercise windowed synthesis, and each of ten measured runs starts with the same
state after one untimed warm-up. Every workload prints the `avg:` line the A/B
runner reads.

`bench_replay encode|hard|soft ambe2400|ambe2450|imbe7200|imbe7100
[repeats=1] [runs=7] [uniform|mixed]` loads standard input into memory, then
allocates the encoder and warms up before timing. Encoding takes mono 8 kHz s16le PCM, zero pads
a partial last frame and rejects partial samples; ProVoice has no encoder.
Decoding takes one rectangular, row-major bit row per line: 96 bits for AMBE,
184 for P25 or 168 for ProVoice, as written by `mbe_quality_reframe` (parameter
rows are not accepted). Soft replay uses confidence 255 or a reproducible mix
of all 256 confidence values, including zero. Use channel-error fixtures to
measure noisy decoding. Inputs are capped at 65,536 frames.

The replay benchmark reports median CPU microseconds per 20 ms frame, excludes
loading, allocation and stream resets, and checks a checksum against its warm-up
on every run. Repeats concatenate the input within a run; each run resets the
encoder, parameter triplet and RNG. Checksumming is included in the timing.
Pin runs to one CPU when possible and alternate matched baseline/candidate
builds. Compare identical compiler, fast-math, LTO and `MBELIB_ENABLE_NATIVE`
settings: a default Release build is tuned for the build host, so pass
`-DMBELIB_ENABLE_NATIVE=OFF` to measure the portable build. `-march=x86-64-v3`
is a separate compile target; portable SIMD builds still require only SSE2 on
x86 (running the dispatched kernels as AVX2 on AVX2 CPUs) and use NEON only
when enabled by the existing target selection. Native ARM
measurements are required before claiming NEON performance improvements.

`bench_tetra [repeats=10] [runs=7]` applies the same in-memory replay method to
TETRA's 30 ms frames. Standard input is CHD format: 138 little-endian 16-bit
words per frame (bad-frame indicator followed by 137 bits). It rejects empty,
truncated or invalid good-frame input, caps input at 65,536 frames, and ignores
bad-frame payloads as the public decoder does. The reported median excludes
loading, allocation, resets and one warm-up run. Repeats concatenate the input
without resetting; every measured run starts from reset and must reproduce the
warm-up checksum. One PCM sample per frame contributes to the timed checksum;
use `mbe_tetra_decode` and byte comparisons to check complete PCM identity.
Include speech, losses and overloaded synthetic streams in comparisons. See
[TETRA SIMD measurements](tetra-acelp.md#simd-optimization) for the measured
speed and bit-exactness checks.

The soft ECC codebooks are immutable 32-bit codewords (32 KiB in all) in the
original enumeration order, with the codeword bit of each data bit. Regenerate
them with Python 3.10+ using `python3 tools/gen_ecc_codebooks.py >
src/internal/mbe_ecc_codebooks.h`. A decode starts from the codeword carrying
the hard decision's data bits. A received word that is that codeword is
returned at once: it costs nothing, and of the codewords that cost nothing it
ranks first. Any other word is searched: the search builds small cost tables
on the stack (data bits by enumeration index, parity bits by value), starts
from that codeword's cost and skips only codewords that cost more, so the
result is exactly that of a full scan, ties included. SIMD builds screen
candidates in groups of four (SSE2) or eight (AVX2), then resolve surviving
candidates in enumeration order with the same exact integer scores and tie
rules. AVX2 gathers pairs of 16-bit parity costs without reading beyond the
cost table. Other targets use the scalar search. The search has no lazy
initialization or shared mutable state; a portable x86-64 build only reads the
run-time dispatch's cached ISA choice, which every thread computes identically.
`test_ecc` checks every clean codeword, at full confidence and under six
confidence patterns (zero confidence among them, where other codewords can
cost as little), then compares 1024 pseudo-random inputs per code against an
independent oracle that scores every codeword, with zero, full, binary,
single-bit, small tied, mostly confident and random confidences. It compares
codewords with one to three different bits flipped at zero confidence, which
are no codewords and so go through the search, against the same oracle, and
checks a case in which the soft decision must override the hard decision
although the received data bits agree with it.

`test_simd_numerics` compares tone fits to direct trigonometric projection and
an independent solver (relative tolerance 1e-9), and, where the target has
vector voiced kernels, voiced output to the scalar formulas (absolute
tolerances 3e-5 for unit-scale windowed components and 2e-4 for
interpolation). It exercises unaligned inputs and scalar tails. Targets
without vector kernels compile the historical voiced loops in `mbelib.c`
unchanged, so their PCM stays bit-identical, under x87 excess precision and
fast-math too; SIMD recurrence batching permits bounded PCM rounding
differences and has separate Debug goldens. Each SIMD interpolation lane
seeds from the phase's own cosine and sine, rotated by its offset in double,
so caller-supplied phases far beyond the decoder's range, up to near the
largest float, follow the exact phase. `mbe_voiced.c` is compiled with IEEE
semantics even in fast-math builds, as `mbe_pcm_convert.c` is, since a
fast-math math library need not reduce huge cos/sin arguments exactly. Phase and RNG updates stay outside
these kernels. Tone refinement
trials remain sequential and double precision; ARM32 uses scalar tone fitting.
Without fast-math, the pitch low-pass filter assigns independent output samples
to SIMD lanes; its 21 taps still accumulate in order. Fast-math builds retain
the compiler-vectorized filter/energy loop. IEEE AMBE codebook searches assign
complete candidate distances to lanes and resolve equal minima by original row
index; fast-math builds retain the compiler-vectorized codebook loop.
Neither kernel reduces partial floating-point sums across lanes. Remaining
analysis reductions retain their existing evaluation order.

`mbe_floattoshort` classifies NaNs in integer SIMD lanes before arithmetic,
converts them to zero, and clips infinities to the signed limit. This avoids a
separate scalar scan and keeps exceptional frames vectorized. Its tests include
arbitrary IEEE bit patterns, signaling NaNs, unaligned inputs/outputs and frame
boundary canaries. Hard/soft bit validation uses compiler-vectorized integer OR
reductions in SIMD builds, with tests for every invalid byte value and tail
length. All three frame encoders use the same hard-bit validator. Reliability
bytes are excluded from validation.

Unvoiced noise generation advances four independent integer LCG lanes (eight
on AVX2), preserving the original sample order and stream state exactly.
`test_simd_numerics` compares every one of the 53,125 states with a scalar
recurrence across three consecutive frames, including the overlap and next
seed. This oracle also runs on ARM32 and AArch64 under QEMU.

These kernels now serve the public `mbe_synthesizeSpeechf()`. The codec paths
render in subframes ([DVSI parity](dvsi-parity.md#synthesis-in-subframes)) with
scalar loops: each pitch pulse's harmonics by rotation, and DVSI's 16-bit noise
generator, whose states a frame draws in eight lanes that each step eight
states at once in 16-bit arithmetic (the same states in the same order:
`test_unvoiced_subframes` compares a frame drawn from every state with the
recurrence). They share the transform, the phase regeneration and the PCM
conversion with the specification's rendering.

On AVX2 compile targets, float kernels process eight lanes, including voiced
oscillators, unvoiced windowing/overlap and PCM conversion. Portable x86-64
builds still require only SSE2; their voiced oscillators, tone fits, soft ECC
search and TETRA adaptive codebook also have an AVX2/FMA copy that runs on
CPUs with the x86-64-v3 set (see
[run-time dispatch](build-installation.md#run-time-avx2fma-dispatch)); the
other kernels keep SSE2 widths. Debug AVX2 has its own float and AMBE sequence
PCM goldens because wider oscillator recurrences change rounding, and the
dispatched AVX2 copy reproduces a Debug `-march=x86-64-v3` build exactly.
Regenerate them in a Debug SIMD build with `gen_golden sse2` and
`gen_golden avx2` (`auto`, the default, prints the probe's choice), and check
`gen_golden avx2` against `gen_golden` from a Debug build with
`-DCMAKE_C_FLAGS=-march=x86-64-v3`.

`tools/bench_compare.sh` runs the encoder, PCM conversion, unvoiced synthesis
and all six steady/gradual/jump mixed/voiced workloads in matching Release
builds with explicit SIMD disabled, enabled (both portable, so the SIMD build
dispatches AVX2 on AVX2 CPUs) and enabled with host tuning. Compilers can still automatically
vectorize C loops when explicit SIMD is disabled. Use `bench_replay` in those
same trees for real speech, all four decoders, tone frames and channel errors;
synthetic synthesis timings alone do not establish end-to-end gains. It also
builds `bench_tetra`; supply local CHD inputs to compare TETRA in the same trees.

The October 2026 SIMD audit compared against `a962d3d` on a Ryzen 7 5800HS,
with GCC 16.2.1, Release, fast-math and LTO enabled on both sides. These are
medians of three alternating comparisons pinned to one CPU; lower is better.
Replay uses the DVSI `dam` speech input. SSE2 and AVX2 are separate compile
targets, each compared with the same target before the changes.

| Workload | Target | Before | After | Time reduction |
| --- | --- | ---: | ---: | ---: |
| PCM conversion, 10 million frames | SSE2 | 0.6210 s | 0.2845 s | 54% |
| PCM conversion, 10 million frames | AVX2 | 0.2800 s | 0.2033 s | 27% |
| Mixed synthesis, steady pitch, 2,000 frames | AVX2 | 9.751 ms | 8.531 ms | 13% |
| Mixed synthesis, pitch jumps, 2,000 frames | AVX2 | 9.612 ms | 7.761 ms | 19% |
| D-STAR hard-frame replay | AVX2 | 3.726 us/frame | 3.466 us/frame | 7% |
| AMBE+2 hard-frame replay | AVX2 | 4.025 us/frame | 3.768 us/frame | 6% |
| IMBE hard-frame replay | AVX2 | 4.101 us/frame | 3.889 us/frame | 5% |

Explicit codebook batching reduced the complete AMBE residual quantizer's time
by 9–17% with GCC and IEEE math. Its fast-math version measured slower than the
compiler-vectorized loop, so explicit batching is enabled only with IEEE math.
Fast-math encoder timings stayed near their previous values. With Clang 23.1.1,
SSE2 and neither fast-math nor LTO, the pitch-filter changes reduced whole-frame
encoding time by about 7% for the AMBE codecs and 10% for IMBE. FFTs, tone
fitting and analysis reductions already had SIMD or effective compiler
vectorization; changes there require their own measurements. ARM numerical
checks run under QEMU establish correctness, not native ARM speedups.

For the same audit, all three encoders produced identical parameter-bit rows
on all 79 cached DVSI inputs (27 speech and 52 tone vectors): 256,707 frames
per build pair. This held for GCC SSE2/AVX2 with fast-math and LTO, GCC SSE2
with IEEE math, Clang AVX2/IEEE SSE2 configurations, and Android NDK Clang 19
IEEE builds for ARM32 and AArch64 under QEMU. Compare each compiler and build
configuration against its own baseline; cross-compiler identity is not the
acceptance criterion.

Decoder comparisons used every cached DVSI frame stream, a fixed seed and the
baseline's alignment against both the original speech and DVSI's decoded PCM.
SSE2 output was byte-identical. Across GCC and Clang AVX2 builds, each changed
int16 sample differed by at most one unit, with no added clipping or nonfinite
samples. The largest absolute change in log-spectral distance was below
0.00018 dB against the input and 0.00007 dB against DVSI; envelope correlation
changed by less than 0.000008. Preserve these per-codec comparisons when
changing SIMD kernels: passing synthetic numerical tests alone does not
establish speech quality. Keep DVSI vectors and generated audio local.

The subsequent noise-generation and frame-validation changes were also
compared directly with the earlier PR head, `c204ce4`. Across GCC SSE2/AVX2
with fast-math and LTO, and Clang 23.1.1 SSE2 with IEEE math, all 79 DVSI inputs
produced identical encoder parameter bits, packed frames and decoded PCM:
770,121 compared frames and 123,105,600 PCM samples, with no nonfinite output.
The established DVSI quality scores therefore remain unchanged by these two
optimizations. Native ARM performance remains a separate measurement.

For this follow-up, five alternating comparisons pinned to one CPU on the
same Ryzen/GCC setup measured the following additional reductions relative
to `c204ce4`. The first three codecs replay the DVSI `dam` stream through
`bench_replay hard <codec> 10 7`; ProVoice uses a converted frame fixture.

| Hard-frame replay | Target | Before | After | Time reduction |
| --- | --- | ---: | ---: | ---: |
| D-STAR | SSE2 | 3.954 us/frame | 3.791 us/frame | 4.1% |
| AMBE+2 | SSE2 | 4.305 us/frame | 4.145 us/frame | 3.7% |
| IMBE | SSE2 | 4.394 us/frame | 4.238 us/frame | 3.6% |
| ProVoice | SSE2 | 4.721 us/frame | 4.568 us/frame | 3.2% |
| D-STAR | AVX2 | 3.493 us/frame | 3.140 us/frame | 10.1% |
| AMBE+2 | AVX2 | 3.802 us/frame | 3.456 us/frame | 9.1% |
| IMBE | AVX2 | 3.911 us/frame | 3.576 us/frame | 8.6% |
| ProVoice | AVX2 | 4.109 us/frame | 3.778 us/frame | 8.1% |

`bench_unvoiced 36 20000 7` took 7.2% less time on SSE2 and 20.8% less on
AVX2. A separate frame-packing microbenchmark cycled 256 deterministic input
bit arrays through two million calls to each `mbe_encode*Frame` API: AMBE
packing took 17–19% less time (about 14–15 ns/frame), and IMBE took 5–9% less
(about 12–25 ns/frame). These packing figures exclude PCM analysis and
parameter quantization; whole-frame decoding above establishes the gain
from batching the noise generator. With Clang 23.1.1 SSE2 and IEEE math,
the same replay comparisons took 6.6–7.9% less time.

### PR #130 synthesis and ECC optimization

The follow-up to the DVSI parity changes compares with PR #130 at `d2293ba`.
It batches the spectral-envelope phase filter across independent harmonics
using SSE2, AVX2 or NEON, keeping each harmonic's ten taps in order. Synthesis
computes the high-pass phase only for voiced oscillators, including harmonics
fading between different model lengths. On SIMD targets the gain loop is
separate from the phase calls, allowing compiler vectorization. The scalar
path retains its combined expression, including x87 rounding. The required
phases still use the double-precision `atan2` law.

AMBE+2's counted C1 path re-encodes the data from its first Golay decode to
count parity errors, removing the second decode and, for soft input, the
second codebook search. It saves the received word before packing output,
preserving callers that put the 49 output bits at the start of the input
frame. Encoder analysis, quantization, FFTs and oscillator kernels retain
the implementations from the earlier SIMD work.

These measurements use a Ryzen 9 5950X, GCC 16.2.1, Release with fast-math
and LTO on both sides. Each figure is the median of five alternating
baseline/candidate comparisons pinned to CPU 6. Each decode invocation uses
`bench_replay hard|soft <codec> 10 5` on DVSI's `dam` stream; ProVoice uses
the same speech encoded to IMBE parameter rows and reframed as ProVoice.
Input loading, allocation and warm-up are outside the timer. Lower is better.

| Complete frame decode | Target | Before | After | Time reduction |
| --- | --- | ---: | ---: | ---: |
| D-STAR, hard | SSE2 | 4.541 us | 3.739 us | 17.7% |
| AMBE+2, hard | SSE2 | 4.991 us | 4.130 us | 17.3% |
| IMBE, hard | SSE2 | 4.957 us | 4.155 us | 16.2% |
| ProVoice, hard | SSE2 | 5.314 us | 4.532 us | 14.7% |
| AMBE+2, soft | SSE2 | 7.502 us | 5.819 us | 22.4% |
| D-STAR, hard | AVX2 | 3.191 us | 2.964 us | 7.1% |
| AMBE+2, hard | AVX2 | 3.519 us | 3.283 us | 6.7% |
| IMBE, hard | AVX2 | 3.574 us | 3.349 us | 6.3% |
| ProVoice, hard | AVX2 | 3.761 us | 3.555 us | 5.5% |
| AMBE+2, soft | AVX2 | 5.948 us | 4.947 us | 16.8% |

The other soft decoders take 8.5–12.8% less time on SSE2 and 2.5–4.7% less
on AVX2. Scalar Release with IEEE math takes 8.1–9.3% less time for hard
decoding. Encoder medians on SSE2 and AVX2 stay within 0.2% of baseline.
With Clang 23.1.1, hard decoding takes 14.5–18.0% less time on SSE2/IEEE
and 17.3–21.5% less on AVX2 with fast-math and LTO. Longer Clang AVX2
encoder runs put AMBE medians within 0.5% of baseline, inside the observed
timing spread, and IMBE 1.8% faster.
Speedups depend on the compiler, target and speech; QEMU correctness checks
do not establish native ARM performance.

For each of six host build pairs, all three encoders produce identical
parameter-bit rows on all 97 cached DVSI inputs: 296,772 encoded frames.
The configurations are GCC scalar/IEEE, GCC SSE2/IEEE, GCC SSE2 and AVX2
with fast-math and LTO, Clang SSE2/IEEE, and Clang AVX2 with fast-math and
LTO. All 321 cached frame streams, including tones, synthetic inputs and
hard/soft channel-error vectors, contribute 55,073,120 compared PCM samples
per pair. Five configurations produce identical PCM. GCC SSE2 with
fast-math changes 298 samples, each by one int16 unit, as the compiler
vectorizes the gain calculation. There are no nonfinite samples, additional
clips, or changed tone, erasure, repeat, silence, mute or error counts.

For the changed SSE2 streams, scoring both outputs at the baseline's lag
changes log-spectral distance by at most 0.000004 dB against the input and
0.000002 dB against DVSI; envelope correlation changes by less than
0.0000004. Tone detector agreement, wrong-tone fraction and steady-tone
window counts are unchanged. The parity PR's existing tone-gate result
remains as recorded in [DVSI parity](dvsi-parity.md#gates). The converted
ProVoice fixture also produces identical PCM in all six host configurations.

ARM32 and AArch64 NEON each compare the same 321 streams under QEMU.
ARM32 PCM is identical. AArch64 changes 6,045 of the 55,073,120 samples,
each by one int16 unit. Its fixed-lag log-spectral distance changes by less
than 0.00004 dB against both the input and DVSI, and envelope correlation
by less than 0.000003. Tone detector agreement, wrong-tone fraction and
steady-tone window counts are unchanged; there are no additional PCM rail hits.

The phase oracle covers every supported model length (1–56), mirrored and
extrapolated envelope edges, zero magnitudes and both voicing classes. ECC
checks cover 512 arbitrary received words with hard and soft decisions,
including zero-reliability ties, and 32 in-place layouts. Scalar, SSE2 and
AVX2 Debug golden hashes are unchanged. Checks include the full Debug and
Release suites, scalar and AVX2 ASan/UBSan, NOTONES, x87 Debug, the six
Release configurations above, ARM32 and AArch64 NEON numerical/API tests
under QEMU with Android NDK Clang 18, and `tools/quality_preflight.sh`.
ARM32 NDK builds disable hardening because that compiler accepts the
stack-clash flag in CMake's probe but rejects it when compiling for ARM32.
ARM32 QEMU runs use a PID namespace to stay within Bionic's 16-bit PID limit.
Gitleaks history and worktree-diff scans are also clean; the history scan
mounts the worktree's common Git directory as well as the worktree.
DVSI inputs, audio, comparison reports and benchmark logs stay local.

### Subframe synthesis

The codec paths' synthesis in 10 ms subframes ([DVSI parity](dvsi-parity.md#synthesis-in-subframes))
renders each pitch pulse's harmonics by rotation in scalar code, where the rendering it replaces
there used the vector oscillator kernels. Measured as above (Ryzen 9
5950X, GCC 16.2.1, Release with fast-math and LTO on both sides, `bench_replay hard <codec> 10 5`
on DVSI's `dam` stream, the median of five alternating comparisons pinned to CPU 6), against main
at `8875729`. Portable is `MBELIB_ENABLE_NATIVE=OFF`, whose AVX2 kernel copies run on this CPU.

| Complete frame decode | Build | Before | After | Time added |
| --- | --- | ---: | ---: | ---: |
| D-STAR, hard | portable | 4.807 us | 5.814 us | 20.9% |
| AMBE+2, hard | portable | 5.174 us | 6.388 us | 23.5% |
| IMBE, hard | portable | 5.117 us | 6.140 us | 20.0% |
| D-STAR, hard | host-tuned | 4.363 us | 5.108 us | 17.1% |
| AMBE+2, hard | host-tuned | 4.662 us | 5.599 us | 20.1% |
| IMBE, hard | host-tuned | 4.623 us | 5.463 us | 18.2% |

### Clean soft blocks and noise generator lanes

Two changes that compute what main computes, against main at `361e0c1`:

- The soft ECC search returns a received word that is already the codeword of
  the hard decision's data bits before it builds its cost tables. On a clean
  stream that is every block.
- A frame draws DVSI's noise in eight lanes of 16-bit states instead of one
  chain of 240 dependent steps (`mbe_dvsi_noise_frame()`).

Measured on a Ryzen 7 5800HS with GCC 16.2.1, Release with fast-math and LTO on
both sides: `bench_replay hard|soft <codec> 10 5` on DVSI's `dam` (ProVoice: the
same speech encoded to IMBE parameter rows and reframed), the median of seven
alternating baseline/candidate invocations pinned to one CPU. Soft replay uses
confidence 255. Portable is `MBELIB_ENABLE_NATIVE=OFF`.

| Complete frame | Build | Before | After | Time reduction |
| --- | --- | ---: | ---: | ---: |
| D-STAR, hard decode | portable | 6.475 us | 6.194 us | 4.3% |
| AMBE+2, hard decode | portable | 7.163 us | 6.919 us | 3.4% |
| IMBE, hard decode | portable | 6.882 us | 6.641 us | 3.5% |
| ProVoice, hard decode | portable | 7.444 us | 7.205 us | 3.2% |
| D-STAR, soft decode | portable | 7.212 us | 6.444 us | 10.6% |
| AMBE+2, soft decode | portable | 7.882 us | 7.135 us | 9.5% |
| IMBE, soft decode | portable | 8.695 us | 6.977 us | 19.8% |
| ProVoice, soft decode | portable | 9.070 us | 7.521 us | 17.1% |

The host-tuned build takes 2.8-5.2% less for a hard decode and 11.1-22.5% less
for a soft decode; Clang 23.1.5 with fast-math and LTO at `-march=x86-64-v3`
takes 5.0-6.0% and 10.8-19.5% less. With errors in the stream fewer blocks
arrive clean: on `dam_e5_hd` with mixed confidences a soft P25 decode takes
6.0% less and a D-STAR one 4.3%. Encoding and TETRA decoding, which neither
change touches, are the controls: within 0.5% with GCC, and within 1.5% with
Clang, whose encode timings scatter more.

Both changes are integer arithmetic and return what main returns for every
input: the search returns the word its scan would have returned, and the lanes
hold the generator's own states, which a frame converts to the same samples
(`test_ecc`, `test_unvoiced_subframes`). What an edit can still move is how a
compiler treats the float code around it, so main and the change were built
alike in eight configurations (GCC portable and host-tuned with fast-math and
LTO, GCC with IEEE math and LTO, Clang with fast-math and LTO at
`-march=x86-64-v3`, and x87 arithmetic in a 32-bit build and in a 64-bit build
with `-mfpmath=387`, each at `-O3` and at `-O2`) and compared in each:

- The code of the shared library (the four builds with LTO) or of the
  library's objects (the others), function by function, with alignment
  padding, addresses and pc-relative displacements set aside, so not the
  constants an instruction loads: the soft search differs, and so do the
  functions that draw or skip a frame's noise and the one that holds the
  unvoiced subframe synthesis, the last in no float arithmetic, only in the
  integer conversions of the noise samples and the noise state. No other
  function does, except that the host-tuned build keeps two spilled vector
  registers of the D-STAR encoder in exchanged stack slots.
- All 321 DVSI frame streams and four ProVoice streams: 313 decoded with hard
  and with soft decisions (confidences running through every value, zero
  included) and DVSI's 12 soft-decision vectors with their own. These 638
  decodes, 695,558 frames and 111,289,280 samples give the same float samples,
  parameter bits, return values and `mbe_process_result` bit for bit, and the
  same 64-bit hash of the three parameter sets after every frame.
- The three encoders give the same parameter rows for all 97 DVSI inputs
  (296,772 frames).

A program that links the static library with LTO is optimized as a whole, and
there an edit can move the compiler's choices in code far from it: the
host-tuned `test_golden_pcm` specializes the FFT's twiddle setup on one side
only, where a division by the transform size, a power of two, becomes a
multiplication. The decode comparison was therefore repeated with the harness
linked that way in the four LTO builds, with the same result, and `gen_golden`
prints the same hashes on both sides in each. Output identity is thus what
these comparisons observed, not a property every toolchain is bound to keep.
The scalar, SSE2 and AVX2 Debug golden hashes, which the suite checks, are
unchanged.

Left out, each measured on the same setup:

- An encoder that sends a frame with no voiced band on its fixed grid taking
  again only the magnitude estimate it reads there, where
  `mbe_analysis_measure_at()` takes every per-harmonic measure (on `dam` 61%
  of rate-33 frames are such frames): 2.7-3.9% less to encode, and the same
  rows on all 296,772 DVSI frames in all eight configurations. But x87
  rounding follows the compiler's inlining and register allocation, which the
  change moves: in the 64-bit `-O2 -mfpmath=387` build the D-STAR encoder
  sends other bits near decisions, on a voiced frame that is never measured
  again among them.
- The TETRA synthesis filter with its history in registers: 13% less for a
  TETRA decode with GCC and 1% more with Clang, whose loop is the faster to
  begin with (4.0 us a frame against 4.8).
- An unvoiced subframe with no noise band skipping its forward transform and
  shaping: 1.3-1.9% less for a hard decode, and the same samples on all the
  streams above with portable GCC, host-tuned GCC and Clang. It does not raise
  the inexact flag main raises there, and a model built by hand with a NaN
  amplitude, which no codec path builds, renders differently.
- The encoder's front-end filters with their state in locals: 3-5% less to
  encode, and identical rows on all 296,772 DVSI frames, yet with portable GCC
  five frames of a quiet 130 Hz tone in noise changed 3 D-STAR, 7 AMBE+2 and 47
  IMBE bits: the vectorized loop on main rounds otherwise, and a decision of
  the noise estimator turns on it.
- One codebook row per SIMD lane in the joint AMBE vector search (7% less to
  encode) and a band fit for pitch refinement without the outputs it discards
  (3-5% less): identical on the DVSI inputs, different rows near decisions.
- A memo of the previous frame's model per thread (6% less for a hard decode):
  it depends on the rounding mode staying the same between calls, and two
  streams decoded alternately on one thread lose it and pay 2.5%.
- The response window of a pitch pulse by rotation instead of a cosine a
  sample, and the low-band response as a rational function without `atan2`:
  3-6% and 7-8% less for a hard decode, with PCM differing by up to 1.3e-5 and
  1.2e-6 of full scale.

Identical rows on a corpus do not show that a change to float code is exact
where the compiler may round one expression differently in two places
(fast-math, x87): a corpus seldom lands within a rounding error of a decision.
The encoder candidates were therefore also compared at such decisions. For
random pairs of three-frame inputs (noise, harmonic series, resonant noise,
voiced onsets) the blend between the two is bisected to where main's bits
change, and both builds encode every blend visited. In the x87 `-O2` build
above, where the corpus rows are identical, that probe finds the fixed-grid
estimates out within its first pair; in the five configurations it was run in
before, it passed them (about 145,000 inputs per encoder in each 64-bit one).
It does not aim at a decision that acts through encoder state some frames
later, and it passed the front-end filters, which an input built for one
decision then broke. So passing it in one configuration says nothing of
another, and passing it in all would not show exactness.

### Phase regeneration and pitch

The codec paths regenerate voiced phases from each frame's decoded `log2Ml`
([DVSI parity](dvsi-parity.md#phase-regeneration)) in a scalar loop of their
own (`mbe_regenerate_phases_log2()`), where they used the vector-batched kernel
that `mbe_synthesizeSpeechf()` keeps; P25 also computes its pitch offset in
double for each frame ([DVSI parity](dvsi-parity.md#p25-pitch)). Measured on a
Ryzen 7 5800HS with GCC 16.2.1, Release with fast-math and LTO on both sides,
portable (`MBELIB_ENABLE_NATIVE=OFF`): `bench_replay hard <codec> 10 5` on
DVSI's `dam`, the median of seven alternating invocations pinned to CPU 6, the
branch's code at `3e0c6ce` (with the review fixes: the D-STAR path's explicit
amplitudes, the overlap guards) against main at `7fb5ca6`. An earlier
measurement of `37f899a`, before those fixes, gave 3.0%, 3.3% and 4.0%.

| Complete frame decode | Build | Before | After | Time added |
| --- | --- | ---: | ---: | ---: |
| D-STAR, hard | portable | 6.066 us | 6.263 us | 3.2% |
| AMBE+2, hard | portable | 6.733 us | 6.964 us | 3.4% |
| IMBE, hard | portable | 6.449 us | 6.667 us | 3.4% |

The seven invocations ran from 6.047 to 6.418, 6.713 to 6.745 and 6.431 to
6.489 us on main and from 6.239 to 6.292, 6.939 to 7.008 and 6.662 to 6.705 us
now. The ranges do not overlap for AMBE+2 and IMBE. For D-STAR they do, by one
invocation of main at 6.418 us, above every invocation of the branch; main's
other six ran from 6.047 to 6.070 us, below all of the branch's.

### Tones and the octave step

The tone model change adds three things to a decode ([DVSI
parity](dvsi-parity.md#november-2026-tones-and-the-octave-step)): the octave
step of the pulse train, a test of every P25 voice frame against the P25 tone
table (`imbe_tone_word()`), and the tone model of the frames that pass. The
tone test first looks the frame's pitch code up in an index of the table's
rows (`imbe_tone_first_row`): a code with no row returns at once, and the scan
covers only the rows of the frame's code. Measured as for the phase
regeneration above, `bench_replay hard <codec> 60 5` on `dam` (which holds no
tone), seven alternating invocations pinned to CPU 6, medians, Release with
SIMD, fast-math and LTO, portable: the head of PR #196 (`fac0111`) against the
change set before the index (`b58d91d`) and with it (`3d84280`).

| Complete frame decode | `fac0111` | `b58d91d` | `3d84280` | Change set |
| --- | ---: | ---: | ---: | ---: |
| D-STAR, hard | 6.389 us | 6.351 us | 6.354 us | -0.5% |
| AMBE+2, hard | 7.053 us | 7.086 us | 7.074 us | +0.3% |
| IMBE, hard | 6.768 us | 6.863 us | 6.777 us | +0.1% |

The seven invocations ran from 6.269 to 6.424, 6.986 to 7.161 and 6.696 to
6.801 us on `fac0111`, from 6.299 to 6.388, 6.997 to 7.198 and 6.776 to 6.914
us on `b58d91d`, and from 6.309 to 6.472, 6.998 to 7.187 and 6.709 to 6.807 us
with the index. With the index every codec's median is within 0.6% of
`fac0111`'s (D-STAR -0.55%, shown rounded in the table) and the ranges
overlap. Before the index IMBE's median was 1.4%
higher, with ranges that also overlapped (6.696 to 6.801 against 6.776 to
6.914 us), so this measurement does not isolate the cost of the table scan.
The index leaves every decode unchanged: the scoreboard's hashes
of all 255 decodes of DVSI's bits on the development, validation, tone and
channel-error partitions, and of the 75 encodings of the speech vectors, are
those of `b58d91d`.

### Mid-frame voicing

The mid-frame voicing change ([DVSI
parity](dvsi-parity.md#november-2026-mid-frame-voicing)) adds the P25 onset
columns to the mid-frame model of a P25 or ProVoice frame after a frame with no
voiced band, and a row lookup to that of rate 33's V/UV rows 13 and 17.
Measured with `bench_replay hard <codec> 60 5` on the `dam` stream of each
codec, seven alternating invocations pinned to CPU 6, a portable Release build
with SIMD, fast-math and LTO, `cc272ef` (the tone model change before its
review's fixes) against the change set on it, before the change set was moved
onto `c2ba220`, on a machine busy with other work (load 1.5 to 2.5). Medians,
with the range of the seven invocations:

| Complete frame decode | `cc272ef` | Mid-frame voicing | Change |
| --- | ---: | ---: | ---: |
| D-STAR, hard | 6.612 us (6.533 to 6.657) | 6.633 us (6.587 to 6.869) | +0.3% |
| AMBE+2, hard | 7.373 us (6.979 to 7.832) | 7.279 us (6.994 to 7.380) | -1.3% |
| IMBE, hard | 6.748 us (6.715 to 7.201) | 6.725 us (6.706 to 6.753) | -0.3% |

The ranges overlap in every codec: no cost is measurable this way.

### Rate 33's mid-frame pitch wave

Rate 33's mid-frame pitch wave ([DVSI
parity](dvsi-parity.md#november-2026-rate-33s-mid-frame-pitch-wave)) adds, to
the mid-frame model between two voiced rate-33 frames of different pitch, one
`log2`, `cos` and `sin` in double for each of the two frames. Measured with
`bench_replay hard <codec> 60 5` on the `dam` stream of each codec, seven
alternating invocations pinned to CPU 6, a portable Release build with SIMD,
fast-math and LTO (GCC 16.2.1, Ryzen 7 5800HS), `dac6d14` against `d08ea39`,
on a machine busy with other work (load 1.7 to 2.0). D-STAR and IMBE, which
the change does not reach, are controls. Medians, with the range of the seven
invocations:

| Complete frame decode | `dac6d14` | Mid-frame pitch wave | Change |
| --- | ---: | ---: | ---: |
| AMBE+2, hard | 7.174 us (7.045 to 7.384) | 7.238 us (7.061 to 7.529) | +0.9% |
| D-STAR, hard | 6.438 us (6.375 to 6.817) | 6.411 us (6.356 to 6.554) | -0.4% |
| IMBE, hard | 6.884 us (6.809 to 7.248) | 6.911 us (6.841 to 7.410) | +0.4% |

The ranges overlap in every codec: no cost is measurable this way. The D-STAR
and IMBE checksums are the same in both builds, the AMBE+2 ones differ.

### Mid-frame voicing by columns

The mid-frame columns change ([DVSI
parity](dvsi-parity.md#november-2026-mid-frame-voicing-by-columns)) samples each
of the two frames' voicing at eight 500 Hz columns and reads each changing band's
two columns, in D-STAR's and P25's renderings. Measured with `bench_replay hard
<codec> 60 5` on the `dam` stream of each codec, seven alternating invocations
pinned to CPU 6, a portable Release build with SIMD, fast-math and LTO (GCC
16.2.1, Ryzen 7 5800HS), `c3c6d42` against `e9cf77c`, on a machine busy with
other work (one-minute load 1.7 to 2.1). AMBE+2, which the change does not
reach, is a control. Medians, with the range of the seven invocations:

| Complete frame decode | `c3c6d42` | Mid-frame columns | Change |
| --- | ---: | ---: | ---: |
| D-STAR, hard | 6.410 us (6.368 to 7.118) | 6.537 us (6.481 to 7.016) | +2.0% |
| IMBE, hard | 7.236 us (7.006 to 7.331) | 7.298 us (7.186 to 7.455) | +0.9% |
| AMBE+2, hard | 7.443 us (7.196 to 7.655) | 7.321 us (7.171 to 7.552) | -1.6% |

The ranges overlap in every codec; the control moves by 1.6%, so a cost of the
D-STAR and IMBE changes is not resolved this way. The AMBE+2 checksums are the
same in both builds, the D-STAR and IMBE ones differ.

### The pulsed state

The pulsed state ([DVSI
parity](dvsi-parity.md#november-2026-the-pulsed-state)) tests each model with no voiced band of a
D-STAR or P25 frame for a reserved code (a range check on the pitch, and the decoder's own pitch
function only when the pitch lies in the codes' range: most such models have no code), in P25 also
for the pitch code that steps the validity counter (the decoder's pitch function for every pitch in
the range of the 208 codes), and renders the pulses and the pulsed bands'
setup where a model has them. Measured with `bench_replay hard <codec> 60 5` on the `dam` stream of each
codec, seven alternating invocations pinned to CPU 6, a portable Release build with SIMD, fast-math
and LTO (GCC 16.2.1, Ryzen 7 5800HS). AMBE+2's output is the same in every build (its checksums
agree), but its timed path runs the changed code too (`prepare_frame()`, `pulsed_frame()`, three
more `prepare_render()` calls and the state store): it is not an untouched control. Medians, with
the range of the seven invocations. First, `74126a3` against `d3fb70f` (one-minute load 1.5 falling
to 1.1):

| Complete frame decode | `74126a3` | Pulsed state | Change |
| --- | ---: | ---: | ---: |
| D-STAR, hard | 6.454 us (6.445 to 6.482) | 6.816 us (6.791 to 6.880) | +5.6% |
| IMBE, hard | 6.864 us (6.850 to 6.959) | 7.246 us (7.211 to 7.307) | +5.6% |
| AMBE+2, hard | 7.032 us (7.013 to 7.056) | 7.092 us (7.088 to 7.241) | +0.9% |

The ranges do not overlap in any of the three; AMBE+2's +0.9% has no examined cause. Under
callgrind, the P25 `dam` stream run once in the same builds executes 265.4 million instructions
before and 272.1 million (+2.5%) after: most of the difference is the rendering of the emitted
pulses and the setup of the pulsed bands' oscillators, and the code tests; why the time rises more
than the instructions was not examined.

After the review, twice, the three builds alternating, `74126a3`, `d3fb70f` and the final code with
P25's validity counter (one-minute load about 1.0); each cell the two runs' medians:

| Complete frame decode | `74126a3` | `d3fb70f` | Final | Change, final |
| --- | ---: | ---: | ---: | ---: |
| D-STAR, hard | 6.442, 6.454 us | 6.695, 6.697 us | 6.682, 6.692 us | +3.7%, +3.7% |
| IMBE, hard | 6.851, 6.849 us | 7.121, 7.146 us | 7.227, 7.180 us | +5.5%, +4.8% |
| AMBE+2, hard | 7.022, 7.014 us | 7.032, 7.014 us | 7.020, 7.016 us | 0.0%, 0.0% |

Against `74126a3` the ranges do not overlap in D-STAR or IMBE and do in AMBE+2, whose +0.9% of the
first run does not recur. Against `d3fb70f` the final code's ranges overlap in D-STAR and AMBE+2;
in IMBE the first run has it 1.5% above with the ranges apart (7.100 to 7.146 against 7.175 to
7.252 us) and the second 0.5% above with them overlapping (7.109 to 7.190 against 7.162 to 7.318
us): these give no bound on the counter's own cost. The checksums of `d3fb70f` and the final code are the
same in all three codecs.

### Rate 33's pulsed rows

Rate 33's pulsed rows ([DVSI parity](dvsi-parity.md#november-2026-rate-33s-pulsed-rows)) read a
frame's word and the last frame's slot, test the band-grid model, and render up to three pulses
where a frame of V/UV rows 18 to 31 has them; every subframe rendering runs `pulsed_rows()` and the
slot's store. Measured with `bench_replay hard <codec> 60 5` on DVSI's `dam` stream of each codec
(`mbe_quality_reframe --from-dvsi`), seven alternating invocations pinned to CPU 6, a portable
Release build with SIMD, fast-math and LTO (GCC 16.2.1), `2b6b6c2` against `1bc657f`, twice, with a
browser open (one-minute load 1.2 to 1.5). Medians, with the range of the seven invocations' medians:

| Complete frame decode | `2b6b6c2` | Pulsed rows | Change |
| --- | ---: | ---: | ---: |
| AMBE+2, hard | 7.308 us (7.162 to 7.592); 7.153 us (7.126 to 7.544) | 7.420 us (7.369 to 7.533); 7.411 us (7.358 to 7.792) | +1.5%; +3.6% |
| D-STAR, hard | 6.847 us (6.801 to 6.890); 6.893 us (6.815 to 7.102) | 6.906 us (6.826 to 7.094); 6.882 us (6.808 to 7.084) | +0.9%; -0.2% |
| IMBE, hard | 7.479 us (7.337 to 7.785); 7.392 us (7.353 to 7.571) | 7.390 us (7.305 to 7.823); 7.366 us (7.308 to 7.493) | -1.2%; -0.4% |

The ranges overlap in every row; D-STAR's and IMBE's checksums are the same in both builds (their
outputs equal), AMBE+2's differ (its decode changes). The runs give no bound tighter than that
spread on the change's cost.

### DVSI's low band

DVSI's low band ([DVSI parity](dvsi-parity.md#november-2026-dvsis-low-band)) gives every voiced and
pulsed harmonic of D-STAR's, P25's and rate 33's renderings the third-order filter's gain and phase
(three double-precision functions per harmonic of each model, against one before) and the cut's
comparison, and every frame the tracker's step. Measured with `bench_replay hard <codec> 60 5` on
DVSI's `dam` stream of each codec, seven alternating invocations pinned to CPU 6, a portable Release
build with SIMD, fast-math and LTO (GCC 16.2.1), `9a101f8` against `0e3463e` (one-minute load 1.6 to
4.2 during the run):

| Complete frame decode | `9a101f8` | Low band | Change |
| --- | ---: | ---: | ---: |
| AMBE+2, hard | 7.345 us (7.327 to 7.682) | 7.650 us (7.582 to 7.970) | +4.2% |
| D-STAR, hard | 6.809 us (6.791 to 7.078) | 7.100 us (7.037 to 7.192) | +4.3% |
| IMBE, hard | 7.279 us (7.231 to 7.574) | 7.528 us (7.492 to 7.560) | +3.4% |

The ranges overlap a little in each row; every checksum differs (every decode changes). Timed again
after the review's fixes (`c449a2d`, the same checksums as `0e3463e`), with the machine loaded by
other work (one-minute load 8 to 13), medians of the seven invocations' medians: AMBE+2 10.804 →
10.516 us (-2.7%), D-STAR 8.734 → 9.266 (+6.1%), IMBE 8.818 → 9.299 (+5.5%), the ranges overlapping
widely; under that load these say nothing finer than the first run.

### AMBE block lengths at 20 harmonics

The block lengths at 20 harmonics ([DVSI parity](dvsi-parity.md#november-2026-the-ambe-block-lengths-at-20-harmonics))
change one row of a table. Measured with `bench_replay hard <codec> 60 5` on DVSI's `dam`
stream of each codec and with `bench_encode`, seven alternating invocations pinned to CPU 6, a portable
Release build with SIMD, fast-math and LTO (GCC 16.2.1), `f725a66` against `59dddb3` (one-minute load 5.4
to 8.8 during the run, another agent's builds running):

| Median per frame | `f725a66` | Block lengths | Change |
| --- | ---: | ---: | ---: |
| AMBE+2 decode, hard | 8.167 us (8.067 to 8.478) | 8.157 us (8.096 to 8.269) | -0.1% |
| D-STAR decode, hard | 7.723 us (7.520 to 9.698) | 7.701 us (7.558 to 9.414) | -0.3% |
| IMBE decode, hard | 8.094 us (7.987 to 12.848) | 8.057 us (7.998 to 11.348) | -0.5% |
| D-STAR encode | 30.32 us (29.61 to 31.71) | 30.16 us (29.66 to 30.44) | -0.5% |
| AMBE+2 encode | 31.55 us (30.67 to 31.92) | 31.42 us (30.48 to 31.79) | -0.4% |
| IMBE encode | 15.18 us (14.77 to 15.25) | 15.16 us (14.75 to 15.19) | -0.1% |

Timed before the merge of the branch below's final head (`e345239`, whose decodes equal `f725a66`'s on every permitted vector). The ranges overlap in every row: no resolvable cost under that load, not a bound (the D-STAR encoder's -0.53% is the largest change); the D-STAR and AMBE+2 checksums differ (the `dam` streams hold frames of
20 harmonics), the IMBE ones are equal. An earlier run under heavier load (10 to 20) scattered by -6.7% to
+1.9%, the unchanged IMBE decode by -5.3%.

### Four rules of DVSI's P25 decoder

The four rules of DVSI's P25 decoder ([DVSI
parity](dvsi-parity.md#november-2026-four-rules-of-dvsis-p25-decoder)) replace the old law's three
double-precision functions per code with an integer quotient, add a few double-precision functions
per mid-frame model of two voiced P25 models, one double-precision logarithm per model for the phase
step, and the column reading of P25's two own models. Measured with `bench_replay hard <codec> 60 5`
on DVSI's `dam` stream of each codec, seven alternating invocations pinned to CPU 6, a portable
Release build with SIMD, fast-math and LTO (GCC 16.2.1), `f725a66` against each rule's commit (rules
2 to 4 measured on their first forms, which differ from the commits only in how `mid_pitch()` nests
its branches), taken under load (the one-minute load averages read during the four runs 10.88,
9.64, 7.71 and 9.67; another agent's builds and tests running):

| Complete frame decode | `f725a66` | Grid | Parabola | Phase step | Own columns |
| --- | ---: | ---: | ---: | ---: | ---: |
| IMBE, hard | 8.558 us | 8.358 (-2.3%) | 9.191 (+7.4%) | 8.213 (-4.0%) | 8.698 (+1.6%) |
| D-STAR, hard (control) | 8.171 us | 8.093 (-1.0%) | 7.774 (-4.9%) | 7.854 (-3.9%) | 8.664 (+6.0%) |
| AMBE+2, hard (control) | 8.478 us | 8.506 (+0.3%) | 8.420 (-0.7%) | 8.601 (+1.5%) | 8.448 (-0.4%) |

The controls' checksums are the same in every build (their outputs equal), IMBE's differ (its decode
changes). The controls, whose decodes do not change, moved by as much as IMBE (D-STAR -4.9% to
+6.0%), and an intermediate IMBE build by +7.4%: under that load the runs do not bound the change's
cost, which is unresolved. A source of extra work worth measuring under controlled load:
`p25_column_band_state()` asks both frames whether they are read by the columns for every harmonic
(a scan of each model's voicing and the tone model's test), though the answers are constant over the
frame; computing them once per frame moved P25's samples in the 32-bit Release build (x87), so it
was reverted. The fix round's additions (the tone model's test through the decoder, the pulse
train's grid test for each interval between two voiced models) were not timed.

## TETRA Conformance

EN 300 395-4 tests a TETRA speech decoder bit for bit: TE1_46.CHD (138 little-endian 16-bit
words per frame: the bad-frame indicator, then B1..B137) decoded must equal TE2_128.PCM. The
files are in `orignl.zip` inside ETSI's en_30039504v010301p0.zip. Keep them, and any build of
ETSI's reference decoder and its output, outside the checkout; never commit them.

```sh
cmake --preset dev-debug -DMBELIB_BUILD_TOOLS=ON
cmake --build --preset dev-debug -j --target mbe_tetra_decode
./build/dev-debug/mbe_tetra_decode TE1_46.CHD /tmp/ours.pcm
python3 -I tools/quality/tetra_conformance.py compare /tmp/ours.pcm TE2_128.PCM
```

`compare` exits 0 only for identical files, as TE1_46.CHD now gives; `--gate` instead takes the
exit status from tolerance thresholds, for runs that may differ (a candidate change, say);
`--json` writes the metrics. For concealment,
`tetra_conformance.py bfi TE1_46.CHD /tmp/lossy.chd --seed 1 --rate 3` (also `--burst START:LENGTH`
and `--first N`) marks frames bad; decode the copy with both decoders and compare. Measured
results are in [docs/tetra-acelp.md](tetra-acelp.md#measured-results).

## Speech Quality Evaluation

Developer-only tooling; no new library runtime dependency. CI builds the two
quality tools and tests their filesystem helper; full quality runs remain opt-in.
The C99 `mbe_quality_eval` uses the installed public API and an independent 256-point FFT.
`mbe_quality_reframe` is a statically linked private fixture generator. A/B runs
additionally require Bash, Python 3, git, gcc/g++ (C++17), `patch`, and
coreutils. Fetching held-out speech requires network access and `ffmpeg`.
Offline reports require NumPy, SciPy, matplotlib, and pystoi; install these in a
local virtual environment rather than making them library dependencies.

Quality-tool outputs are created with owner read/write permissions on POSIX
(`0600`, further restricted by the caller's umask). Reopening an output truncates
its contents while retaining its existing permissions. Windows uses CRT
read/write permissions together with the directory's inherited ACL.

### Windows

With `MBELIB_BUILD_TOOLS=ON`, `mbe_quality_eval` and `mbe_quality_reframe`
build with MSVC or MinGW for local single measurements and fixture generation.
The evaluator always links the shared library; static linking would embed one
implementation and prevent measuring a different library at runtime.

The A/B runner, calibration, and OP25 encoder build remain Linux-only. The
comparison workflow uses `LD_TRACE_LOADED_OBJECTS` to verify which shared
library was actually loaded; Windows DLL placement does not reproduce that
provenance check. Windows evaluator paths must be ordinary ASCII drive paths;
non-ASCII names, UNC, `\\?\` and `\\.\` namespaces, alternate data streams and
leaves ending in a dot or space are rejected.

### Reproducible comparison workflow

```sh
cmake --preset dev-debug -DMBELIB_BUILD_TOOLS=ON
cmake --build --preset dev-debug --parallel
tools/quality/build_op25_encoder.sh
tools/quality/fetch_speech_corpus.sh --set legacy
python3 -m venv build/quality/venv
build/quality/venv/bin/python -m pip install numpy scipy matplotlib pystoi

# Preserve source revisions in external worktrees; do not reuse a different revision.
git worktree add --detach /tmp/mbelib-quality-original 5fd3f3dd738b76bc9a8446ff4d7ae99e565e35e0
git worktree add --detach /tmp/mbelib-quality-checkpoint 8da251d105e50815514b1eee994c2c8471469654
cmake -S /tmp/mbelib-quality-original -B build/quality/original-debug \
  -DCMAKE_BUILD_TYPE=Debug -DMBELIB_BUILD_TESTS=OFF -DMBELIB_BUILD_EXAMPLES=OFF
cmake --build build/quality/original-debug --parallel
cmake -S /tmp/mbelib-quality-checkpoint -B build/quality/checkpoint-debug \
  -DCMAKE_BUILD_TYPE=Debug -DMBELIB_BUILD_TESTS=OFF -DMBELIB_BUILD_EXAMPLES=OFF
cmake --build build/quality/checkpoint-debug --parallel

tools/quality/run_quality_ab.sh build/quality/original-debug build/quality/original-debug \
  --out-dir build/quality/runs/identity-debug
tools/quality/run_quality_ab.sh build/quality/original-debug build/dev-debug \
  --out-dir build/quality/runs/original-to-final-debug
tools/quality/run_quality_ab.sh build/quality/checkpoint-debug build/dev-debug \
  --calibration-libdir build/quality/original-debug \
  --out-dir build/quality/runs/checkpoint-to-final-debug
build/quality/venv/bin/python tools/quality/analyze_quality.py report \
  --run-dir build/quality/runs/original-to-final-debug
```

The runner accepts `baseline_libdir candidate_libdir [corpus_dir]`, followed by
`--out-dir`, `--evaluator`, `--reframer`, `--calibration-libdir`, and
`--gain-profile`. Without `--out-dir`, it creates a unique directory under
`build/quality/runs/`; an explicit nonempty output directory is rejected.
Historical `build/quality/out/` and `build/quality/frames/` are not overwritten.
Current cache entries live under `build/quality/frame-cache-v3/`.

Both libraries must have compatible installed signatures/state layouts and
matched compiler/DSP options. Loader tracing checks the actual selected library,
including inode identity, and rejects `LD_PRELOAD` and unexpected extra libmbe
libraries. A sanitizer evaluator must use its instrumented library for both
operands, never loader-swap an uninstrumented baseline into that process.

Each schema-2 `manifest.json` binds library/tool/build identities, seed, flushing,
references, encoded/reframed bits, calibration records, WAVs, metric JSON, command
logs, and tables. Inputs and reports are copied into the run directory and indexed
by SHA256; analysis does not depend on the mutable frame cache. Source revision
labels require a captured matching binary hash; a checkout or directory name
alone does not prove binary provenance. Unknown source identities remain explicit.
`build/quality/checkpoint-capture.json`, when present, supplies these SHA-bound
comparison points. Same-library runs require exact PCM, frame-count, and complete
metric/null/support equality. Original-only identity controls can pass that gate
without establishing the later corrected-tone/adaptive regression gates.

### Encoder provenance, frame paths, and calibration

- Codec2 revision: `310777b1c6f1af0bc7c72f5b32f80f6fd9136962`
  (<https://github.com/drowe67/codec2>). The named **legacy/development** set is
  `hts1a`, `hts2a`, `hts1`, `kristoff`, and `ve9qrp_10s`. `hts1a` is exactly the
  first 24,000 samples of `hts1`; both remain available, but the duplicate is
  excluded from independent aggregate statistics when its parent is present.
- OP25 revision: `71abcd0ead32f86f51615ea6cc8a6a4dba4c949a`
  (<https://github.com/boatbod/op25>). Its GPL-3.0-or-later encoder and bundled
  mbelib remain in a separate executable. Never mix their incompatible
  `mbe_parms` with mbelib-neo. `op25_mbelib_gain.patch` changes disposable build
  copies only. Its AMBE2400 bridge is a decoder-targeted fixture, not an
  interoperable D-STAR transmitter or independent wire-conformance proof.
- IMBE7200 is encoded/calibrated once, then its canonical 88-bit parameters are
  reframed into both 184-bit IMBE7200 and 168-bit IMBE7100 rectangular frames.
  AMBE2450's 49-bit parameters become 96-bit frames; OP25's existing 96-bit
  AMBE2400 output is retained. Reports exercise all four actual `Framef` APIs.
  Canonical IMBE `Dataf` is an exact-PCM equality control, not a fifth independent
  codec. Reframed IMBE7100 is not a native transmitter/capture validation.
- `--gain-adjust` is log2 encoder attenuation, default `log2(7)`. Per-clip
  calibration uses only the calibration library, targeting decoded active RMS
  −3 ± 0.5 dB relative to reference, with at most six passes and lag 0…640.
  Schema-2 speech support, finite alignment/level, full reference coverage,
  and zero nonfinite float samples are mandatory. Invalid calibration fails;
  silent controls bypass calibration rather than selecting a fallback gain.
- `--flush-frames N` accepts integers 0…50, default 5. After the final nonempty
  input block, the encoder submits exactly N zero PCM blocks through its normal
  output path; ordinary partial-block padding remains. Empty input fails.
  Original sample counts and flushing are separately recorded and cache-bound.
  Five blocks cover alignment without discarding the reference tail; this is
  encoder fixture completion, not added decoder latency.

### Evaluator interface and measurements

```sh
build/quality/op25_encode --mode imbe7200 \
  --in build/quality/corpus/hts1a.raw --out build/quality/example.data.txt \
  --gain-adjust 3.4 --flush-frames 5
build/dev-debug/mbe_quality_reframe --codec imbe7100 \
  --in build/quality/example.data.txt --out build/quality/example.frames.txt
build/dev-debug/mbe_quality_eval --codec imbe7100 \
  --frames build/quality/example.frames.txt --out build/quality/example.wav \
  --ref build/quality/corpus/hts1a.raw --seed 0x12345678 \
  --json build/quality/example.json
build/dev-debug/mbe_quality_eval --decoded build/quality/example.wav \
  --ref build/quality/corpus/hts1a.raw --lag-samples 456 \
  --json build/quality/example.pcm.json
```

Use baseline calibration or a frozen profile for comparisons, not the example
gain. Frame rows are literal ASCII 0/1, including rectangular padding:
IMBE7200 accepts 88-bit Dataf or 184-bit Framef; IMBE7100 only 168-bit Framef;
either AMBE mode accepts 49-bit Dataf or 96-bit Framef. One input file uses one
selected public API, reported as `input_api`. The reframer accepts parameter
rows only and validates everything before changing its destination.

`--decoded input.raw|input.wav` is measurement-only: it rejects `--codec`,
`--frames`, `--out`, and `--seed`, and never rewrites PCM. Raw/WAV inputs must
be mono 8-kHz s16le; WAV chunks are parsed rather than assuming a fixed header.
Output WAVs preserve pre-normalization PCM. Invalid input, aliases (including
hard links and normalized not-yet-created output paths), or fewer than 256
aligned samples fail with exit 2 **before output creation**. A negative codec
status or nonfinite float PCM exits 3 with the frame index.

The 16-bit PCM written and measured comes from the codec's own public `short`
entry points, run over the same rows from a second, separately initialized and
seeded state after the float pass, so it is what a caller of those entry points
gets (the D-STAR and AMBE 3600x2450 ones round to nearest, the IMBE ones
truncate; [DVSI parity](dvsi-parity.md#16-bit-output)). The float pass keeps
supplying the statistics, the float diagnostics below and the `--params` dump.
The two passes must agree on every frame's status, flags, error counts and
parameter bits, or the tool exits 3. Up to 2.6.8 the tool converted the float
output with `mbe_floattoshort()`.

Schema-2 metrics preserve the original supported formulas:

- `lag_samples` is the lag actually used. Optional `--lag-samples` accepts
  −160…800. Automatic alignment always runs independently: `auto_lag_samples`,
  nullable `alignment_corr`, and `alignment_at_limit` describe its estimate.
  It correlates log-RMS envelopes from centered 160-sample windows at 8-sample
  hops. Correlations within `1e-12` are tied, favoring the smallest absolute lag.
  Positive lag means decoded speech is delayed. The runner measures the
  candidate at the baseline lag, preserving common reference support/activity
  and the synthesis-join grid; differing automatic lags generate a separately
  labelled auto-alignment sensitivity report, not a claimed decoder-delay change.
- `reference_samples`, `decoded_samples`, `aligned_samples`, both signals'
  `*_trim_head`/`*_trim_tail`, `active_frames`, `spectral_frames`, and `join_count`
  expose support. Active complete 160-sample reference frames have positive
  energy within 40 dB of maximum. Boundary alignment or incomplete reference
  coverage invalidates a comparison.
- `reference_state` is exactly `no_reference`, `silent_reference`,
  `silent_decoded`, or `speech`. Speech requires positive active energy in both
  signals. Undefined reference metrics/correlation serialize as JSON `null`,
  never favorable zero LSD or zero crest error. Pure silence succeeds as an
  explicitly inapplicable speech-quality measurement.
- `level_offset_db` is decoded/reference active RMS before normalization.
  Reference-based measurements then use active-RMS-matched int16-scale signals.
  `lsd_db` is mean active-frame log-spectral distance: symmetric 256-point Hann,
  80-sample hop, bins 2…118, reference-derived power floor
  `max(1e-6 * max(reference power), 1e-20)`.
- `waveform_corr` is the normalized correlation of the two waveforms over the
  active frames, the best within 32 samples either side of the lag used
  (`waveform_shift_samples`), since that lag comes from envelopes at 8-sample
  hops and can land 16 samples from the waveforms' common timing. The other
  metrics see relative harmonic phase only indirectly (crest factor,
  short-window envelopes and spectra); this one sees the waveform itself, so two
  decoders with equal spectra and unrelated pulse positions score near 0, while
  a pulse train shifted by up to 32 samples as a whole still scores high.
  `waveform_corr_unshifted` is the same correlation with the two files sample
  against sample (the lag undone, no search), reported when the lag used is
  within 32 samples. Two decoders of the same bits share their timing, so
  between them a sample of delay, or a noise sequence one sample off, shows
  here and not in `waveform_corr`.
- `env_corr` correlates log 40-sample RMS envelopes, floored 60 dB below each
  maximum. `crest_ref_db`, `crest_dec_db`, and `crest_delta_db` are mean 20-ms
  peak/RMS values and decoded-minus-reference difference. `band_delta_db_*`
  measures active STFT power ratios in 0–500, 500–1000, 1000–2000, 2000–3000,
  and 3000–4000 Hz, plus overlapping 0–250 Hz and 3500–4000 Hz bands; the
  top one covers what `lsd_db` (up to 3687.5 Hz) leaves out. Without reference, decoded absolute
  `band_db_*` remains available where supported.
- `join_dec_db`/`join_ref_db` compare squared across-join derivatives with the
  eight neighboring derivatives on the same reference mask and actual retained
  synthesis grid; `join_excess_db` is their difference. `boundary_index_db`
  compares unshifted boundary-region and mid-frame derivative energy. These are
  **diagnostics, not standardized click detectors**: a continuous 75-Hz sine
  produces boundary index +16.61 dB and zero reference-matched join excess.
- `float_nonfinite_samples`, `float_peak`, and `float_clip_samples` inspect
  **post-library-limiter, pre-int16** PCM, using bitwise finite classification
  even under fast-math and the existing `(32767 * 0.95) / 7` threshold. They do
  not reveal inaccessible pre-limiter peaks. PCM-input mode reports them null.
  `pcm_peak`, `pcm_rail_samples` (`abs(sample) >= 31128`), and
  `pcm_max_rail_run` are measured before normalization. Tone/erasure/repeat/mute/
  silence frame counts and C0/protected/C4/total error sums retain process-result
  context; clips with any of them are not clean speech for the checkpoint gate.

### Held-out speech and fixed operating points

```sh
tools/quality/fetch_speech_corpus.sh --set heldout
tools/quality/run_quality_ab.sh build/quality/original-debug build/dev-debug \
  build/quality/corpus-heldout \
  --gain-profile build/quality/runs/original-to-final-debug/gain-profile.json \
  --out-dir build/quality/runs/heldout-debug
build/quality/venv/bin/python tools/quality/analyze_quality.py report \
  --run-dir build/quality/runs/heldout-debug
```

The fetcher preserves `arctic_a0001.wav`…`arctic_a0010.wav` for each of
`bdl`, `rms`, `clb`, and `slt` from
<http://festvox.org/cmu_arctic/cmu_arctic/>: 40 US-English recordings from two
male and two female speakers. Original mono16-kHz/s16 WAVs and each `COPYING`
remain under `build/quality/arctic-src/`; manifests bind source/converted hashes,
conversion argv and ffmpeg version. There is no loudness normalization.
`build/quality/corpus-heldout-minus12db/` contains
`round(sample * 10**(-12/20))` copies with preserved notices and native hashes.
Missing/invalid assets fail preparation; development speech is not substituted.

For a complete original-main native legacy calibration run, `report` emits
`gain-profile.json`: schema 1, calibration library SHA256, twelve source-manifest
hashes, and the per-mode median `gain_adjust_log2` across
`hts1`, `hts2a`, `kristoff`, `ve9qrp_10s`. Fixed-profile runs validate and preserve
these exact bytes, bypass per-clip calibration, and use the same gains across
native/−12 dB levels and Debug/Release. IMBE7100 inherits IMBE7200. Held-out
failures are findings, not permission to retune gains or phase constants.

Repeat the command with the attenuated corpus and a distinct output directory.
Release uses `--evaluator build/dev-release/mbe_quality_eval` and
`--reframer build/dev-release/mbe_quality_reframe`, matched Release library
operands, and the same profile; encoded-frame hashes must match Debug.

### Controls, acceptance, and listening

```sh
build/quality/venv/bin/python tools/quality/analyze_quality.py controls \
  --evaluator build/dev-debug/mbe_quality_eval --out-dir build/quality/control-run
build/quality/venv/bin/python tools/quality/analyze_quality.py benchmark \
  --baseline build/quality/original-release/bench_synth \
  --candidate build/dev-release/bench_synth --out build/quality/benchmark-original-final.json
```

Controls exercise the real PCM-input evaluator: identity, exact ×2 gain, ±64
sample delays, silence/silence, speech/zero, too-short input, aliases, clipping,
continuous 75-Hz sine, and increasing alternating frame-boundary steps. NumPy
independently cross-checks LSD on controls and real speech. The benchmark
alternates seven process runs per variant, records stdout/build/binary/library
identities, and compares median reported `avg` CPU seconds; the cost gate is
≤10% increase, not an algorithmic-delay measurement.

`report` validates immutable snapshots and writes `analysis.json`,
`acceptance.json`, per-clip/grouped CSV tables, shared-support waveform and
reference-scaled spectrogram figures, harmonic detail, band differences, and
zoomed joins for every flagged clip. STOI/ESTOI use the complete aligned overlap,
not concatenated active snippets, and are intelligibility diagnostics—not MOS
or naturalness estimates.

Acceptance separates `correctness`, `objective_nonregression`, and `perceptual`,
each `pass`, `fail`, or `not_established`, for every actual mode:

- Correctness requires supported/full-coverage measurements, no nonfinite PCM,
  exact canonical/frame equivalence, same-library identity where applicable,
  and candidate-bound native regression evidence for corrected tones,
  nonnegative attenuation, and warmed FFT/WOLA. Checkpoint comparisons also
  require unchanged clean-speech PCM.
- Objective screens use independent per-mode clip aggregates:
  mean absolute candidate-minus-baseline LSD change ≤0.3 dB; mean envelope
  correlation decrease ≤0.01; absolute mean per-band change ≤0.5 dB; candidate
  mean join excess ≤`max(0, baseline mean) + 0.1` dB; and **mean per-clip absolute
  crest error must not increase**. Absolute signed-mean crest error is reported
  separately. Every per-clip exceedance, worst clip, rail rate/run and increase
  remains visible even when aggregate screens pass.
- Listening material is reference/A/B with shared timing, equal active RMS and
  one common extra attenuation if needed to avoid clipping. Opaque presentation
  filenames and deterministic seed `0x12345678` randomize labels. The answer key
  is under `listening/private/`; `listening/listener.csv` contains item IDs and
  blank `listener_id,preference,naturalness_a,naturalness_b,intelligibility_a,
  intelligibility_b,artifact_notes` fields. Do not expose the key to listeners.
- Each held-out run prepares twelve balanced items per codec (each speaker's
  utterances 1,2,3); legacy per-clip exceedances are separate diagnostic items.
  At least five actual listeners must each rate all twelve items before a mode
  claim. Preferences are `A`, `B`, or `tie`; ratings are integers 1…5.
  Run `report --run-dir DIR --listening-results FILE` after collection.
  Ten thousand seeded bootstrap replicates independently resample listeners
  and clips, with clip resampling stratified by speaker. Preference CI must be
  strictly positive and intelligibility CI not wholly negative for a pass;
  wholly negative preference/intelligibility means fail, otherwise preference
  is not established. Shared IMBE7100/7200 PCM does not inflate independent
  sample size. This is a small blinded engineering comparison, not an
  ITU-compliant MOS/MUSHRA study. Playback alone is not a judgment.

### DVSI reference vectors

DVSI publishes AMBE-3000 test vectors on
<https://www.dvsinc.com/dlapps/appsoft.shtml>: input speech, DVSI-encoded bit
streams and DVSI-decoded speech for D-STAR, P25 and AMBE+2 (rate 33), plus
channel-error variants. They are the closest available stand-in for real
hardware, and they let decoder output be compared with DVSI's own decoder on
identical bits.

The [site terms](https://www.dvsinc.com/terms_of_use.shtml) grant no license to
copy, modify, transfer or mirror its materials. Local measurement is the tool's
intended use; it does not establish permission. Never commit the vectors,
attach them to issues or pull requests, or use DVSI audio in published
listening material. Only URLs, SHA-256 digests and aggregate numbers belong in
this repository.

```sh
# about 223 MB of the 1.9 GB archive; a directory outside the checkout survives
# rebuilds and can be shared by several checkouts
python3 tools/quality/fetch_dvsi_vectors.py --out ../dvsi-vectors
build/dev-debug/mbe_quality_reframe --codec ambe2400 \
  --from-dvsi ../dvsi-vectors/tv-rc/dstar/dam.bit --out dam.frames
```

- The fetcher needs `curl`. It reads only the needed members of
  `get-usb/tv-rc.zip` through HTTPS range requests. It verifies each CRC-32 and writes
  `dvsi-manifest.json` with SHA-256 digests, retrieval date, source/terms URLs,
  notice hashes and restricted redistribution metadata into the output directory
  (default `../dvsi-vectors`). Reports retain this evidence and the terms notice.
  See [reference-data handling](../legal/reference-data.md) for corpus notices,
  listening-pack checks and recording user-acquired ETSI inputs.
- `tv-rc/*.pcm` are the inputs (s16le, 8 kHz). Each mode directory holds
  `X.bit` (DVSI encoder) and `X.pcm` (DVSI decoder). The fetcher takes the modes
  this library decodes: `dstar`, `p25`, `p25_nofec`, `r33` and `r34`.
  `cmprc.txt` records the rate words: D-STAR uses
  `0x0130 0x0763 0x4000 0x0000 0x0000 0x0048`.
- `tv-rc/<from>/<to>/` (only `alltone`, `dam` and `sine0_4k`) hold DVSI's rate
  conversions between those modes, made with `-rc -rd <from> -re <to>`: DVSI's
  decoded model of the `<from>` bits, re-quantized into `<to>`. They are lossy
  probes of DVSI's model, not an exact readout of it.
- `--from-dvsi` takes 9-byte D-STAR DV data (LSB first), 18-byte P25 frames, or
  9-byte rate-33 frames (MSB first, DMR interleave). It writes Framef rows for
  `mbe_quality_eval`. Bits are only deinterleaved, so the `_eN` channel-error
  variants still exercise FEC and repeat handling. A truncated final frame is
  rejected.
- The `_eN_hd` and `_eN_sd` variants carry the same channel errors, at bit
  error rates of about 1.1, 2.2, 4.9 and 10.6% for `e1`, `e2`, `e5` and `e10`.
  `_sd.bit` (DVSI's `-sdbits 4`) holds one 4-bit value per bit of the `_hd.bit`
  file, two per byte, high nibble first, in that file's MSB-first bit order;
  values 8 to 15 decide 1, so the `_hd` file is exactly their hard decisions.
  `--from-dvsi` recognizes the `_sd.bit` name and writes one lowercase hex digit
  per frame cell; `mbe_quality_eval --soft` decodes those rows through the
  `*SoftFramef` APIs, with each value's distance from 7.5 as the reliability
  (16 to 240).
- Copy the inputs and DVSI's decoded files to `.raw` names before using them as
  `--ref` or `--decoded` operands. The evaluator selects the format by
  extension.

### DVSI scoreboard

`tools/quality/dvsi_scoreboard.py` runs every speech vector of the D-STAR, P25
and rate-33 sets through the tools above and writes `scoreboard.json` and
`scoreboard.md` (needs NumPy). A change is judged by scoring a baseline build
and the candidate, each aligning freely, then both again at their joint lags,
the candidate against the baseline:

```sh
SB="python3 tools/quality/dvsi_scoreboard.py --vectors ../dvsi-vectors --decode-tools build/tools-debug"
Q=build/quality/scoreboard
$SB --tools build/baseline --out $Q/baseline-free
$SB --tools build/tools-debug --out $Q/candidate-free
JOINT=$Q/baseline-free/scoreboard.json,$Q/candidate-free/scoreboard.json
$SB --tools build/baseline --out $Q/baseline --lags-from $JOINT
$SB --tools build/tools-debug --out $Q/candidate --lags-from $JOINT \
  --compare $Q/baseline/scoreboard.json --gates pitch-law
```

- **Speech:** ours vs the input, DVSI's decoder vs the input, and ours vs
  DVSI's decoder, with the `mbe_quality_eval` metrics above (absolute level,
  RMS-matched spectra and, against DVSI's decoder, `waveform_corr`).
- **Pitch:** a harmonic-peak f0 estimator measures our output, DVSI's output
  and the input on stable voiced frames and reports each against the pitch our
  decoder coded, by f0 band. A decoder that plays sharp or flat shows up as a
  ratio away from 1.
- **Encoder:** our encoding of the input (D-STAR, IMBE and AMBE+2 encoders)
  against DVSI's encoding of it, both decoded by the `--decode-tools` build
  (default `--tools`): pitch error
  against the input, octave errors, per-1 kHz-band voicing agreement, the speech
  metrics against the input, the same metrics of our decoded encoding against
  DVSI's (`encoder.vs_dvsi_bits.*`, with `encoder.level_gap_db` the unsigned
  level difference), rail samples beyond DVSI's encoding and the fraction of
  tone or silence frames. DVSI's encoding is the reference there (its active
  frames and level); DVSI's encoders delay the speech about 200 samples more
  than ours (125 for P25), past the evaluator's -160 lag limit, so ours is
  first delayed by 320 samples (240 for P25). It cannot show how DVSI hardware
  decodes our bits. Each encoding's rows are hashed (`rows_sha256`); with
  `--compare` against a baseline that recorded them, `encoder.rows_changed_vectors`
  counts the vectors whose encoder bits changed. A build whose
  `mbe_quality_encode` does not offer a codec skips that encoder, so this
  script can score a baseline built from an older checkout. Give the baseline
  and the candidate run the same `--decode-tools`, so that only the encoders
  differ (where a change moves both, decode with the candidate's decoder: it
  stands for the decoders that will receive the bits). The run records that
  decoder's identity (`decoder_identity`, a SHA-256 of its `mbe_quality_eval`
  and the library beside it), and `--compare` refuses encoder deltas between
  runs decoded by different or unrecorded decoders.
- **Encoder tones (rate 33 and D-STAR):** our encoding of each tone vector
  against DVSI's bits, frame by frame at the best frame shift (of equally good
  shifts, the one where DVSI's encoder runs a frame behind ours, as on every
  speech vector): the share of
  DVSI's tone frames we send with the same tone index
  (`tone.encoder_agreement`), the share of frames we send as a tone where DVSI
  sends voice (`tone.encoder_extra_rate`) and the mean level difference where
  both agree, in AD steps for rate 33 (`tone.encoder_level_error_ad`) and in
  steps of D-STAR's 8-bit volume (`tone.encoder_level_error_volume`).
  The tones group reports the worst vector of each and, for these, also the
  mean over vectors (`_mean`), since DVSI itself is inconsistent on its
  deliberately malformed tone vectors (`cpvbad`, `dtmfvbad`).
- Every comparison must find an alignment inside the evaluator's lag search.
  A decoder comparison without one (a silent decode, or a lag at a limit)
  fails the run. An encoding without one (one that is mostly silence frames,
  say) is left unscored, counted in `encoder.alignment_failures`, and fails
  any gate on its partition.
- **Channel errors** (`--partitions channel`): `dam`'s `_eN_hd` and `_eN_sd`
  vectors, decoded as hard and as soft decisions, against DVSI's decode of the
  same bits. A loud frame (RMS 300 or more in DVSI's decode of the clean bits)
  is muted by DVSI when its decode of the errored bits is 30 dB or more below
  that; ours is muted when the frame carries `MBE_PROCESS_FLAG_MUTE`. The
  section reports, per vector, how often our mute decisions agree with DVSI's
  on loud frames, the false and missed mute rates, the first loud frame each
  mutes, our repeat rate, and the log-spectral distance from DVSI's output on
  the loud frames with no muted frame at or next to them (a frame's spectra
  reach into its neighbours), with our output aligned and level-matched by the
  clean bits' comparison. Each run mutes its own frames, so `--compare`
  measures that distance for both runs on the frames that count in both
  (`channel.hd.lsd_common_db`) and reports their share of all the loud frames
  DVSI plays (`coverage`), so frames a run did not measure count as lost. Soft
  vectors need a baseline build whose tools read them (`--from-dvsi` of an
  `_sd.bit` file and `mbe_quality_eval --soft`): for an older library, build it
  with this checkout's `tools/quality` sources.
- **Partitions** are by source recording. Several vectors are one recording
  through different filters, noise or levels (found by matching 3 s blocks of
  the inputs at any lag): `dam` and its variants (`lvl` among them),
  `fambm22a`, `tambf22a`, `ucarm15a`, `ucarf15a` and `ucarm20a`, `tambf22b` and
  `tambf32b`, and `tia11` and `p01mirs`. Each such group sits in one partition,
  and the bootstrap intervals resample whole groups. `development` and
  `validation` were both examined during earlier changes; `final` (`irstia`
  and `mark`) is scored once per change set (`--partitions final`). Not scored
  as speech: `xfer`, the DTX and channel-error variants and the synthetic
  pulse trains and sweeps. The `tones` partition covers the DTMF, single,
  KNOX, alert and call-progress vectors.
- **Tones** are compared with one acoustic detector run on both outputs, on
  windows where DVSI's output is steady: level, frequencies, balance, how
  often both agree on the tone playing, and the fraction of windows playing a
  different tone. Single tones must agree within 10 Hz (they are 31.25 Hz
  apart); each dual-tone component within 25 Hz, since DVSI plays dual tones
  up to 20 Hz off nominal and neighbouring DTMF frequencies are 73 Hz or more
  apart. Since 2.6.13 this library plays the frequencies DVSI's decoders play
  for a rate-33 or D-STAR tone, and plays the P25 tones it recognises as DVSI's
  decoder does, so a tone can also be compared with DVSI's decode sample for
  sample, frame for frame at lag 0: the correlation of each steady tone frame
  (a tone frame between two tone frames of the same tone, DVSI's output above
  a level) and the residual of a tone's first frame and of the frame after a
  tone. The scoreboard does not compute those; the figures in
  [DVSI parity](dvsi-parity.md#november-2026-tones-and-the-octave-step) were
  measured so with a local script, whose frame selection is stated there. Its
  P25 tone frames are those its own implementation of the library's tone test
  passes (both limits, 0.375 and 0.0235). That implementation is not a copy of
  the library's classifier (it lacks the positive-normal test of the tone
  magnitudes and adds 1e-12 to a denominator), but its selection agrees with
  the library's test on every one of the 26,200 frames of the P25 tone
  vectors. A selection made with another classifier counts other frames.
- `--compare` adds paired per-vector deltas with bootstrap 95% intervals and the
  worst vector, and counts the baseline's speech vectors whose PCM changed
  (`speech.pcm_changed_vectors`). Each tone vector's and each channel-error
  vector's decode of DVSI's bits is hashed by its samples (`pcm_sha256`), and
  each mode's tones and channel groups count the baseline's vectors whose
  decode changed (`tone.pcm_changed_vectors`, `channel.pcm_changed_vectors`).
  In all three counts a vector missing from the run, or without a hash in
  either run, counts as changed, so a baseline scored by a version of the
  script that recorded no hashes cannot pass an identity gate. A tone or
  channel gate set's `speech.pcm_changed_vectors` check speaks only for the
  speech vectors the run scored, counted against the baseline: a run whose
  two reports leave out the same speech vectors passes it. That the speech
  partitions are complete comes from evaluating the speech gate sets in the
  same run, as the protocol does. Both runs are
  measured at one lag per comparison, so a change that moves the evaluator's alignment
  search does not also move the samples measured: `--lags-from` takes the lags
  a scoreboard recorded (`lags` in each vector's entry, with the frame shift of
  the encoder-tone comparison), or, given two, their joint lags. Where the two
  free runs' searches disagree, the joint lag is that of the run whose
  alignment reached the higher correlation (`alignment`), since a change can
  let a decoder align where the old one could not: main's D-STAR gain law
  played onsets up to 16 dB low and pulled its alignment with DVSI's output to
  40-96 samples, where the corrected decoder aligns at -16 to 0 (P25 and rate
  33 at -16 to +16 in both builds). Without `--lags-from`, `--compare` uses
  the baseline's lags. A vector
  whose recorded lags still differ from the baseline's is not paired.
  `--gates` evaluates gate sets (a comma list) from
  `tools/quality/dvsi_gates.json`, declared before the change it judges, and
  exits non-zero when a gate fails. A gate set fails unless the run was scored
  at fixed lags and both runs scored every vector of its partition in every
  mode it gates, and a paired check fails when a metric was measured for a
  vector in only one of the two runs.
  `all_no_increase` and `all_no_decrease` hold every vector to the limit, not
  just the mean, and a check's `min_coverage` fails it when its common frames
  are too few. A metric that neither run reports for a vector (as
  `waveform_corr_unshifted` at a lag past 32 samples) leaves that vector out
  of a check; `every_vector` fails the check instead.
  A paired check needs its metric from both runs on every vector it covers,
  so check the baseline's coverage of a metric before declaring a gate on it:
  the P25 checks of `tone-model-tones` (the gate sets of the tone model change are
  `tone-model-development`, `-validation`, `-tones`, `-channel` and `-final`)
  fail on `knox_9`, whose P25 decode by the baseline has no steady tone window
  the detector recognises, while the candidate's has
  ([DVSI parity](dvsi-parity.md#november-2026-tones-and-the-octave-step)).
  The gate sets of the mid-frame voicing change (`mid-voicing-development`,
  `-validation`, `-tones`, `-channel` and `-final`) were declared after that
  check: the description of `mid-voicing-tones` records, from the baseline's
  own scoreboard, which tone metrics each vector can be measured on (`cp1`,
  `cpvbad` and `dtmfvbad` have no steady tone window there, and `alert`, a
  single tone with steady windows, has no balance; the description wrongly
  puts `alert` with the three)
  ([DVSI parity](dvsi-parity.md#november-2026-mid-frame-voicing)).
  The gate sets of rate 33's mid-frame pitch wave (`r33-pitch-wave-development`,
  `-validation`, `-tones`, `-channel` and `-final`, declared in `8942924`)
  record their coverage the same way, `alert` with its single tone among the
  vectors with no balance only, and check the tone measures in rate 33 alone,
  D-STAR's and P25's tone decodes by their hashes
  ([DVSI parity](dvsi-parity.md#november-2026-rate-33s-mid-frame-pitch-wave)).
  The gate sets of the mid-frame columns change (`mid-columns-development`,
  `-validation`, `-tones`, `-channel` and `-final`, declared in `3593e8c`)
  require every vector in each paired check whose metric the baseline's own
  scoreboard measured on every vector of the partition, record the coverage of
  the others in their descriptions, and gate no channel-error vector on waveform
  correlation (P25's falls with the change on the four hard-decision vectors and
  rises on the four soft-decision ones, while a spectral distance on the loud
  changed frames falls with hard decisions; the cause of the fall was not
  established). A check on a metric that applies to part of a partition cannot
  require every vector, so a vector that lost the metric in both runs would
  pass it unnoticed; the results section gives this run's coverage
  ([DVSI parity](dvsi-parity.md#november-2026-mid-frame-voicing-by-columns)).
  The gate sets of the row-17 train change (`r33-row17-train-development`, `-validation`,
  `-tones`, `-channel` and `-final`, declared in `6a32fb8`) gate rate 33's speech, tone and
  channel measures, D-STAR's and P25's decodes by their hashes, and every encoder's rows; every
  paired check of the development, validation and final sets requires every vector, and the
  tones and channel sets record the coverage of the others; where their descriptions say more
  than the checks (the channel distance's frame selection, the tones' encoder level population,
  the encoder-row counter's skipping of a vector the baseline has no hash for, and what
  validation's description says of a gain) is in the results section's notes
  ([DVSI parity](dvsi-parity.md#november-2026-the-pulse-train-through-a-row-17-frame)). An
  identity counter that skips a vector the baseline has no hash for (`encoder_row_changes()`)
  needs its coverage checked with every scoring. Mutations of the change checked to fail the
  tests: dropping the new condition fails `test_waveform_through_changes` (rate-33 rendering),
  `test_train_through_row17`, `test_row17_mark_in_codec_states` and
  `test_marks_act_in_their_own_rendering`; honouring the mark in every rendering fails
  `test_row17_mark_in_codec_states` (in D-STAR's path and in P25's alone),
  `test_marks_act_in_their_own_rendering` and `test_foreign_marks_in_history`; returning from
  `mbe_advance_subframes()` at once after a rate-33 previous model marked 6, making row 15 start
  early in the rate-33 rendering, and loading the phase of a previous model marked 6 without its
  pulses each fail `test_train_through_row17` run alone (the last also
  `test_waveform_through_changes`).
  The gate sets of the pulsed state (`pulsed-state-development`, `-validation`, `-tones`,
  `-channel` and `-final`, declared in `525d017`; `-final` scored once, on 2026-11-27, after the
  two reviews and their fixes: passed, 18 of 18) gate D-STAR's and P25's speech, tone and channel
  measures (a rise of the waveform correlation and a fall of the distance from DVSI's decode on the
  development, validation and final partitions, the distance from the input allowed to rise within
  a limit, the change's measured cost), rate 33's decodes by their hashes and every encoder's rows;
  every paired check of the development, validation and final sets requires every vector, and the
  tones and channel sets record the coverage of the others. Where their descriptions say more than
  the checks or the measurements (the rule's exceptions for the mid-frame model, where the tone
  decodes change, the alternatives they attribute to scratch library builds, what the channel and
  tone selections exclude, and that they predate P25's validity counter) is in the results
  section's notes ([DVSI parity](dvsi-parity.md#november-2026-the-pulsed-state)).

`mbe_quality_eval --params frames.jsonl` writes one JSON record per frame:
`bits` (the frame's parameter bits after FEC), `decoded` (the frame's own
model from the public `mbe_decode*Parms` decoder,
with its status), `history` (the prediction state after the frame) and `synth`
(the enhanced model synthesized). Each holds `w0` (radians per sample), `L` and
per-band `Vl`, `Ml` (linear amplitude as synthesized; for AMBE unvoiced bands
it includes the unvoiced scale factor) and `log2Ml` (the codec's log2 spectral
amplitude before that factor), plus `gamma`. `synth` also holds `PHIl` and
`PSIl`, the model's harmonic phases for the frame. Synthesis renders each voiced
harmonic with `PHIl` plus the output high-pass's phase at its frequency
(`mbe_low_band_phase()`: 180 degrees toward DC, 90 at 100 Hz, under 9 above
1 kHz), so phase readbacks can be checked against the decoder itself. Dumping
does not change the PCM.

### Encoder evaluation

`mbe_quality_encode` (built with `MBELIB_BUILD_TOOLS=ON`) encodes 8 kHz s16le
speech with the D-STAR (`--codec ambe2400`), AMBE+2 (`ambe2450`) or IMBE
(`imbe7200`) encoder and writes the 49- or 88-bit rows `mbe_quality_eval`
decodes, so encoder output is scored against the input with the same metrics
as decoder output. Comparing with DVSI's own encoding of the same input:

```sh
build/dev-debug/mbe_quality_encode --codec ambe2400 --in dam.raw --out dam.rows
build/dev-debug/mbe_quality_eval --codec ambe2400 --frames dam.rows --out ours.wav --ref dam.raw --json ours.json
build/dev-debug/mbe_quality_reframe --codec ambe2400 --from-dvsi ../dvsi-vectors/tv-rc/dstar/dam.bit --out dvsi.frames
build/dev-debug/mbe_quality_eval --codec ambe2400 --frames dvsi.frames --out dvsi.wav --ref dam.raw --json dvsi.json
```

Both bit streams then go through the same decoder, which itself matches DVSI's
output on DVSI's bits, so the band, LSD and envelope figures compare the
encoders. One zero flush frame (`--flush-frames`, default 1) covers the
encoders' analysis delay.

### Attribution and limitations

Shared phase regeneration and corrected WOLA were already present at checkpoint
`8da251d`; interpolation, enhancement, adaptive smoothing and most excitation
state handling predate this quality-completion work. The additional production
corrections are limited to clamping the **applied** adaptive attenuation to
nonnegative values while preserving signed threshold recurrence, and dispatching
accepted AMBE2400 tone IDs 5/6 consistently with IDs 7…122 while setting `TONE`.
Tone mapping, NOTONES behavior, limiter, codec concealment policy, public headers,
state layout and exports are unchanged. Regressions cover positive attenuation,
fade/recovery, tones and rejection thresholds, warm replay/interleaving, all
160 WOLA reconstruction samples, all 256 IMBE fundamentals, and corrected ECC.

Historical three-mode figures under the old tooling were absolute **signed-mean**
crest errors 1.225/1.941/2.086 → 0.301/0.205/0.317 dB, not mean per-clip absolute
errors. That 27-second corpus included the duplicated `hts1a`; historic tooling
hashes do not bind every report to the current evaluator/calibrator. Preserve
those artifacts as historical input, not regenerated four-mode evidence.
Spectral/envelope improvements alone never establish listener preference.
Without the specified actual human judgments, perceptual acceptance remains
`not_established`; an objective failure is a completed negative assessment,
not grounds for claiming “all modes improved.”

The completed fixed-profile assessment is preserved under
`build/quality/runs/{original-to-final-debug,checkpoint-to-final-debug,heldout-debug,heldout-minus12db-debug,heldout-release,heldout-minus12db-release}/`.
Each comparison contains its own immutable manifest and generated acceptance,
analysis, tables, figures, and listening material. `build/quality/evidence-index.json`
indexes the delivered runs, controls, compatibility proof, and benchmark results.

On the 40 held-out references (118.032 seconds), Debug and matched Release agree
on the following aggregate decisions; each level uses the same frozen gains:

| Actual frame mode | Correctness | Native objective | −12 dB objective | Perceptual |
| --- | --- | --- | --- | --- |
| IMBE7200 | pass | fail | fail | not established |
| IMBE7100 | pass | fail | fail | not established; shared IMBE PCM |
| AMBE2450 | pass | fail | fail | not established |
| AMBE2400 | pass | pass | pass | not established |

Debug IMBE mean per-clip absolute crest error increases 1.072→1.145 dB at native
level and 1.007→1.023 dB at −12 dB. AMBE2450 mean join excess increases
0.332→0.564 dB and 0.043→0.303 dB respectively, exceeding its +0.1 dB allowance.
AMBE2400 passes the aggregate screens but retains 17 native and 15 attenuated
per-clip exceptions. Mean STOI/ESTOI decrease slightly in all three independent
speech modes; they do not establish perceptual harm or preference by themselves.
There are no candidate rail-count increases in these held-out runs.
The two additional production fixes preserve all legacy clean checkpoint PCM;
the held-out screen failures concern the pre-existing branch synthesis changes,
not evidence for retuning either correction. The result does **not** establish
“all modes improved.”

After the automated assessment, the maintainer completed all 48 native-level
blinded presentations as an informal single-listener preference session:
14 baseline preferences, 14 candidate preferences, and 20 ties. These are
presentation counts, not 48 independent pairs: the two IMBE modes repeat the
same PCM comparisons. Only the first two items received full numerical ratings;
the remaining scores were left blank. The strong second-played preference
(A: 2, B: 26, ties: 20) is a playback-order concern, not proof of its cause.
The formal multi-listener perceptual result remains `not_established`.
The maintainer explicitly chose to retain the candidate despite the mixed
objective and personal listening results; that decision does not alter the
recorded measurement gates.
