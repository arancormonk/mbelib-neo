// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_TONE_H
#define MBELIB_NEO_INTERNAL_MBE_TONE_H

struct mbe_tone_frequency {
    float freq1;
    float freq2;
};

static inline int
mbe_tone_lookup_freqs(int tone_id, float* freq1, float* freq2) {
    static const struct mbe_tone_frequency dual_tones[36] = {
        {1336.0f, 941.0f}, {1209.0f, 697.0f}, {1336.0f, 697.0f}, {1477.0f, 697.0f}, {1209.0f, 770.0f},
        {1336.0f, 770.0f}, {1477.0f, 770.0f}, {1209.0f, 852.0f}, {1336.0f, 852.0f}, {1477.0f, 852.0f},
        {1633.0f, 697.0f}, {1633.0f, 770.0f}, {1633.0f, 852.0f}, {1633.0f, 941.0f}, {1209.0f, 941.0f},
        {1477.0f, 941.0f}, {1162.0f, 820.0f}, {1052.0f, 606.0f}, {1162.0f, 606.0f}, {1279.0f, 606.0f},
        {1052.0f, 672.0f}, {1162.0f, 672.0f}, {1279.0f, 672.0f}, {1052.0f, 743.0f}, {1162.0f, 743.0f},
        {1279.0f, 743.0f}, {1430.0f, 606.0f}, {1430.0f, 672.0f}, {1430.0f, 743.0f}, {1430.0f, 820.0f},
        {1052.0f, 820.0f}, {1279.0f, 820.0f}, {440.0f, 350.0f},  {480.0f, 440.0f},  {620.0f, 480.0f},
        {490.0f, 350.0f},
    };

    *freq1 = 0.0f;
    *freq2 = 0.0f;

    if (tone_id == 5) {
        *freq1 = 156.25f;
        *freq2 = *freq1;
        return 1;
    }
    if (tone_id == 6) {
        *freq1 = 187.5f;
        *freq2 = *freq1;
        return 1;
    }
    if ((tone_id >= 7) && (tone_id <= 122)) {
        *freq1 = 31.25f * (float)tone_id;
        *freq2 = *freq1;
        return 1;
    }
    if ((tone_id >= 128) && (tone_id <= 163)) {
        const struct mbe_tone_frequency tone = dual_tones[tone_id - 128];
        *freq1 = tone.freq1;
        *freq2 = tone.freq2;
        return 1;
    }

    return 0;
}

/*
 * D-STAR tone index to frequencies. Single tones share AMBE+2's indices.
 * DTMF is 128 + 4 * column + row (rows 697..941 Hz, columns 1209..1633 Hz),
 * and 144..147 are AMBE+2's call-progress tones 160..163; DVSI's D-STAR
 * encoder sends KNOX tones as voice. Established from the tones DVSI's
 * AMBE-3000 plays for its D-STAR tone vectors.
 */
static inline int
mbe_tone_lookup_dstar_freqs(int tone_id, float* freq1, float* freq2) {
    static const float rows[4] = {697.0f, 770.0f, 852.0f, 941.0f};
    static const float columns[4] = {1209.0f, 1336.0f, 1477.0f, 1633.0f};
    if ((tone_id >= 128) && (tone_id <= 143)) {
        *freq1 = rows[(tone_id - 128) & 3];
        *freq2 = columns[(tone_id - 128) >> 2];
        return 1;
    }
    if ((tone_id >= 144) && (tone_id <= 147)) {
        return mbe_tone_lookup_freqs(tone_id + 16, freq1, freq2);
    }
    return (tone_id <= 122) && mbe_tone_lookup_freqs(tone_id, freq1, freq2);
}

/*
 * The frequencies DVSI's decoders play for a tone frame, which are not Table
 * 9's nominal ones: DVSI's AMBE-3000 plays a tone through its voiced
 * synthesizer as a model of one or two harmonics of a fundamental, the tone
 * frame parameters of TIA-102.BABA-1 Annex J (f0, l1, l2 per tone index; clause
 * 8 turns a tone frame into such a model). Measured on DVSI's rate-33 and
 * D-STAR decodes of its tone vectors (both codecs play the same values), the
 * fundamentals lie on a grid of 2^-19 cycle per sample to 0.003 of a step.
 * Below, a component's frequency is counted in those steps: Annex J's f0 on
 * that grid, round(f0 * 2^19 / 8000), for 34 of the 36 dual tones; IDs 145 and
 * 147 are DVSI's own (9883 = 150.80 Hz and 5667 = 86.47 Hz, where the text
 * prints 150.89 and 86.5). So DTMF "1" (ID 129) plays 693.909 + 1214.340 Hz,
 * where Table 9 has 697 + 1209.
 */
#define MBE_TONE_GRID_BITS 19

struct mbe_tone_harmonics {
    unsigned short f0; /* fundamental, in 2^-19 cycle per sample */
    unsigned char l1;  /* the lower component's harmonic number */
    unsigned char l2;  /* the upper component's */
};

/*
 * The played components of an AMBE+2 tone ID, in 2^-19 cycle per sample,
 * `upper` and `lower` (equal for a single tone); 0 for an ID with no tone.
 * Dual tones (128..163): l2 f0 and l1 f0 of the table. Single tones (5..122):
 * harmonic l = ID / 16 + 1 of the fundamental floor(2048 ID / l), so that the
 * fundamental stays under 500 Hz (DVSI's decodes: ID 30 plays 937.5000 Hz, ID
 * 57 1781.2500 and ID 46 1437.470, three times 31402 steps; back-to-back
 * single tones keep the pulse train's phase only for these harmonic numbers),
 * within 0.1 Hz of the nominal 31.25 ID.
 */
static inline int
mbe_tone_played_steps(int tone_id, unsigned* upper, unsigned* lower) {
    static const struct mbe_tone_harmonics dual_tones[36] = {
        {5145, 12, 17}, {11369, 4, 7},  {4588, 10, 19}, {5702, 8, 17}, {7206, 7, 11},  {12562, 4, 7},
        {4599, 11, 21}, {4657, 12, 17}, {7968, 7, 11},  {13894, 4, 7}, {7629, 6, 14},  {6301, 8, 17},
        {4653, 12, 23}, {15352, 4, 7},  {8807, 7, 9},   {8805, 7, 11}, {4478, 12, 17}, {9883, 4, 7},
        {4445, 9, 17},  {5667, 7, 15},  {6278, 7, 11},  {10939, 4, 7}, {4437, 10, 19}, {4898, 10, 14},
        {6940, 7, 11},  {6080, 8, 14},  {6655, 6, 14},  {5506, 8, 17}, {4445, 11, 21}, {6704, 8, 14},
        {7668, 7, 9},   {7700, 7, 11},  {5753, 4, 5},   {4642, 6, 7},  {7995, 4, 5},   {4588, 5, 7},
    };

    *upper = 0u;
    *lower = 0u;
    if ((tone_id >= 5) && (tone_id <= 122)) {
        const unsigned l = (unsigned)tone_id / 16u + 1u;
        *upper = l * ((2048u * (unsigned)tone_id) / l);
        *lower = *upper;
        return 1;
    }
    if ((tone_id >= 128) && (tone_id <= 163)) {
        const struct mbe_tone_harmonics tone = dual_tones[tone_id - 128];
        *upper = (unsigned)tone.l2 * tone.f0;
        *lower = (unsigned)tone.l1 * tone.f0;
        return 1;
    }
    return 0;
}

/*
 * The played components of a D-STAR tone index, in the order
 * mbe_tone_lookup_dstar_freqs() gives the nominal ones (DTMF: the row tone,
 * then the column tone; call progress and single tones as AMBE+2's), through
 * the same index mapping. 0 for an index with no tone.
 */
static inline int
mbe_tone_dstar_played_steps(int tone_id, unsigned* first, unsigned* second) {
    /* D-STAR DTMF 128 + 4 * column + row as the AMBE+2 ID of the same pair. */
    static const unsigned char dtmf[16] = {129, 132, 135, 142, 130, 133, 136, 128,
                                           131, 134, 137, 143, 138, 139, 140, 141};
    unsigned upper, lower;
    if ((tone_id >= 128) && (tone_id <= 143)) {
        (void)mbe_tone_played_steps(dtmf[tone_id - 128], &upper, &lower);
        *first = lower;
        *second = upper;
        return 1;
    }
    *first = 0u;
    *second = 0u;
    if ((tone_id >= 144) && (tone_id <= 147)) {
        return mbe_tone_played_steps(tone_id + 16, first, second);
    }
    return (tone_id <= 122) && mbe_tone_played_steps(tone_id, first, second);
}

/* D-STAR 8-bit tone volume: bits 12..16, 44, 45, 17 from most to least significant. */
static inline int
mbe_tone_dstar_volume_bit(int i) {
    static const unsigned char order[8] = {12, 13, 14, 15, 16, 44, 45, 17};
    return order[i & 7];
}

static inline int
mbe_tone_dstar_volume(const char* ambe_d) {
    int volume = 0;
    for (int i = 0; i < 8; i++) {
        volume = (volume << 1) | (ambe_d[mbe_tone_dstar_volume_bit(i)] & 1);
    }
    return volume;
}

static inline void
mbe_tone_dstar_set_volume(char* ambe_d, int volume) {
    for (int i = 0; i < 8; i++) {
        ambe_d[mbe_tone_dstar_volume_bit(i)] = (char)((volume >> (7 - i)) & 1);
    }
}

/*
 * Tone component levels, in dB relative to an RMS of 32768 on the int16
 * scale. AMBE+2 follows TIA-102.BABA-1 7.2: +3.17 dBm0 at AD 127 down to
 * -87.13 dBm0 at AD 0, 0.711 dB per step; DVSI's decoder plays +3.17 dBm0 at
 * -2.52 dB per component. D-STAR's 8-bit volume steps by half as much; its law
 * is fitted to DVSI's D-STAR tone vectors (within 0.04 dB from volume 181 to
 * 207). Dual tones carry both components at the level.
 */
#define MBE_TONE_AMBE2450_DB_AT_AD127 (-2.523)
#define MBE_TONE_AMBE2450_DB_PER_STEP ((3.17 + 87.13) / 127.0)
#define MBE_TONE_DSTAR_DB_AT_VOL255   3.540
#define MBE_TONE_DSTAR_DB_PER_STEP    0.35495

static inline double
mbe_tone_ambe2450_level_db(int ad) {
    return MBE_TONE_AMBE2450_DB_AT_AD127 + (MBE_TONE_AMBE2450_DB_PER_STEP * (double)(ad - 127));
}

static inline double
mbe_tone_dstar_level_db(int volume) {
    return MBE_TONE_DSTAR_DB_AT_VOL255 + (MBE_TONE_DSTAR_DB_PER_STEP * (double)(volume - 255));
}

static inline int
mbe_tone_id_is_valid(int tone_id) {
    float freq1, freq2;
    return mbe_tone_lookup_freqs(tone_id, &freq1, &freq2);
}

#endif /* MBELIB_NEO_INTERNAL_MBE_TONE_H */
