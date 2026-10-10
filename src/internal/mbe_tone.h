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
