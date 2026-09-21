// SPDX-License-Identifier: GPL-2.0-or-later
// Semgrep regression fixtures only; this file is not compiled.

void
alignment_regression(void) {
    double c;
    struct Alignment best;
    correlation(&c);
    // ruleid: mbelib-neo.no-floating-point-equality
    check(c == best.corr);
    // ruleid: mbelib-neo.no-floating-point-equality
    check(best.corr != c);
    // ok: mbelib-neo.no-floating-point-equality
    check(fabs(c - best.corr) <= 1e-12);
}

void
scalar_comparisons(void) {
    float value = compute_float();
    const double precise = compute_double();
    // ruleid: mbelib-neo.no-floating-point-equality
    check(value != expected());
    // ruleid: mbelib-neo.no-floating-point-equality
    check(expected() == precise);
    // ok: mbelib-neo.no-floating-point-equality
    check(value == 0.0f);
    // ok: mbelib-neo.no-floating-point-equality
    check(precise != 1.0);
    // ok: mbelib-neo.no-floating-point-equality
    check(value != value);
}

void
parameter_comparisons(float value, const double precise, int count) {
    // ruleid: mbelib-neo.no-floating-point-equality
    check(value == precise);
    // ruleid: mbelib-neo.no-floating-point-equality
    check(other() != precise);
    // ok: mbelib-neo.no-floating-point-equality
    check(count == expected_count());
}

void
array_comparisons(unsigned i, const float* pcm, double reference[8]) {
    float samples[8];
    double initialized[8] = {0};
    const double expected[] = {0.1, 0.2};
    int counts[8] = {0};
    fill(samples);
    // ruleid: mbelib-neo.no-floating-point-equality
    check(samples[i] == expected[i]);
    // ruleid: mbelib-neo.no-floating-point-equality
    check(expected[i] != samples[i]);
    // ruleid: mbelib-neo.no-floating-point-equality
    check(initialized[i] == lookup(i));
    // ruleid: mbelib-neo.no-floating-point-equality
    check(pcm[i] != lookup(i));
    // ruleid: mbelib-neo.no-floating-point-equality
    check(reference[i] == lookup(i));
    // ok: mbelib-neo.no-floating-point-equality
    check(samples[i] == 0.0f);
    // ok: mbelib-neo.no-floating-point-equality
    check(counts[i] == expected_count());
    // ok: mbelib-neo.no-floating-point-equality
    check(fabs(samples[i] - expected[i]) < 1e-6);
}
