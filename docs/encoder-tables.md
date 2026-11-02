# Fixed encoder coefficients

The private lookup tables replace repeated evaluations of fixed expressions. They preserve
the single-precision values used by the reference build, including small cosine residuals
near zero. The AMBE DCT table also supplies the existing shared decoder caches; their DCT
arithmetic is unchanged.

| Header in `src/internal/` | Source expression or domain |
| --- | --- |
| `mbe_analysis_tables.h` | Kaiser windows, low-pass taps, and the reconstruction-window transform |
| `ambe_dct_lut.h` | Existing AMBE DCT cosine expressions, in cache initialization order |
| `mbe_pitch_model_lut.h` | AMBE pitch laws for indices 0–123; IMBE harmonic counts for 0–207 |
| `mbe_tone_hann.h` | Existing 160-sample tone observation window |
| `imbe_enc_dct_lut.h` | Existing float-angle DCT expressions for lengths 1–10 and packed offsets |
| `mbe_analysis_window_index.h` | Nearest-even window indices for half-lengths 84–110 |

The window-index table applies only to round-to-nearest operation. Other floating-point
rounding modes use the existing runtime expression. The offline generator uses exact
integer ratios for this finite domain.

## Check the values

`tools/generate_encoder_tables.c` contains the original reference formulas. It is a
standalone development tool, not part of the library build. It intentionally retains the
original intermediate precision so it can check the stored results. The Python wrapper
requires only Python 3's standard library and compares every stored float by its binary32
bits, rather than by an approximate numeric tolerance. It also checks the integer tables.

From the repository root, build the reference program in a scratch directory. In an MSVC
Developer Command Prompt:

```bat
mkdir table-check
cl /nologo /O2 /fp:precise /MD tools/generate_encoder_tables.c /Fe:table-check/reference.exe /Fo:table-check/reference.obj
python tools/generate_encoder_tables.py --reference table-check/reference.exe
```

For a C99 compiler on a Unix-like system:

```sh
mkdir -p table-check
cc -std=c99 -O2 -ffp-contract=off tools/generate_encoder_tables.c -lm -o table-check/reference
python3 tools/generate_encoder_tables.py --reference table-check/reference
```

The checked-in values were produced with MSVC 19.29 on x64 under `/fp:precise`. Other
compilers or math libraries can produce different final bits for transcendental functions.
A failed check on another toolchain identifies that difference; it does not imply the
checked-in table is corrupt. Do not enable fast-math or change the arithmetic ordering
when checking the reference.

## Regenerate after an intentional formula change

Update the reference formula, affected table dimensions, and wrapper counts together, then pass `--write`
to the same Python command. This replaces differing numeric literals in the headers while
preserving their declarations and comments. The default command is read-only. Format the
changed headers with the repository's `clang-format` configuration and rerun the check.
Treat changed coefficient bits as a codec behavior change and validate the encoder output
before accepting them. Routine builds do not need either generator.
