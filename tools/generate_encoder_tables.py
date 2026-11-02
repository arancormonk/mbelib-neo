#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check or regenerate encoder constants from the offline C reference.

Only the Python standard library is needed. See docs/encoder-tables.md.
"""

import argparse
from collections import defaultdict
from pathlib import Path
import re
import struct
import subprocess


FLOAT_TOKEN = re.compile(r"(?<![\w.])[-+]?(?:[0-9]+\.[0-9]*|\.[0-9]+|[0-9]+)(?:[eE][-+]?[0-9]+)?f\b")
COMMENTS = re.compile(r"/\*.*?\*/|//[^\n]*", re.DOTALL)
FLOAT_TABLES = {
    "analysis": ("mbe_analysis_tables.h", 1443),
    "ambe": ("ambe_dct_lut.h", 1849),
    "pitch": ("mbe_pitch_model_lut.h", 248),
    "hann": ("mbe_tone_hann.h", 160),
    "imbe": ("imbe_enc_dct_lut.h", 385),
}


def float_bits(token):
    return struct.unpack("<I", struct.pack("<f", float(token[:-1])))[0]


def float_literal(bits):
    value = struct.unpack("<f", struct.pack("<I", bits))[0]
    result = format(value, ".9g")
    if "." not in result and "e" not in result:
        result += ".0"
    return result + "f"


def round_ratio_even(numerator, denominator):
    """Round an exact rational to nearest even, including negative values."""
    quotient, remainder = divmod(abs(numerator), denominator)
    if 2 * remainder > denominator or (2 * remainder == denominator and quotient % 2):
        quotient += 1
    return -quotient if numerator < 0 else quotient


def replace_tokens(text, matches, values):
    # Work backwards so replacing a literal cannot move the remaining spans.
    for match, value in reversed(list(zip(matches, values))):
        text = text[:match.start()] + value + text[match.end():]
    return text


def check_float_table(path, expected, write):
    text = path.read_text(encoding="utf-8")
    # Mask comments while retaining offsets into the original header.
    code = COMMENTS.sub(lambda match: " " * len(match.group()), text)
    matches = list(FLOAT_TOKEN.finditer(code))
    if len(matches) != len(expected):
        raise ValueError(f"{path.name}: expected {len(expected)} float literals, found {len(matches)}")
    actual = [float_bits(match.group()) for match in matches]
    count = sum(a != b for a, b in zip(actual, expected))
    if count and write:
        text = replace_tokens(text, matches, [float_literal(value) for value in expected])
        path.write_text(text, encoding="utf-8")
    return count


def check_integer_table(path, name, expected, write):
    text = path.read_text(encoding="utf-8")
    body = re.search(r"\b" + re.escape(name) + r"\s*\[[^]]+\]\s*=\s*\{([^}]+)\}", text)
    if body is None:
        raise ValueError(f"{path.name}: array {name} was not found")
    code = COMMENTS.sub(lambda match: " " * len(match.group()), text)
    matches = list(re.compile(r"[0-9]+").finditer(code, body.start(1), body.end(1)))
    if len(matches) != len(expected):
        raise ValueError(f"{name}: expected {len(expected)} integers, found {len(matches)}")
    actual = [int(match.group()) for match in matches]
    count = sum(a != b for a, b in zip(actual, expected))
    if count and write:
        text = replace_tokens(text, matches, [str(value) for value in expected])
        path.write_text(text, encoding="utf-8")
    return count


def reference_values(executable):
    output = subprocess.run([str(executable.resolve())], check=True, capture_output=True, text=True)
    groups = defaultdict(list)
    for line in output.stdout.splitlines():
        group, value = line.split()
        groups[group].append(int(value, 10 if group == "harmonics" else 16))
    expected_counts = {name: count for name, (_, count) in FLOAT_TABLES.items()}
    expected_counts["harmonics"] = 208
    if {name: len(values) for name, values in groups.items()} != expected_counts:
        raise ValueError("The reference executable returned unexpected table names or sizes")
    return groups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path, help="compiled generate_encoder_tables.c executable")
    parser.add_argument("--write", action="store_true", help="replace differing literals (default: check only)")
    args = parser.parse_args()
    headers = Path(__file__).resolve().parents[1] / "src" / "internal"
    groups = reference_values(args.reference)
    differences = 0
    for name, (filename, _) in FLOAT_TABLES.items():
        count = check_float_table(headers / filename, groups[name], args.write)
        differences += count
        print(f"{filename}: {count} differing float values")

    offsets = [0, 0]
    for length in range(1, 10):
        offsets.append(offsets[-1] + length * length)
    indices = [110 + round_ratio_even(n * 110, half)
               for half in range(84, 111) for n in range(-half, half + 1)]
    integer_tables = (
        ("mbe_pitch_model_lut.h", "imbe_pitch_harmonics", groups["harmonics"]),
        ("imbe_enc_dct_lut.h", "imbe_enc_dct_offset", offsets),
        ("mbe_analysis_window_index.h", "mbe_window_indices", indices),
    )
    for filename, name, values in integer_tables:
        count = check_integer_table(headers / filename, name, values, args.write)
        differences += count
        print(f"{name}: {count} differing integer values")
    if differences:
        print(f"{differences} values {'updated' if args.write else 'differ; files unchanged'}")
    else:
        print("All encoder table values match the reference.")
    return 0 if args.write or not differences else 1


if __name__ == "__main__":
    raise SystemExit(main())
