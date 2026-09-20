#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Exercise the examples through nonseekable binary stdin/stdout pipes."""

import io
import struct
import subprocess
import sys
import wave


def chunk(tag, payload):
    return tag + struct.pack("<I", len(payload)) + payload + b"\0" * (len(payload) & 1)


def riff(*chunks):
    body = b"WAVE" + b"".join(chunks)
    return b"RIFF" + struct.pack("<I", len(body)) + body


def pcm_format(channels=1):
    return struct.pack("<HHIIHH", 1, channels, 8000, 16000 * channels, 2 * channels, 16)


def invoke(command, data, *args):
    result = subprocess.run(
        [*command, *args], input=data, capture_output=True,
        timeout=60,
    )
    return result


def success(command, data):
    result = invoke(command, data)
    assert result.returncode == 0, result.stderr
    return result.stdout


def check_wav(data, frames):
    assert data[:4] == b"RIFF" and data[8:12] == b"WAVE"
    assert struct.unpack_from("<I", data, 4)[0] == len(data) - 8
    assert struct.unpack_from("<I", data, 40)[0] == frames * 320
    assert len(data) == 44 + frames * 320
    with wave.open(io.BytesIO(data), "rb") as wav:
        assert (wav.getnchannels(), wav.getsampwidth(), wav.getframerate()) == (1, 2, 8000)
        assert wav.getnframes() == frames * 160


def verify(encode, decode):
    frames = 12
    pcm = b"".join(struct.pack("<h", ((i * 997) % 24001) - 12000) for i in range(frames * 160))
    fmt = chunk(b"fmt ", pcm_format())
    data = chunk(b"data", pcm)
    tail = chunk(b"LIST", b"INFO" + chunk(b"INAM", b"pipe fixture\0"))
    wav = riff(fmt, data, tail)
    encoded = success(encode, wav)
    assert len(encoded) == 12 * (frames + 1)
    assert all(encoded[i:i + 3] == b"\x55\x2d\x16" for i in range(0, len(encoded), 12))
    check_wav(success(decode, encoded), frames + 1)
    print(f"PASS pipes + trailing LIST: {frames} input frames, {len(encoded)} encoded bytes, "
          f"{160 * (frames + 1)} decoded samples")

    # Chunk extensions and odd unknown-chunk padding must never enter the PCM.
    extended = riff(chunk(b"JUNK", b"\r\n\x1a"), chunk(b"fmt ", pcm_format() + b"x"), data, tail)
    assert success(encode, extended) == encoded
    partial = invoke(encode, riff(fmt, chunk(b"data", pcm + b"\x01\x00"), tail))
    assert partial.returncode == 0 and partial.stdout == encoded
    assert b"1 trailing samples ignored" in partial.stderr
    assert success(encode, riff(fmt, data)) == encoded
    print("PASS chunk padding, extended fmt, data_size boundary, partial-frame warning")

    invalid = {
        "odd data": riff(fmt, chunk(b"data", pcm + b"x")),
        "stereo": riff(chunk(b"fmt ", pcm_format(2)), data),
        "data before fmt": riff(data, fmt),
        "short fmt": riff(chunk(b"fmt ", pcm_format()[:15]), data),
        "wrong rate": riff(chunk(b"fmt ", struct.pack("<HHIIHH", 1, 1, 16000, 32000, 2, 16)), data),
        "non PCM": riff(chunk(b"fmt ", struct.pack("<HHIIHH", 3, 1, 8000, 16000, 2, 16)), data),
        "wrong bits": riff(chunk(b"fmt ", struct.pack("<HHIIHH", 1, 1, 8000, 8000, 1, 8)), data),
        "chunk exceeds RIFF": riff(fmt, b"data" + struct.pack("<I", 0xFFFFFFFE)),
        "padding exceeds RIFF": riff(b"JUNK" + struct.pack("<I", 0xFFFFFFFF)),
        "truncated header": wav[:11],
        "truncated fmt": wav[:30],
        "truncated PCM": wav[:44 + len(pcm) - 1],
        "truncated LIST": wav[:-1],
        "trailing chunk exceeds RIFF": riff(fmt, data, b"LIST" + struct.pack("<I", 100)),
        "partial trailing chunk header": riff(fmt, data, b"LIST"),
    }
    for label, malformed in invalid.items():
        result = invoke(encode, malformed)
        assert result.returncode == 1, (label, result.returncode, result.stderr)
        assert result.stderr, label
        assert len(result.stdout) % 12 == 0 and not result.stdout.startswith(b"RIFF"), label
    for malformed in (encoded[:-1], encoded + b"x", b"\x55"):
        result = invoke(decode, malformed)
        assert result.returncode == 1 and result.stdout == b"" and result.stderr
    for command in (encode, decode):
        for args in (("unused",), ("input", "output")):
            result = invoke(command, b"", *args)
            assert result.returncode == 1 and result.stdout == b"" and b"usage:" in result.stderr
    print(f"PASS {len(invalid)} malformed WAVs, truncated containers, and path arguments: exit 1")

    silence = success(encode, riff(fmt, chunk(b"data", bytes(len(pcm)))))
    # Fixed b0=127/tone-index=128 payload, independent of encoder startup.
    silence_frame = bytes.fromhex("55 2d 16 e3 0b 2c c3 42 18 20 40 11")
    assert silence[-12:] == silence_frame
    silent_wav = success(decode, silence_frame * (frames + 1))
    check_wav(silent_wav, frames + 1)
    assert not any(silent_wav[44:])
    burst = bytes(len(pcm) - 64) + b"".join(
        struct.pack("<h", 24000 if i % 20 < 10 else -24000) for i in range(32)
    )
    flushed = success(encode, riff(fmt, chunk(b"data", burst)))
    assert len(flushed) == len(encoded) and flushed[-12:] != silence[-12:]
    flushed_wav = success(decode, flushed)
    check_wav(flushed_wav, frames + 1)
    assert any(flushed_wav[-320:])
    print("PASS exact silence mute and last-32-sample burst: non-silent flush frame")

    mismatched = invoke(decode, b"\x1a\r\n" + encoded[3:])
    assert mismatched.returncode == 0
    check_wav(mismatched.stdout, frames + 1)
    assert b"1 of 13 frames had sync word mismatches" in mismatched.stderr
    check_wav(success(decode, b""), 0)
    print("PASS sync-word warning and empty-container WAV")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: test_dstar_examples.py ENCODER DECODER")
    verify([sys.argv[1]], [sys.argv[2]])
