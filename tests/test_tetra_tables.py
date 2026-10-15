#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check the TETRA table generators and the committed tables they produced.

The codebook extractor runs on a synthetic archive (no ETSI data); the committed codebook
file is checked for its shape and provenance; the computed constants are regenerated from
their formulas and compared with the committed header.
"""

import contextlib
import hashlib
import importlib.util
import io
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
EXTRACTOR = ROOT / "tools" / "tetra" / "extract_acelp_tables.py"
CONSTANTS = ROOT / "tools" / "tetra" / "gen_acelp_constants.py"
CODEBOOKS_C = ROOT / "src" / "tetra" / "tetra_acelp_codebooks.c"
CODEBOOKS_H = ROOT / "src" / "internal" / "tetra_acelp_codebooks.h"
SHAPES = {"tetra_lsp_cb1": (256, 3), "tetra_lsp_cb2": (512, 3), "tetra_lsp_cb3": (512, 4), "tetra_gain_cb": (64, 2)}


def load(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    assert spec is not None and spec.loader is not None, path
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic_archive(path, gain_count=128):
    def array(name, values):
        body = ",\n".join(values)
        return f"/* {name} */\nstatic Word16 {name}[] = {{\n{body}\n}};\n"

    lsp = array("dico1_clsp", [str(i - 300) for i in range(768)])
    lsp += array("dico2_clsp", [f"{i} * 2" for i in range(1536)])
    lsp += "/* a comment with dico3_clsp[] = { 9 }; inside */\n"
    lsp += array("dico3_clsp", [str(-i) for i in range(2048)])
    gain = array("t_qua_ener", [str(i * 3) for i in range(gain_count)])
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("C-CODE/CLSP_334.TAB", lsp)
        archive.writestr("C-CODE/ENER_QUA.TAB", gain)
        archive.writestr("C-CODE/OTHER.C", "not read")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expect_error(extractor, function, *args):
    try:
        function(*args)
    except extractor.ExtractError:
        return
    raise AssertionError(f"{function.__name__} did not fail")


def test_extractor(root):
    extractor = load(EXTRACTOR)
    archive = root / "archive.zip"
    digest = synthetic_archive(archive)
    known = {digest: "synthetic.zip (test)"}

    members, description = extractor.read_members(archive, known)
    assert description == "synthetic.zip (test)"
    assert set(members) == {"C-CODE/CLSP_334.TAB", "C-CODE/ENER_QUA.TAB"}
    text = extractor.render(members, known)
    assert text.startswith("// SPDX-License-Identifier: LicenseRef-ETSI-EN-300-395-2-data\n")
    assert '#include "tetra_acelp_codebooks.h"' in text
    assert f"synthetic.zip (test): {digest}" in text
    assert "const int16_t tetra_lsp_cb1[256][3] = {\n    {  -300,   -299,   -298}," in text
    assert "const int16_t tetra_lsp_cb2[512][3] = {\n    {     0,      2,      4}," in text
    assert "    {   465,    466,    467},\n};" in text  # last row of cb1
    assert "const int16_t tetra_gain_cb[64][2] = {\n    {     0,      3}," in text
    assert "static" not in text

    # Unknown archives, unexpected elements, wrong counts and 16-bit overflow are refused.
    expect_error(extractor, extractor.read_members, archive, {"0" * 64: "other"})
    bad = root / "bad.zip"
    digest_bad = synthetic_archive(bad, gain_count=127)
    members_bad, _ = extractor.read_members(bad, {digest_bad: "bad"})
    expect_error(extractor, extractor.render, members_bad, {digest_bad: "bad"})
    expect_error(extractor, extractor.parse_array, "static Word16 a[] = { 40000 };", "a", 1)
    expect_error(extractor, extractor.parse_array, "static Word16 a[] = { 1 + 2 };", "a", 1)
    expect_error(extractor, extractor.parse_array, "static Word16 b[] = { 1 };", "a", 1)
    assert extractor.parse_array("static Word16 a[2] = { -5, 16383 * 2 };", "a", 2) == [-5, 32766]

    # Write, then --check: a match exits 0, an edited file exits 1.
    setattr(extractor, "KNOWN_ARCHIVES", known)
    out = root / "codebooks.c"
    with contextlib.redirect_stdout(io.StringIO()):
        assert extractor.main([str(archive), "--output", str(out)]) == 0
        assert extractor.main([str(archive), "--output", str(out), "--check"]) == 0
    out.write_text(out.read_text(encoding="ascii").replace("-300", "-301", 1), encoding="ascii")
    with contextlib.redirect_stderr(io.StringIO()):
        assert extractor.main([str(archive), "--output", str(out), "--check"]) == 1
        assert extractor.main([str(root / "missing.zip"), "--output", str(out)]) == 2


def test_committed_codebooks():
    extractor = load(EXTRACTOR)
    text = CODEBOOKS_C.read_text(encoding="ascii")
    assert text.startswith("// SPDX-License-Identifier: LicenseRef-ETSI-EN-300-395-2-data\n")
    for digest, name in extractor.KNOWN_ARCHIVES.items():
        assert f" *   {name}: {digest}\n" in text, name
    header = CODEBOOKS_H.read_text(encoding="ascii")
    for name, (rows, columns) in SHAPES.items():
        assert f"extern const int16_t {name}[{rows}][{columns}];" in header, name
        match = re.search(rf"const int16_t {name}\[{rows}\]\[{columns}\] = {{\n(.*?)\n}};", text, re.S)
        assert match, name
        lines = match.group(1).split("\n")
        assert len(lines) == rows, (name, len(lines))
        for line in lines:
            values = [int(v) for v in line.strip().strip("{},").split(",")]
            assert len(values) == columns and all(-32768 <= v <= 32767 for v in values), line
    # The LSP sub-vectors are ordered (cosines decrease) within each codebook entry.
    for name in ("tetra_lsp_cb1", "tetra_lsp_cb2", "tetra_lsp_cb3"):
        rows, _ = SHAPES[name]
        match = re.search(rf"{name}\[{rows}\]\[\d\] = {{\n(.*?)\n}};", text, re.S)
        assert match, name
        body = match.group(1)
        for line in body.split("\n"):
            values = [int(v) for v in line.strip().strip("{},").split(",")]
            assert all(a > b for a, b in zip(values, values[1:])), (name, line)


def test_constants():
    result = subprocess.run([sys.executable, "-I", str(CONSTANTS), "--check"], capture_output=True, timeout=60)
    assert result.returncode == 0, result.stderr


def main():
    with tempfile.TemporaryDirectory() as tmp:
        test_extractor(Path(tmp))
    test_committed_codebooks()
    test_constants()
    return 0


if __name__ == "__main__":
    sys.exit(main())
