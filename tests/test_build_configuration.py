#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Configure-only checks of per-configuration build settings.

Configures the source tree with Ninja and Ninja Multi-Config, using the
parent build's PFFFT sources, and reads the generated compile and link
commands:

  (a) Ninja Multi-Config, no build type, LTO on: Release compiles with LTO
      and SIMD; Debug compiles with SIMD but without LTO.
  (b) Release: the libraries (and the tests that select their kernels) get
      the native flag when the configure log reports one; Debug never does.
  (c) a target ISA in CMAKE_C_FLAGS (-march=x86-64 and the like): no
      native flag.
  (d) MBELIB_ENABLE_NATIVE=OFF: no native flag.

Where the configure log reports the run-time AVX2 dispatch, the AVX2 kernel
copies compile with its flags and MBELIB_AVX2_KERNELS in untuned
configurations (Debug, MBELIB_ENABLE_NATIVE=OFF) and without them, empty, in
tuned ones.

Exits 77 (skip) when Ninja is not available.
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

NATIVE_RE = re.compile(r"MBELIB_ENABLE_NATIVE: release-like configurations use (\S+)")
DISPATCH_RE = re.compile(r"AVX2 dispatch: (.+?) kernels chosen at run time")
AVX2_UNIT = "src/core/mbe_voiced_avx2.c"
LTO_FLAGS = ("-flto", "/GL")


def configure(args, build, generator, *options):
    command = [
        args.cmake, "-S", args.source, "-B", str(build), "-G", generator,
        f"-DCMAKE_MAKE_PROGRAM={args.ninja}",
        f"-DCMAKE_C_COMPILER={args.c_compiler}",
        f"-DFETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT={args.pffft_source}",
        "-DFETCHCONTENT_FULLY_DISCONNECTED=ON",
        "-DMBELIB_BUILD_EXAMPLES=OFF",
        *options,
    ]
    result = subprocess.run(command, capture_output=True, text=True, timeout=600)
    log = result.stdout + result.stderr
    if result.returncode != 0:
        sys.exit(f"configure failed: {' '.join(command)}\n{log}")
    return log


def commands(args, build, target, ninja_file=None):
    command = [args.ninja, "-C", str(build)]
    if ninja_file:
        command += ["-f", ninja_file]
    command += ["-t", "commands", target]
    return subprocess.run(command, capture_output=True, text=True, check=True, timeout=120).stdout.splitlines()


def compile_line(lines, source):
    """The compile command of one source file (forward or back slashes)."""
    pattern = re.compile(r"[\\/]" + re.escape(source).replace("/", r"[\\/]") + r"(\s|$)")
    matches = [line for line in lines if pattern.search(line) and (" -c " in line or " /c " in line)]
    if not matches:
        sys.exit(f"no compile command for {source}")
    return matches[0]


def check(condition, message):
    if not condition:
        sys.exit(f"FAIL: {message}")
    print(f"PASS: {message}")


def has_lto(line):
    return any(token.startswith(LTO_FLAGS) for token in line.split())


def scenario_multi_config(args, work):
    build = work / "multi"
    log = configure(args, build, "Ninja Multi-Config", "-DMBELIB_ENABLE_LTO=ON", "-DMBELIB_BUILD_TESTS=OFF")
    release = compile_line(commands(args, build, "mbe-shared", "build-Release.ninja"), "src/core/mbelib.c")
    debug = compile_line(commands(args, build, "mbe-shared", "build-Debug.ninja"), "src/core/mbelib.c")
    check("MBELIB_ENABLE_SIMD=1" in release, "(a) multi-config Release compiles with SIMD")
    check("MBELIB_ENABLE_SIMD=1" in debug, "(a) multi-config Debug compiles with SIMD")
    if "IPO/LTO not supported" in log:
        print("note: the compiler does not support IPO; LTO checks skipped")
    else:
        check(has_lto(release), "(a) multi-config Release compiles with LTO")
    check(not has_lto(debug), "(a) multi-config Debug compiles without LTO")
    tuned = native_flag_of(log) is not None
    check_dispatch(args, build, log, not tuned, "(a) multi-config Release", "build-Release.ninja")
    check_dispatch(args, build, log, True, "(a) multi-config Debug", "build-Debug.ninja")


def dispatch_flags_of(log):
    match = DISPATCH_RE.search(log)
    return match.group(1).split() if match else None


def check_dispatch(args, build, log, expect_copies, name, ninja_file=None):
    flags = dispatch_flags_of(log)
    if flags is None:
        print(f"note: no run-time AVX2 dispatch for this target; {name} dispatch checks skipped")
        return
    line = compile_line(commands(args, build, "mbe-shared", ninja_file), AVX2_UNIT)
    tokens = line.split()
    present = all(flag in tokens for flag in flags) and any("MBELIB_AVX2_KERNELS=1" in t for t in tokens)
    absent = not any(flag in tokens for flag in flags) and not any("MBELIB_AVX2_KERNELS" in t for t in tokens)
    if expect_copies:
        check(present, f"{name}: the AVX2 kernel copies compile with {' '.join(flags)}")
    else:
        check(absent, f"{name}: the AVX2 kernel copies compile empty")


def native_flag_of(log):
    match = NATIVE_RE.search(log)
    return match.group(1) if match else None


def library_line(args, build):
    return compile_line(commands(args, build, "mbe-shared"), "src/core/mbelib.c")


def scenario_release(args, work):
    build = work / "release"
    log = configure(args, build, "Ninja", "-DCMAKE_BUILD_TYPE=Release", "-DMBELIB_BUILD_TESTS=ON")
    flag = native_flag_of(log)
    if flag is None:
        reason = re.search(r"MBELIB_ENABLE_NATIVE: not tuning for the build host: (.*)", log)
        check(reason is not None, "(b) Release reports why it does not tune for the host")
        print(f"note: no native tuning on this host ({reason.group(1)})")
        return None
    check(flag in library_line(args, build).split(), f"(b) Release library compiles with {flag}")
    check_dispatch(args, build, log, False, "(b) tuned Release")
    golden = compile_line(commands(args, build, "test_golden_pcm"), "tests/test_golden_pcm.c")
    check(flag in golden.split(), f"(b) Release test_golden_pcm compiles with {flag}")
    api = compile_line(commands(args, build, "test_api"), "tests/test_api.c")
    check(flag not in api.split(), f"(b) Release test_api compiles without {flag}")
    debug = work / "debug"
    debug_log = configure(args, debug, "Ninja", "-DCMAKE_BUILD_TYPE=Debug", "-DMBELIB_BUILD_TESTS=OFF",
                          "-DMBELIB_ENABLE_SIMD=ON")
    check(flag not in library_line(args, debug).split(), f"(b) Debug library compiles without {flag}")
    check_dispatch(args, debug, debug_log, True, "(b) Debug")
    return flag


def scenario_no_native(args, work, flag, name, *options):
    build = work / name
    log = configure(args, build, "Ninja", "-DCMAKE_BUILD_TYPE=Release", "-DMBELIB_BUILD_TESTS=OFF", *options)
    check(native_flag_of(log) is None, f"{name}: configure reports no native flag")
    if flag:
        check(flag not in library_line(args, build).split(), f"{name}: library compiles without {flag}")
    return build, log


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cmake", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--pffft-source", required=True)
    parser.add_argument("--c-compiler", required=True)
    parser.add_argument("--user-isa", default="", help="a target ISA option the compiler accepts")
    parser.add_argument("--ninja", default="")
    args = parser.parse_args()
    args.ninja = args.ninja or shutil.which("ninja") or shutil.which("ninja-build") or ""
    if not args.ninja:
        print("ninja not found; skipping")
        return 77
    with tempfile.TemporaryDirectory(prefix="mbelib-build-configuration-") as directory:
        work = Path(directory)
        scenario_multi_config(args, work)
        flag = scenario_release(args, work)
        if args.user_isa:
            scenario_no_native(args, work, flag, "(c) user target ISA", f"-DCMAKE_C_FLAGS={args.user_isa}")
        else:
            print("note: no known target ISA option for this compiler; (c) skipped")
        build, log = scenario_no_native(args, work, flag, "(d) MBELIB_ENABLE_NATIVE=OFF",
                                        "-DMBELIB_ENABLE_NATIVE=OFF")
        check_dispatch(args, build, log, True, "(d) MBELIB_ENABLE_NATIVE=OFF")
    return 0


if __name__ == "__main__":
    sys.exit(main())
