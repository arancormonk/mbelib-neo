#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Checks of the per-configuration build settings and the host-tuning rules.

The rules of cmake/MbeNativeArch.cmake run first, in CMake script mode: when
host tuning must not apply (cross-compiling, a compiler target, an Apple
device SDK, foreign or several OSX architectures, a target ISA in the compiler
command, the C flags or MSVC's CL/_CL_ environment), and the MinGW/Cygwin GCC
assembler requirement.

Then the source tree is configured, configure-only, with Ninja and Ninja
Multi-Config, using the parent build's PFFFT sources, and the generated
compile commands are read:

  (a) Ninja Multi-Config, LTO on, with or without CMAKE_BUILD_TYPE: Release
      compiles with LTO and (64-bit targets) SIMD; Debug compiles with SIMD
      but without LTO.
  (b) Release: the libraries (and the tests that select their kernels) get
      the native flag when the configure log reports one; Debug never does.
  (c) a target ISA in CMAKE_C_FLAGS or the compiler command: no native flag.
  (d) MBELIB_ENABLE_NATIVE=OFF: no native flag.
  (e) Ninja Multi-Config with only Debug targeting AVX2: Release still builds
      the AVX2 kernel copies.

Where the configure log reports the run-time AVX2 dispatch, the AVX2 kernel
copies compile with its flags and MBELIB_AVX2_KERNELS in untuned
configurations and without them, empty, in tuned ones.

Build settings a developer may export (CFLAGS, CL, CMAKE_GENERATOR, ...) are
removed from the environment of every configure. The configure-only scenarios
exit 77 (skip) without Ninja, and are skipped when the forwarded compiler
command already selects a target ISA.
"""

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

NATIVE_RE = re.compile(r"MBELIB_ENABLE_NATIVE: release-like configurations use (.+?) \(")
DISPATCH_RE = re.compile(r"AVX2 dispatch: (.+?) kernels chosen at run time")
AVX2_UNIT = "src/core/mbe_voiced_avx2.c"
LTO_FLAGS = ("-flto", "/GL")
USER_ISA = re.compile(r"(?:^|\s)[\"']?(?:-march=|-mcpu=|[-/]arch:|-arch\s)", re.I)
ISOLATED_ENV = ("CFLAGS", "CXXFLAGS", "LDFLAGS", "CL", "_CL_", "CMAKE_TOOLCHAIN_FILE", "CMAKE_OSX_ARCHITECTURES",
                "CMAKE_BUILD_TYPE", "CMAKE_CONFIGURATION_TYPES", "CMAKE_GENERATOR")
ENV = {name: value for name, value in os.environ.items() if name not in ISOLATED_ENV}

RULES_SCRIPT = r'''
include("${MBELIB_SOURCE_DIR}/cmake/MbeNativeArch.cmake")
mbelib_native_arch_skip_reason(reason)
message("REASON=[${reason}]")
if(DEFINED COMPOSE_FLAG)
    mbelib_native_arch_compose(flags why "${COMPOSE_FLAG}")
    message("FLAGS=[${flags}]")
    message("WHY=[${why}]")
endif()
'''


def check(condition, message):
    if not condition:
        sys.exit(f"FAIL: {message}")
    print(f"PASS: {message}")


def run_rules(args, environment=None, **variables):
    with tempfile.TemporaryDirectory(prefix="mbelib-native-rules-") as directory:
        script = Path(directory) / "rules.cmake"
        script.write_text(RULES_SCRIPT, encoding="utf-8")
        command = [args.cmake, f"-DMBELIB_SOURCE_DIR={Path(args.source).as_posix()}"]
        command += [f"-D{name}={value}" for name, value in variables.items()]
        command += ["-P", str(script)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=120,
                                env={**ENV, **(environment or {})})
        log = result.stdout + result.stderr
        if result.returncode != 0:
            sys.exit(f"rules script failed: {' '.join(command)}\n{log}")
        return {key: value for key, value in re.findall(r"^(REASON|FLAGS|WHY)=\[(.*)\]$", log, re.M)}


def check_rules(args):
    cases = [
        ({}, r"^$"),
        ({"CMAKE_CROSSCOMPILING": "TRUE"}, r"^cross-compiling$"),
        ({"CMAKE_SYSTEM_NAME": "Android"}, r"^Android target$"),
        ({"CMAKE_SYSTEM_NAME": "iOS"}, r"^iOS target$"),
        ({"CMAKE_SYSTEM_NAME": "Emscripten"}, r"^Emscripten target$"),
        ({"CMAKE_C_COMPILER_TARGET": "aarch64-linux-gnu"}, r"^CMAKE_C_COMPILER_TARGET selects aarch64-linux-gnu$"),
        ({"CMAKE_SYSTEM_NAME": "Darwin", "CMAKE_HOST_SYSTEM_PROCESSOR": "arm64",
          "CMAKE_OSX_ARCHITECTURES": "x86_64;arm64"}, r"several architectures"),
        ({"CMAKE_SYSTEM_NAME": "Darwin", "CMAKE_HOST_SYSTEM_PROCESSOR": "arm64",
          "CMAKE_OSX_ARCHITECTURES": "x86_64"}, r"is not the host architecture \(arm64\)$"),
        ({"CMAKE_SYSTEM_NAME": "Darwin", "CMAKE_HOST_SYSTEM_PROCESSOR": "x86_64",
          "CMAKE_OSX_ARCHITECTURES": "x86_64"}, r"^$"),
        # CMake ignores CMAKE_OSX_ARCHITECTURES off Apple platforms.
        ({"CMAKE_SYSTEM_NAME": "Linux", "CMAKE_OSX_ARCHITECTURES": "x86_64;arm64"}, r"^$"),
        ({"CMAKE_SYSTEM_NAME": "Darwin", "CMAKE_HOST_SYSTEM_PROCESSOR": "x86_64",
          "CMAKE_OSX_SYSROOT": "iphonesimulator"}, r"^CMAKE_OSX_SYSROOT selects iphonesimulator$"),
        ({"CMAKE_SYSTEM_NAME": "Darwin", "CMAKE_OSX_SYSROOT":
          "/Applications/Xcode.app/Contents/Developer/Platforms/AppleTVOS.platform/Developer/SDKs/AppleTVOS.sdk"},
         r"^CMAKE_OSX_SYSROOT selects "),
        ({"CMAKE_SYSTEM_NAME": "Darwin", "CMAKE_OSX_SYSROOT": "macosx"}, r"^$"),
        ({"CMAKE_C_FLAGS": "-O2 -march=x86-64-v2"}, r"^CMAKE_C_FLAGS already selects a target ISA$"),
        ({"CMAKE_C_FLAGS": "-pipe -mcpu=cortex-a72"}, r"^CMAKE_C_FLAGS already selects"),
        ({"CMAKE_C_FLAGS": "-O2\t-march=x86-64"}, r"^CMAKE_C_FLAGS already selects"),
        ({"CMAKE_C_FLAGS": '-O2 "-march=x86-64"'}, r"^CMAKE_C_FLAGS already selects"),
        ({"CMAKE_C_FLAGS": "-O2 '-mcpu=cortex-a72'"}, r"^CMAKE_C_FLAGS already selects"),
        ({"CMAKE_C_FLAGS_RELEASE": "-O2 -arch arm64"}, r"^CMAKE_C_FLAGS_RELEASE already"),
        ({"CMAKE_C_FLAGS_RELWITHDEBINFO": "-O2 -g -march=znver3"}, r"^CMAKE_C_FLAGS_RELWITHDEBINFO already"),
        ({"CMAKE_C_FLAGS_MINSIZEREL": "/O1 /arch:AVX2"}, r"^CMAKE_C_FLAGS_MINSIZEREL already"),
        ({"CMAKE_C_FLAGS_MINSIZEREL": "/O1 -ARCH:AVX2"}, r"^CMAKE_C_FLAGS_MINSIZEREL already"),
        ({"CMAKE_C_COMPILER_ARG1": " -march=x86-64"}, r"^CMAKE_C_COMPILER_ARG1 already selects"),
        ({"CMAKE_C_COMPILER_ARG1": " gcc"}, r"^$"),
        # Only release-like configurations are tuned; neither -mtune nor an
        # option merely containing the text selects a target ISA.
        ({"CMAKE_C_FLAGS_DEBUG": "-march=x86-64"}, r"^$"),
        ({"CMAKE_C_FLAGS": "-O2 -mtune=generic"}, r"^$"),
        ({"CMAKE_C_FLAGS": "-I/opt/x-march=native"}, r"^$"),
    ]
    for variables, expected in cases:
        reason = run_rules(args, **variables)["REASON"]
        check(re.search(expected, reason) is not None, f"rule {variables}: {reason!r}")
    env_cases = [
        ({"CL": "/arch:AVX"}, {"MSVC": "TRUE"}, r"^the CL environment variable already selects"),
        ({"_CL_": "/O2 /arch:AVX2"}, {"CMAKE_C_COMPILER_FRONTEND_VARIANT": "MSVC"},
         r"^the _CL_ environment variable already selects"),
        ({"CL": "/arch:AVX"}, {}, r"^$"),
    ]
    for environment, variables, expected in env_cases:
        reason = run_rules(args, environment, **variables)["REASON"]
        check(re.search(expected, reason) is not None, f"rule env {environment} {variables}: {reason!r}")
    compose_cases = [
        ({"CMAKE_C_COMPILER_ID": "GNU"}, r"^-march=native$"),
        ({"WIN32": "1", "CMAKE_C_COMPILER_ID": "Clang"}, r"^-march=native$"),
        ({"WIN32": "1", "CMAKE_C_COMPILER_ID": "GNU", "MBELIB_C_HAS_UNALIGNED_VECTOR_MOVE": "1"},
         r"^-march=native;-Wa,-muse-unaligned-vector-move$"),
        ({"CYGWIN": "1", "CMAKE_C_COMPILER_ID": "GNU", "MBELIB_C_HAS_UNALIGNED_VECTOR_MOVE": "1"},
         r"^-march=native;-Wa,-muse-unaligned-vector-move$"),
        ({"WIN32": "1", "CMAKE_C_COMPILER_ID": "GNU", "MBELIB_C_HAS_UNALIGNED_VECTOR_MOVE": "0"}, r"^$"),
    ]
    for variables, expected in compose_cases:
        result = run_rules(args, COMPOSE_FLAG="-march=native", **variables)
        check(re.search(expected, result["FLAGS"]) is not None, f"compose {variables}: {result}")
        if expected == r"^$":
            check("-muse-unaligned-vector-move" in result["WHY"], f"compose {variables} says why: {result['WHY']!r}")


def configure(args, build, generator, *options, source=None):
    command = [
        args.cmake, "-S", source or args.source, "-B", str(build), "-G", generator,
        f"-DCMAKE_MAKE_PROGRAM={args.ninja}",
        f"-DCMAKE_C_COMPILER={args.c_compiler}",
        f"-DFETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT={args.pffft_source}",
        "-DFETCHCONTENT_FULLY_DISCONNECTED=ON",
        "-DMBELIB_BUILD_EXAMPLES=OFF",
    ]
    # CC can be a launcher ("ccache gcc") or carry arguments ("gcc -pipe");
    # CMake keeps the rest of that command in ARG1.
    if args.c_compiler_arg1.strip():
        command.append(f"-DCMAKE_C_COMPILER_ARG1={args.c_compiler_arg1}")
    command += options
    result = subprocess.run(command, capture_output=True, text=True, timeout=600, env=ENV)
    log = result.stdout + result.stderr
    if result.returncode != 0:
        sys.exit(f"configure failed: {' '.join(command)}\n{log}")
    return log


def commands(args, build, target, ninja_file=None):
    command = [args.ninja, "-C", str(build)]
    if ninja_file:
        command += ["-f", ninja_file]
    command += ["-t", "commands", target]
    return subprocess.run(command, capture_output=True, text=True, check=True, timeout=120,
                          env=ENV).stdout.splitlines()


def compile_line(lines, source):
    """The compile command of one source file (forward or back slashes, quoted or not)."""
    pattern = re.compile(r"[\\/]" + re.escape(source).replace("/", r"[\\/]") + r"[\"']?(\s|$)")
    matches = [line for line in lines if pattern.search(line) and (" -c " in line or " /c " in line)]
    if not matches:
        sys.exit(f"no compile command for {source}")
    return matches[0]


def tokens_of(line):
    try:
        return shlex.split(line)
    except ValueError:
        return line.split()


def has_lto(line):
    return any(token.startswith(LTO_FLAGS) for token in tokens_of(line))


def pointer_size(build):
    for path in (build / "CMakeFiles").glob("*/CMakeCCompiler.cmake"):
        match = re.search(r'set\(CMAKE_C_SIZEOF_DATA_PTR "(\d+)"\)', path.read_text(encoding="utf-8"))
        if match:
            return int(match.group(1))
    sys.exit(f"no CMAKE_C_SIZEOF_DATA_PTR in {build}")


def dispatch_flags_of(log):
    match = DISPATCH_RE.search(log)
    return match.group(1).split() if match else None


def check_dispatch(args, build, log, expect_copies, name, ninja_file=None):
    flags = dispatch_flags_of(log)
    if flags is None:
        print(f"note: no run-time AVX2 dispatch for this target; {name} dispatch checks skipped")
        return
    tokens = tokens_of(compile_line(commands(args, build, "mbe-shared", ninja_file), AVX2_UNIT))
    present = all(flag in tokens for flag in flags) and any("MBELIB_AVX2_KERNELS=1" in t for t in tokens)
    absent = not any(flag in tokens for flag in flags) and not any("MBELIB_AVX2_KERNELS" in t for t in tokens)
    if expect_copies:
        check(present, f"{name}: the AVX2 kernel copies compile with {' '.join(flags)}")
    else:
        check(absent, f"{name}: the AVX2 kernel copies compile empty")


def native_flags_of(log):
    match = NATIVE_RE.search(log)
    return match.group(1).split() if match else None


def library_tokens(args, build, ninja_file=None):
    return tokens_of(compile_line(commands(args, build, "mbe-shared", ninja_file), "src/core/mbelib.c"))


def scenario_multi_config(args, work, name, *options, source=None):
    build = work / name
    log = configure(args, build, "Ninja Multi-Config", "-DMBELIB_ENABLE_LTO=ON", "-DMBELIB_BUILD_TESTS=OFF",
                    *options, source=source)
    release = compile_line(commands(args, build, "mbe-shared", "build-Release.ninja"), "src/core/mbelib.c")
    debug = compile_line(commands(args, build, "mbe-shared", "build-Debug.ninja"), "src/core/mbelib.c")
    simd = pointer_size(build) == 8
    expect = "with" if simd else "without"
    check(("MBELIB_ENABLE_SIMD=1" in release) == simd, f"({name}) multi-config Release compiles {expect} SIMD")
    check(("MBELIB_ENABLE_SIMD=1" in debug) == simd, f"({name}) multi-config Debug compiles {expect} SIMD")
    if "IPO/LTO not supported" in log:
        print("note: the compiler does not support IPO; LTO checks skipped")
    else:
        check(has_lto(release), f"({name}) multi-config Release compiles with LTO")
    check(not has_lto(debug), f"({name}) multi-config Debug compiles without LTO")
    tuned = native_flags_of(log) is not None
    check_dispatch(args, build, log, not tuned, f"({name}) multi-config Release", "build-Release.ninja")
    check_dispatch(args, build, log, True, f"({name}) multi-config Debug", "build-Debug.ninja")


def scenario_debug_avx2_only(args, work):
    """(e) Only Debug targets AVX2: portable Release still dispatches."""
    if not args.user_isa.startswith("-march=x86-64"):
        print("note: not an x86-64 GNU-style compiler; (e) skipped")
        return
    build = work / "debug-avx2-only"
    log = configure(args, build, "Ninja Multi-Config", "-DMBELIB_BUILD_TESTS=OFF", "-DMBELIB_ENABLE_NATIVE=OFF",
                    "-DCMAKE_C_FLAGS_DEBUG=-O0 -g -mavx2 -mfma")
    if pointer_size(build) != 8:
        print("note: not a 64-bit target; (e) skipped")
        return
    check(dispatch_flags_of(log) is not None, "(e) the AVX2 dispatch stays on")
    check_dispatch(args, build, log, True, "(e) multi-config Release", "build-Release.ninja")


def scenario_release(args, work):
    build = work / "release"
    log = configure(args, build, "Ninja", "-DCMAKE_BUILD_TYPE=Release", "-DMBELIB_BUILD_TESTS=ON")
    flags = native_flags_of(log)
    if flags is None:
        reason = re.search(r"MBELIB_ENABLE_NATIVE: not tuning for the build host: (.*)", log)
        check(reason is not None, "(b) Release reports why it does not tune for the host")
        print(f"note: no native tuning on this host ({reason.group(1)})")
        return None
    flag = flags[0]
    check(all(f in library_tokens(args, build) for f in flags), f"(b) Release library compiles with {flags}")
    check_dispatch(args, build, log, False, "(b) tuned Release")
    golden = tokens_of(compile_line(commands(args, build, "test_golden_pcm"), "tests/test_golden_pcm.c"))
    check(flag in golden, f"(b) Release test_golden_pcm compiles with {flag}")
    api = tokens_of(compile_line(commands(args, build, "test_api"), "tests/test_api.c"))
    check(flag not in api, f"(b) Release test_api compiles without {flag}")
    debug = work / "debug"
    debug_log = configure(args, debug, "Ninja", "-DCMAKE_BUILD_TYPE=Debug", "-DMBELIB_BUILD_TESTS=OFF",
                          "-DMBELIB_ENABLE_SIMD=ON")
    check(flag not in library_tokens(args, debug), f"(b) Debug library compiles without {flag}")
    check_dispatch(args, debug, debug_log, True, "(b) Debug")
    return flag


def scenario_no_native(args, work, flag, name, *options):
    build = work / re.sub(r"\W+", "-", name)
    log = configure(args, build, "Ninja", "-DCMAKE_BUILD_TYPE=Release", "-DMBELIB_BUILD_TESTS=OFF", *options)
    check(native_flags_of(log) is None, f"{name}: configure reports no native flag")
    if flag:
        check(flag not in library_tokens(args, build), f"{name}: library compiles without {flag}")
    return build, log


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cmake", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--pffft-source", required=True)
    parser.add_argument("--c-compiler", required=True)
    parser.add_argument("--c-compiler-arg1", default="")
    parser.add_argument("--user-isa", default="", help="a baseline target ISA option the compiler accepts")
    parser.add_argument("--ninja", default="")
    args = parser.parse_args()
    check_rules(args)
    args.ninja = args.ninja or shutil.which("ninja") or shutil.which("ninja-build") or ""
    if not args.ninja:
        print("ninja not found; configure scenarios skipped")
        return 77
    if USER_ISA.search(args.c_compiler_arg1):
        # The scenarios look for host-tuning flags, which a compiler command
        # with its own ISA option (CC="gcc -march=x86-64") never gets.
        print("the forwarded compiler command already selects a target ISA; configure scenarios skipped")
        return 77
    with tempfile.TemporaryDirectory(prefix="mbelib-build-configuration-") as directory:
        work = Path(directory)
        scenario_multi_config(args, work, "a")
        # Multi-config generators ignore CMAKE_BUILD_TYPE.
        scenario_multi_config(args, work, "a-build-type", "-DCMAKE_BUILD_TYPE=Debug")
        # Ninja quotes paths with spaces; the source may live under one.
        spaced = work / "source with spaces"
        try:
            spaced.symlink_to(Path(args.source).resolve(), target_is_directory=True)
        except OSError:
            print("note: cannot create a symbolic link here; spaced-path scenario skipped")
        else:
            scenario_multi_config(args, work, "a-spaced-source", source=str(spaced))
        flag = scenario_release(args, work)
        if args.user_isa:
            scenario_no_native(args, work, flag, "(c) user target ISA", f"-DCMAKE_C_FLAGS={args.user_isa}")
            scenario_no_native(args, work, flag, "(c) ISA in the compiler command",
                               f"-DCMAKE_C_COMPILER_ARG1={args.c_compiler_arg1} {args.user_isa}")
        else:
            print("note: no known baseline target ISA option for this compiler; (c) skipped")
        build, log = scenario_no_native(args, work, flag, "(d) MBELIB_ENABLE_NATIVE=OFF",
                                        "-DMBELIB_ENABLE_NATIVE=OFF")
        check_dispatch(args, build, log, True, "(d) MBELIB_ENABLE_NATIVE=OFF")
        scenario_debug_avx2_only(args, work)
    return 0


if __name__ == "__main__":
    sys.exit(main())
