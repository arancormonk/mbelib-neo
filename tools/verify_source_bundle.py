#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Extract, build, test, install and consume a source bundle (run with network disabled)."""

import argparse
from pathlib import Path
import subprocess
import tempfile

from release_sbom import package_files, validate_package


def run(*command, cwd):
    subprocess.run(command, cwd=cwd, check=True)


def check_dependency_preferences(source, work):
    dependency = next((source / "third_party").glob("pffft-*.tar.gz"))
    original = dependency.read_bytes()
    extracted = source / "build/dev-debug/_deps/mbelib_pffft-src"
    try:
        dependency.write_bytes(b"intentionally corrupted dependency test")
        command = ["cmake", "-S", str(source), "-B", str(work / "reject-corrupt")]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode == 0 or "Bundled PFFFT SHA256" not in result.stderr:
            raise ValueError("CMake did not reject the corrupted bundled dependency")
        run("cmake", "-S", str(source), "-B", str(work / "explicit-override"),
            f"-DFETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT={extracted}", cwd=work)
    finally:
        dependency.write_bytes(original)


def verify(archive, work):
    files = package_files(archive)
    _, _, identity, _ = validate_package(files, "source")
    source = work / "source"
    source.mkdir()
    modes = {f["path"]: f["mode"] for f in identity["files"]}
    for name, data in files.items():
        target = source / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        target.chmod(modes.get(name, 0o644))
    install = work / "install"
    run("cmake", "--preset", "dev-debug", "-DFETCHCONTENT_FULLY_DISCONNECTED=ON", cwd=source)
    run("cmake", "--build", "--preset", "dev-debug", "--parallel", "4", cwd=source)
    run("ctest", "--preset", "dev-debug", "--output-on-failure", cwd=source)
    check_dependency_preferences(source, work)
    run("cmake", "--install", "build/dev-debug", "--prefix", str(install), cwd=source)
    consumer = work / "consumer"
    consumer.mkdir()
    (consumer / "main.c").write_text(
        '// SPDX-License-Identifier: GPL-2.0-or-later\n#include <mbelib-neo/mbelib.h>\n'
        'int main(void) { mbe_parms a, b, c; mbe_initMbeParms(&a, &b, &c); return a.L < 9; }\n')
    (consumer / "CMakeLists.txt").write_text(
        '# SPDX-License-Identifier: GPL-2.0-or-later\ncmake_minimum_required(VERSION 3.20)\n'
        'project(consumer C)\nfind_package(mbe-neo CONFIG REQUIRED)\n'
        'add_executable(consumer main.c)\ntarget_link_libraries(consumer PRIVATE mbe_neo::mbe_shared)\n'
        'add_executable(consumer_static main.c)\ntarget_link_libraries(consumer_static PRIVATE mbe_neo::mbe_static)\n')
    run("cmake", "-S", ".", "-B", "build", f"-DCMAKE_PREFIX_PATH={install}", cwd=consumer)
    run("cmake", "--build", "build", cwd=consumer)
    run(str(consumer / "build/consumer"), cwd=consumer)
    run(str(consumer / "build/consumer_static"), cwd=consumer)
    print(f"Verified offline source build and installed consumers: {identity['tag']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="mbelib-source-verify-") as directory:
        verify(args.archive.resolve(), Path(directory))


if __name__ == "__main__":
    main()
