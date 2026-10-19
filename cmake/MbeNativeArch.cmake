# SPDX-License-Identifier: GPL-2.0-or-later
# Host CPU tuning for release-like builds (MBELIB_ENABLE_NATIVE).
#
# mbelib_native_arch_flags(<flags_var> <reason_var>) sets <flags_var> to the
# compiler options that tune code generation for the build host, or to an
# empty list, and <reason_var> to why. It never tunes when the result could
# run on another machine or would override a target the user chose:
#
# - cross-compiling, or targeting Android, iOS (and the other Apple device
#   platforms) or Emscripten;
# - CMAKE_OSX_ARCHITECTURES naming more than one architecture, or one other
#   than the host's;
# - -march=, -mcpu= or /arch: already in CMAKE_C_FLAGS or a release-like
#   CMAKE_C_FLAGS_<CONFIG> (distribution flags such as makepkg's
#   -march=x86-64 included);
# - a compiler that rejects the option.
#
# GNU-style compilers (GCC, Clang, AppleClang) get -march=native on x86 and
# -mcpu=native on ARM. MSVC and clang-cl have no native option; on an x86 or
# x64 host they get the highest of /arch:AVX512, /arch:AVX2 and /arch:AVX that
# a CPUID and XGETBV probe run at configure time finds supported by the CPU and
# the OS. The probe result is cached for the build directory.
include_guard(GLOBAL)

include(CheckCCompilerFlag)
include(CheckCSourceCompiles)

# Compile-only probe of the active compiler target.
function(_mbelib_native_target_is out_var source)
    set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
    check_c_source_compiles("${source}" ${out_var})
    set(${out_var} "${${out_var}}" PARENT_SCOPE)
endfunction()

function(_mbelib_native_user_isa out_var)
    set(_found "")
    foreach(
        _flag_var
        IN
        ITEMS
            CMAKE_C_FLAGS
            CMAKE_C_FLAGS_RELEASE
            CMAKE_C_FLAGS_RELWITHDEBINFO
            CMAKE_C_FLAGS_MINSIZEREL
    )
        if("${${_flag_var}}" MATCHES "(-march=|-mcpu=|[/-]arch:)")
            set(_found "${_flag_var} already selects a target ISA")
            break()
        endif()
    endforeach()
    set(${out_var} "${_found}" PARENT_SCOPE)
endfunction()

function(_mbelib_native_osx_mismatch out_var)
    set(_mismatch "")
    if(APPLE AND CMAKE_OSX_ARCHITECTURES)
        list(LENGTH CMAKE_OSX_ARCHITECTURES _count)
        set(_host "${CMAKE_HOST_SYSTEM_PROCESSOR}")
        if(_host STREQUAL "aarch64")
            set(_host arm64)
        endif()
        if(_count GREATER 1)
            set(_mismatch
                "CMAKE_OSX_ARCHITECTURES names several architectures (${CMAKE_OSX_ARCHITECTURES})"
            )
        elseif(NOT CMAKE_OSX_ARCHITECTURES STREQUAL _host)
            set(_mismatch
                "CMAKE_OSX_ARCHITECTURES (${CMAKE_OSX_ARCHITECTURES}) is not the host architecture (${_host})"
            )
        endif()
    endif()
    set(${out_var} "${_mismatch}" PARENT_SCOPE)
endfunction()

# Highest /arch: level the build host supports: AVX512, AVX2, AVX or none.
# AVX2 needs the whole x86-64-v3 set because clang-cl's /arch:AVX2 also emits
# FMA, F16C, MOVBE, BMI1/2 and LZCNT; AVX512 needs F, CD, BW, DQ and VL as
# well, and every level needs the OS to save the registers it uses (XCR0).
function(_mbelib_native_msvc_probe out_var)
    if(DEFINED MBELIB_NATIVE_MSVC_ARCH_PROBED)
        set(${out_var} "${MBELIB_NATIVE_MSVC_ARCH}" PARENT_SCOPE)
        return()
    endif()
    set(_dir "${CMAKE_BINARY_DIR}/CMakeFiles/MbeNativeArch")
    file(MAKE_DIRECTORY "${_dir}")
    file(
        WRITE "${_dir}/msvc_arch_probe.c"
        [=[
#include <intrin.h>

#if defined(__clang__) || defined(__GNUC__)
__attribute__((target("xsave")))
#endif
static unsigned long long
probe_xgetbv(void) {
    return _xgetbv(0);
}

static int
bit(int reg, int b) {
    return (int)(((unsigned)reg >> b) & 1u);
}

int
main(void) {
    int r[4];
    __cpuid(r, 0);
    const int max_leaf = r[0];
    if (max_leaf < 1) {
        return 0;
    }
    __cpuid(r, 1);
    const int ecx = r[2];
    /* OSXSAVE, AVX, SSE3, SSSE3, SSE4.1, SSE4.2, POPCNT */
    if (!bit(ecx, 27) || !bit(ecx, 28) || !bit(ecx, 0) || !bit(ecx, 9) || !bit(ecx, 19) || !bit(ecx, 20)
        || !bit(ecx, 23)) {
        return 0;
    }
    const unsigned long long xcr0 = probe_xgetbv();
    if ((xcr0 & 0x6u) != 0x6u) {
        return 0;
    }
    if (max_leaf < 7) {
        return 1;
    }
    __cpuidex(r, 7, 0);
    const int ebx = r[1];
    __cpuid(r, (int)0x80000000u);
    int lzcnt = 0;
    if ((unsigned)r[0] >= 0x80000001u) {
        __cpuid(r, (int)0x80000001u);
        lzcnt = bit(r[2], 5);
    }
    /* FMA, MOVBE, F16C, BMI1, AVX2, BMI2, LZCNT */
    if (!bit(ecx, 12) || !bit(ecx, 22) || !bit(ecx, 29) || !bit(ebx, 3) || !bit(ebx, 5) || !bit(ebx, 8) || !lzcnt) {
        return 1;
    }
    /* AVX-512 F, DQ, CD, BW, VL with opmask, ZMM0-15 upper and ZMM16-31 state */
    if (bit(ebx, 16) && bit(ebx, 17) && bit(ebx, 28) && bit(ebx, 30) && bit(ebx, 31) && (xcr0 & 0xE0u) == 0xE0u) {
        return 3;
    }
    return 2;
}
]=]
    )
    try_run(
        _run_result
        _compile_result
        "${_dir}/build"
        "${_dir}/msvc_arch_probe.c"
        COMPILE_OUTPUT_VARIABLE _compile_output
    )
    set(_arch "")
    if(_compile_result AND _run_result MATCHES "^[0-3]$")
        if(_run_result EQUAL 3)
            set(_arch "/arch:AVX512")
        elseif(_run_result EQUAL 2)
            set(_arch "/arch:AVX2")
        elseif(_run_result EQUAL 1)
            set(_arch "/arch:AVX")
        endif()
    endif()
    set(MBELIB_NATIVE_MSVC_ARCH
        "${_arch}"
        CACHE INTERNAL
        "Highest /arch: level found by the configure-time CPU probe"
    )
    set(MBELIB_NATIVE_MSVC_ARCH_PROBED
        ON
        CACHE INTERNAL
        "The MSVC native-arch CPU probe has run"
    )
    set(${out_var} "${_arch}" PARENT_SCOPE)
endfunction()

# The target architecture family, probed once: MBELIB_NATIVE_TARGET_X86 or
# MBELIB_NATIVE_TARGET_ARM (both false elsewhere).
_mbelib_native_target_is(MBELIB_NATIVE_TARGET_X86
  [=[
  #if !(defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)) || defined(_M_ARM64EC)
  #error not an x86 target
  #endif
  int mbelib_native_target_x86;
]=]
)
_mbelib_native_target_is(MBELIB_NATIVE_TARGET_ARM
  [=[
  #if !(defined(__aarch64__) || defined(__arm__) || defined(_M_ARM64) || defined(_M_ARM) || defined(_M_ARM64EC))
  #error not an ARM target
  #endif
  int mbelib_native_target_arm;
]=]
)

function(mbelib_native_arch_flags flags_var reason_var)
    set(${flags_var} "" PARENT_SCOPE)
    if(CMAKE_CROSSCOMPILING)
        set(${reason_var} "cross-compiling" PARENT_SCOPE)
        return()
    endif()
    if(
        EMSCRIPTEN
        OR CMAKE_SYSTEM_NAME
            MATCHES
            "^(Android|iOS|tvOS|watchOS|visionOS|Emscripten)$"
    )
        set(${reason_var} "${CMAKE_SYSTEM_NAME} target" PARENT_SCOPE)
        return()
    endif()
    _mbelib_native_osx_mismatch(_osx)
    if(_osx)
        set(${reason_var} "${_osx}" PARENT_SCOPE)
        return()
    endif()
    _mbelib_native_user_isa(_user)
    if(_user)
        set(${reason_var} "${_user}" PARENT_SCOPE)
        return()
    endif()

    if(MSVC)
        if(NOT MBELIB_NATIVE_TARGET_X86)
            set(${reason_var}
                "MSVC-style compiler: no host probe for this target architecture"
                PARENT_SCOPE
            )
            return()
        endif()
        _mbelib_native_msvc_probe(_arch)
        if(NOT _arch)
            set(${reason_var}
                "the host CPU probe found no /arch: level above the default"
                PARENT_SCOPE
            )
            return()
        endif()
        set(${flags_var} "${_arch}" PARENT_SCOPE)
        set(${reason_var} "host CPU probe" PARENT_SCOPE)
        return()
    endif()

    if(NOT CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
        set(${reason_var}
            "no native tuning option for ${CMAKE_C_COMPILER_ID}"
            PARENT_SCOPE
        )
        return()
    endif()
    if(MBELIB_NATIVE_TARGET_X86)
        set(_flag -march=native)
        set(_check MBELIB_C_HAS_MARCH_NATIVE)
    elseif(MBELIB_NATIVE_TARGET_ARM)
        set(_flag -mcpu=native)
        set(_check MBELIB_C_HAS_MCPU_NATIVE)
    else()
        set(${reason_var}
            "no native tuning option for this target architecture"
            PARENT_SCOPE
        )
        return()
    endif()
    check_c_compiler_flag("${_flag}" ${_check})
    if(NOT ${_check})
        set(${reason_var}
            "${CMAKE_C_COMPILER_ID} rejects ${_flag}"
            PARENT_SCOPE
        )
        return()
    endif()
    set(${flags_var} "${_flag}" PARENT_SCOPE)
    set(${reason_var} "host CPU" PARENT_SCOPE)
endfunction()
