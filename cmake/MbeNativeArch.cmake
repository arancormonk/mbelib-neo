# SPDX-License-Identifier: GPL-2.0-or-later
# Host CPU tuning for release-like builds (MBELIB_ENABLE_NATIVE).
#
# mbelib_native_arch_flags(<flags_var> <reason_var>) sets <flags_var> to the
# compiler options that tune code generation for the build host, or to an
# empty list, and <reason_var> to why. It never tunes when the result could
# run on another machine or would override a target the user chose
# (mbelib_native_arch_skip_reason()):
#
# - cross-compiling, CMAKE_C_COMPILER_TARGET set, or targeting Android, iOS
#   (and the other Apple device platforms, also when only CMAKE_OSX_SYSROOT
#   selects their SDK) or Emscripten;
# - CMAKE_OSX_ARCHITECTURES naming more than one architecture, or one other
#   than the host's;
# - -march=, -mcpu=, /arch: or -arch already in the compiler command
#   (CMAKE_C_COMPILER_ARG1, as for CC="gcc -march=x86-64"), CMAKE_C_FLAGS or a
#   release-like CMAKE_C_FLAGS_<CONFIG> (distribution flags such as makepkg's
#   -march=x86-64 included), or for MSVC-style compilers in the CL and _CL_
#   environment variables;
# - a compiler that rejects the option, or MinGW/Cygwin GCC without the
#   assembler's -muse-unaligned-vector-move (GCC does not align the Windows
#   x64 stack for 32-byte spills, GCC PR 54412).
#
# GNU-style compilers (GCC, Clang, AppleClang) get -march=native on x86 and
# -mcpu=native on ARM. MSVC and clang-cl have no native option; on an x86 or
# x64 host they get the highest of /arch:AVX512, /arch:AVX2 and /arch:AVX that
# a CPUID and XGETBV probe run at configure time finds supported by the CPU and
# the OS. The probe result is cached for the build directory.
#
# mbelib_native_arch_skip_reason() and mbelib_native_arch_compose() only read
# variables, so tests run them in script mode (cmake -P).
include_guard(GLOBAL)

include(CheckCCompilerFlag)
include(CheckCSourceCompiles)

# Option forms that pick a target ISA, also when quoted, at the start or after
# a space or tab. A user choice wins over host tuning.
set(_MBELIB_NATIVE_USER_ISA_REGEX
    "(^|[ \t])[\"']?(-march=|-mcpu=|[-/][Aa][Rr][Cc][Hh]:|-arch[ \t])"
)

# mbelib_native_arch_skip_reason(<out_var>)
#
# Sets <out_var> to why host tuning must not apply, or to "" when it may.
function(mbelib_native_arch_skip_reason out_var)
    set(_reason "")
    if(CMAKE_CROSSCOMPILING)
        set(_reason "cross-compiling")
    elseif(
        EMSCRIPTEN
        OR CMAKE_SYSTEM_NAME
            MATCHES
            "^(Android|iOS|tvOS|watchOS|visionOS|Emscripten)$"
    )
        set(_reason "${CMAKE_SYSTEM_NAME} target")
    elseif(CMAKE_C_COMPILER_TARGET)
        set(_reason
            "CMAKE_C_COMPILER_TARGET selects ${CMAKE_C_COMPILER_TARGET}"
        )
    endif()

    # An iOS-family SDK can be selected through the sysroot alone, without
    # CMAKE_SYSTEM_NAME; resolved SDK paths are mixed case (iPhoneOS.sdk).
    string(TOLOWER "${CMAKE_OSX_SYSROOT}" _sysroot)
    if(
        _reason STREQUAL ""
        AND _sysroot MATCHES "(iphone|appletv|watch|xr)(os|simulator)"
    )
        set(_reason "CMAKE_OSX_SYSROOT selects ${CMAKE_OSX_SYSROOT}")
    endif()

    if(
        _reason STREQUAL ""
        AND (APPLE OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        AND CMAKE_OSX_ARCHITECTURES
    )
        list(LENGTH CMAKE_OSX_ARCHITECTURES _count)
        set(_host "${CMAKE_HOST_SYSTEM_PROCESSOR}")
        if(_host STREQUAL "aarch64")
            set(_host arm64)
        endif()
        if(_count GREATER 1)
            set(_reason
                "CMAKE_OSX_ARCHITECTURES names several architectures (${CMAKE_OSX_ARCHITECTURES})"
            )
        elseif(NOT CMAKE_OSX_ARCHITECTURES STREQUAL _host)
            set(_reason
                "CMAKE_OSX_ARCHITECTURES (${CMAKE_OSX_ARCHITECTURES}) is not the host architecture (${_host})"
            )
        endif()
    endif()

    if(_reason STREQUAL "")
        # CC="gcc -march=x86-64" keeps the option in CMAKE_C_COMPILER_ARG1.
        foreach(
            _var
            IN
            ITEMS
                CMAKE_C_COMPILER_ARG1
                CMAKE_C_FLAGS
                CMAKE_C_FLAGS_RELEASE
                CMAKE_C_FLAGS_RELWITHDEBINFO
                CMAKE_C_FLAGS_MINSIZEREL
        )
            if("${${_var}}" MATCHES "${_MBELIB_NATIVE_USER_ISA_REGEX}")
                set(_reason "${_var} already selects a target ISA")
                break()
            endif()
        endforeach()
    endif()

    # MSVC and clang-cl also take options from the CL and _CL_ environment
    # variables, which CMake does not copy into the flags above.
    if(
        _reason STREQUAL ""
        AND (MSVC OR CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    )
        foreach(_env IN ITEMS CL _CL_)
            if("$ENV{${_env}}" MATCHES "${_MBELIB_NATIVE_USER_ISA_REGEX}")
                set(_reason
                    "the ${_env} environment variable already selects a target ISA"
                )
                break()
            endif()
        endforeach()
    endif()
    set(${out_var} "${_reason}" PARENT_SCOPE)
endfunction()

# mbelib_native_arch_compose(<flags_var> <reason_var> <flag>)
#
# The GNU-style host-tuning options for an accepted <flag>. MinGW and Cygwin
# GCC also need the assembler to turn aligned vector moves into unaligned ones
# (MBELIB_C_HAS_UNALIGNED_VECTOR_MOVE); without it, tuning is declined.
function(mbelib_native_arch_compose flags_var reason_var flag)
    set(_flags "${flag}")
    set(_reason "host CPU")
    if((WIN32 OR CYGWIN) AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
        if(MBELIB_C_HAS_UNALIGNED_VECTOR_MOVE)
            list(APPEND _flags -Wa,-muse-unaligned-vector-move)
            set(_reason
                "host CPU; unaligned vector moves for the Windows stack"
            )
        else()
            set(_flags "")
            set(_reason
                "MinGW/Cygwin GCC needs the assembler option -muse-unaligned-vector-move (binutils 2.38 or later) for wide vectors on the Windows x64 stack"
            )
        endif()
    endif()
    set(${flags_var} "${_flags}" PARENT_SCOPE)
    set(${reason_var} "${_reason}" PARENT_SCOPE)
endfunction()

# Compile-only probe of the active compiler target.
function(_mbelib_native_target_is out_var source)
    set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
    check_c_source_compiles("${source}" ${out_var})
    set(${out_var} "${${out_var}}" PARENT_SCOPE)
endfunction()

# mbelib_native_probe_target(): the target architecture family, probed once
# per build directory: MBELIB_NATIVE_TARGET_X86 or MBELIB_NATIVE_TARGET_ARM
# (both false elsewhere).
macro(mbelib_native_probe_target)
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
endmacro()

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

function(mbelib_native_arch_flags flags_var reason_var)
    set(${flags_var} "" PARENT_SCOPE)
    mbelib_native_arch_skip_reason(_skip)
    if(_skip)
        set(${reason_var} "${_skip}" PARENT_SCOPE)
        return()
    endif()
    mbelib_native_probe_target()

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
    if((WIN32 OR CYGWIN) AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
        check_c_compiler_flag(
            "-Wa,-muse-unaligned-vector-move"
            MBELIB_C_HAS_UNALIGNED_VECTOR_MOVE
        )
    endif()
    mbelib_native_arch_compose(_flags _reason "${_flag}")
    set(${flags_var} "${_flags}" PARENT_SCOPE)
    set(${reason_var} "${_reason}" PARENT_SCOPE)
endfunction()
