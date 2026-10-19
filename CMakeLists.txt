# SPDX-License-Identifier: GPL-2.0-or-later
cmake_minimum_required(VERSION 3.20)

project(mbelib-neo VERSION 2.6.0 LANGUAGES C)
if(POLICY CMP0083)
    cmake_policy(SET CMP0083 NEW)
endif()

# Language and build defaults
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCCompilerFlag)
include(CheckCSourceCompiles)
include(CheckLinkerFlag)

function(_mbelib_check_c_source_compiles_static source out_var)
    set(_mbelib_had_try_compile_target_type OFF)
    if(DEFINED CMAKE_TRY_COMPILE_TARGET_TYPE)
        set(_mbelib_had_try_compile_target_type ON)
        set(_mbelib_saved_try_compile_target_type
            "${CMAKE_TRY_COMPILE_TARGET_TYPE}"
        )
    endif()

    set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
    check_c_source_compiles("${source}" ${out_var})

    if(_mbelib_had_try_compile_target_type)
        set(CMAKE_TRY_COMPILE_TARGET_TYPE
            "${_mbelib_saved_try_compile_target_type}"
        )
    else()
        unset(CMAKE_TRY_COMPILE_TARGET_TYPE)
    endif()

    set(${out_var} "${${out_var}}" PARENT_SCOPE)
endfunction()

# Options
option(NOTONES "Disable AMBE+2 Tones" OFF)
option(MBELIB_BUILD_TESTS "Build mbelib tests" ON)
option(MBELIB_BUILD_EXAMPLES "Build examples" ON)
option(MBELIB_BUILD_DOCS "Generate API docs with Doxygen" OFF)
option(MBELIB_ENABLE_ASAN "Enable AddressSanitizer in Debug builds" OFF)
option(
    MBELIB_ENABLE_UBSAN
    "Enable UndefinedBehaviorSanitizer in Debug builds"
    OFF
)
option(MBE_ENABLE_DEBUG_LOGS "Enable debug logging in codec sources" OFF)
option(MBELIB_ENABLE_WARNINGS "Enable common warning flags" ON)
option(MBELIB_WARNINGS_AS_ERRORS "Treat warnings as errors" ON)
option(
    MBELIB_ENABLE_FAST_MATH
    "Enable fast-math optimizations (may relax IEEE semantics)"
    OFF
)
option(
    MBELIB_ENABLE_LTO
    "Enable Link Time Optimization (IPO/LTO) for Release builds"
    OFF
)
option(
    MBELIB_ENABLE_HARDENING
    "Enable Release-like compiler and linker hardening when supported"
    ON
)
# Release builds of 64-bit targets default to SIMD. 32-bit targets stay
# opt-in: x86 would require SSE2 and ARM needs an explicit -mfpu=neon. A
# multi-config generator (Ninja Multi-Config, Visual Studio, Xcode) picks the
# configuration at build time and ignores CMAKE_BUILD_TYPE, so it gets the
# Release default; its Debug configuration is the hash-covered Debug + SIMD
# build.
get_property(_mbelib_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
set(_mbelib_simd_default OFF)
if(
    CMAKE_SIZEOF_VOID_P EQUAL 8
    AND (
        _mbelib_multi_config
        OR CMAKE_BUILD_TYPE MATCHES "Release|RelWithDebInfo|MinSizeRel"
    )
)
    set(_mbelib_simd_default ON)
endif()
option(
    MBELIB_ENABLE_SIMD
    "Enable SIMD optimized routines when available (SSE2/NEON)"
    ${_mbelib_simd_default}
)
unset(_mbelib_simd_default)
# Release-like configurations of the libraries are tuned for the build host's
# CPU unless that could misfire (cmake/MbeNativeArch.cmake lists the cases).
# Packagers and CI that ship binaries pass -DMBELIB_ENABLE_NATIVE=OFF.
option(
    MBELIB_ENABLE_NATIVE
    "Tune release-like library builds for the build host's CPU when safe"
    ON
)
option(MBELIB_BUILD_BENCHMARKS "Build micro-benchmarks (not run in CI)" OFF)
option(MBELIB_BUILD_TOOLS "Build developer quality-analysis tools" OFF)
option(
    MBELIB_TEST_REQUIRE_AVX2
    "Fail, rather than skip, the AVX2 kernel tests when the build or CPU cannot run AVX2 kernels"
    OFF
)
# Warnings
if(MBELIB_ENABLE_WARNINGS)
    if(MSVC)
        foreach(
            _mbelib_msvc_flag_var
            IN
            ITEMS
                CMAKE_C_FLAGS
                CMAKE_C_FLAGS_DEBUG
                CMAKE_C_FLAGS_RELEASE
                CMAKE_C_FLAGS_RELWITHDEBINFO
                CMAKE_C_FLAGS_MINSIZEREL
        )
            string(
                REGEX REPLACE "(^| )/W[0-4]( |$)"
                " "
                ${_mbelib_msvc_flag_var}
                "${${_mbelib_msvc_flag_var}}"
            )
        endforeach()
        # MSVC reports large legacy DSP tables and vendored PFFFT scalar math as
        # noisy conversion diagnostics. Keep /WX for the rest of /W4.
        add_compile_options(/W4 /wd4244 /wd4305)
        if(MBELIB_WARNINGS_AS_ERRORS)
            add_compile_options(/WX)
        endif()
    elseif(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
        add_compile_options(
            -Wall
            -Wextra
            -Wpedantic
            -Wformat=2
            -Wvla
            -Wstrict-prototypes
            -Wmissing-prototypes
            -Wold-style-definition
        )
        if(MBELIB_WARNINGS_AS_ERRORS)
            add_compile_options(-Werror)
        endif()
    endif()
endif()

# Debug builds round every multiply and add separately, in the library and in
# the tests and tools that build its golden fixtures, so the exact golden PCM
# hashes hold for every compiler: Clang fuses a * b + c into an FMA by default,
# even at -O0, wherever the target has FMA (x86-64-v3, AArch64), while GCC
# contracts only when optimizing.
if(
    CMAKE_C_COMPILER_ID MATCHES "Clang|GNU"
    AND NOT CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC"
)
    add_compile_options($<$<CONFIG:Debug>:-ffp-contract=off>)
endif()

set(_mbelib_release_like_config
    "$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>,$<CONFIG:MinSizeRel>>"
)

function(_mbelib_flags_define_fortify out_var)
    set(_mbelib_fortify_defined OFF)
    foreach(
        _mbelib_flag_var
        IN
        ITEMS
            CMAKE_C_FLAGS
            CMAKE_C_FLAGS_RELEASE
            CMAKE_C_FLAGS_RELWITHDEBINFO
            CMAKE_C_FLAGS_MINSIZEREL
    )
        if(
            DEFINED ${_mbelib_flag_var}
            AND "${${_mbelib_flag_var}}" MATCHES "_FORTIFY_SOURCE"
        )
            set(_mbelib_fortify_defined ON)
            break()
        endif()
    endforeach()
    set(${out_var} "${_mbelib_fortify_defined}" PARENT_SCOPE)
endfunction()

if(MBELIB_ENABLE_HARDENING AND (CMAKE_C_COMPILER_ID MATCHES "Clang|GNU"))
    check_c_compiler_flag(
        "-fstack-protector-strong"
        MBELIB_C_HAS_STACK_PROTECTOR_STRONG
    )
    if(NOT APPLE)
        check_c_compiler_flag(
            "-fstack-clash-protection"
            MBELIB_C_HAS_STACK_CLASH_PROTECTION
        )
    endif()
    if(UNIX AND NOT APPLE)
        check_linker_flag(C "-Wl,-z,relro" MBELIB_LD_HAS_RELRO)
        check_linker_flag(C "-Wl,-z,now" MBELIB_LD_HAS_NOW)
        include(CheckPIESupported)
        check_pie_supported(LANGUAGES C)
    endif()
endif()

function(_mbelib_apply_hardening target)
    if(NOT MBELIB_ENABLE_HARDENING)
        return()
    endif()

    get_target_property(_mbelib_target_type ${target} TYPE)

    if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
        _mbelib_flags_define_fortify(_mbelib_c_flags_define_fortify)
        if(NOT _mbelib_c_flags_define_fortify)
            target_compile_definitions(
                ${target}
                PRIVATE "$<${_mbelib_release_like_config}:_FORTIFY_SOURCE=2>"
            )
        endif()

        if(MBELIB_C_HAS_STACK_PROTECTOR_STRONG)
            target_compile_options(
                ${target}
                PRIVATE
                    "$<${_mbelib_release_like_config}:-fstack-protector-strong>"
            )
        endif()
        if(NOT APPLE AND MBELIB_C_HAS_STACK_CLASH_PROTECTION)
            target_compile_options(
                ${target}
                PRIVATE
                    "$<${_mbelib_release_like_config}:-fstack-clash-protection>"
            )
        endif()

        if(
            UNIX
            AND NOT APPLE
            AND NOT _mbelib_target_type STREQUAL "STATIC_LIBRARY"
            AND NOT _mbelib_target_type STREQUAL "OBJECT_LIBRARY"
        )
            if(MBELIB_LD_HAS_RELRO)
                target_link_options(
                    ${target}
                    PRIVATE "$<${_mbelib_release_like_config}:-Wl,-z,relro>"
                )
            endif()
            if(MBELIB_LD_HAS_NOW)
                target_link_options(
                    ${target}
                    PRIVATE "$<${_mbelib_release_like_config}:-Wl,-z,now>"
                )
            endif()
        endif()
    elseif(
        MSVC
        AND NOT _mbelib_target_type STREQUAL "STATIC_LIBRARY"
        AND NOT _mbelib_target_type STREQUAL "OBJECT_LIBRARY"
    )
        target_link_options(
            ${target}
            PRIVATE
                "$<${_mbelib_release_like_config}:/DYNAMICBASE>"
                "$<${_mbelib_release_like_config}:/NXCOMPAT>"
        )
        if(CMAKE_SIZEOF_VOID_P EQUAL 8)
            target_link_options(
                ${target}
                PRIVATE "$<${_mbelib_release_like_config}:/HIGHENTROPYVA>"
            )
        endif()
    endif()
endfunction()

# Generated headers (version)
include(CMakePackageConfigHelpers)
set(MBELIB_GENERATED_INCLUDE_DIR "${PROJECT_BINARY_DIR}/include")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/mbelib-neo/version.h.in"
    "${MBELIB_GENERATED_INCLUDE_DIR}/mbelib-neo/version.h"
    @ONLY
)

# Pinned, checksum-verified PFFFT sources.
include(cmake/Pffft.cmake)

# Explicit source lists (avoids FILE(GLOB ...))
set(MBELIB_FIRST_PARTY_SOURCES
    src/ambe/ambe3600x2400.c
    src/ambe/ambe3600x2400_enc.c
    src/ambe/dstar_frame.c
    src/ambe/ambe3600x2450.c
    src/ambe/ambe3600x2450_enc.c
    src/ambe/ambe_common.c
    src/ambe/ambe_encoder.c
    src/ecc/ecc_soft.c
    src/ecc/ecc_soft_search.c
    src/ecc/ecc_soft_search_avx2.c
    src/ecc/ecc.c
    src/ecc/ecc_const.c
    src/imbe/imbe7100x4400.c
    src/imbe/imbe7200x4400.c
    src/imbe/imbe7200x4400_enc.c
    src/core/mbelib.c
    src/core/mbe_cpu.c
    src/core/mbe_voiced.c
    src/core/mbe_voiced_avx2.c
    src/core/mbe_pcm_convert.c
    src/core/mbe_adaptive.c
    src/core/mbe_encoder.c
    src/core/mbe_speech_analysis.c
    src/core/mbe_tone_detect.c
    src/core/mbe_tone_fit.c
    src/core/mbe_tone_fit_avx2.c
    src/core/mbe_unvoiced_fft.c
    src/tetra/tetra_acelp_codebooks.c
    src/tetra/tetra_acelp_decoder.c
    src/tetra/tetra_acelp_dsp.c
    src/tetra/tetra_acelp_avx2.c
)
set(MBELIB_SOURCES ${MBELIB_FIRST_PARTY_SOURCES} ${PFFFT_SOURCES})

add_library(mbe-shared SHARED ${MBELIB_SOURCES})
add_library(mbe-static STATIC ${MBELIB_SOURCES})
add_library(
    mbe-codeql-first-party
    OBJECT
    EXCLUDE_FROM_ALL
    ${MBELIB_FIRST_PARTY_SOURCES}
)

add_library(mbe_neo::mbe_shared ALIAS mbe-shared)
add_library(mbe_neo::mbe_static ALIAS mbe-static)

target_compile_definitions(mbe-shared PRIVATE MBE_BUILDING=1)
foreach(_target mbe-shared mbe-static mbe-codeql-first-party)
    target_compile_definitions(
        ${_target}
        PRIVATE PFFFT_STATIC_DEFINE $<$<CONFIG:Debug>:PFFFT_DISABLE_FMA=1>
    )
endforeach()

target_include_directories(
    mbe-shared
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${MBELIB_GENERATED_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/src/internal
        ${PFFFT_INCLUDE_DIR}
)
target_include_directories(
    mbe-static
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${MBELIB_GENERATED_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/src/internal
        ${PFFFT_INCLUDE_DIR}
)
target_include_directories(
    mbe-codeql-first-party
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${MBELIB_GENERATED_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/src/internal
        ${PFFFT_INCLUDE_DIR}
)

# Mark static consumption for Windows to avoid dllimport decorations in headers
target_compile_definitions(mbe-static PUBLIC MBE_STATIC=1)
_mbelib_apply_hardening(mbe-shared)
_mbelib_apply_hardening(mbe-static)

# Optional verbose debug logging in library sources
if(MBE_ENABLE_DEBUG_LOGS)
    target_compile_definitions(
        mbe-shared
        PRIVATE MBE_DEBUG IMBE_DEBUG AMBE_DEBUG
    )
    target_compile_definitions(
        mbe-static
        PRIVATE MBE_DEBUG IMBE_DEBUG AMBE_DEBUG
    )
endif()

# Target-scoped compile definitions
if(NOTONES)
    target_compile_definitions(mbe-shared PRIVATE DISABLE_AMBE_TONES)
    target_compile_definitions(mbe-static PRIVATE DISABLE_AMBE_TONES)
    target_compile_definitions(
        mbe-codeql-first-party
        PRIVATE DISABLE_AMBE_TONES
    )
endif()
# Documentation (Doxygen)
if(MBELIB_BUILD_DOCS)
    find_package(Doxygen QUIET)
    if(DOXYGEN_FOUND)
        add_custom_target(
            docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM
        )
    else()
        message(STATUS "Doxygen not found; docs target will not be available.")
    endif()
endif()

if(MSVC)
    target_compile_definitions(mbe-shared PRIVATE _USE_MATH_DEFINES)
    target_compile_definitions(mbe-static PRIVATE _USE_MATH_DEFINES)
endif()

if(NOT WIN32)
    target_link_libraries(mbe-static m)
    target_link_libraries(mbe-shared m)
endif()

# Detect 32-bit x86 targets so the SIMD option can ensure SSE2 codegen when
# needed instead of relying on the host processor string.
set(MBELIB_TARGET_X86_32 OFF)
set(_mbelib_x86_32_probe_source
    [=[
    #if !defined(__i386__) && !defined(_M_IX86)
    #error Not a 32-bit x86 target
    #endif
    int mbelib_target_x86_32;
  ]=]
)
set(_mbelib_x86_32_sse2_probe_source
    [=[
    #if !defined(__i386__) && !defined(_M_IX86)
    #error Not a 32-bit x86 target
    #endif
    #if !defined(__SSE2__) && !(defined(_M_IX86_FP) && (_M_IX86_FP >= 2))
    #error 32-bit x86 target is not compiling for SSE2
    #endif
    int mbelib_target_x86_32_sse2;
  ]=]
)
if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # Query the active compiler target instead of the host processor string so
    # multilib builds (for example -m32 on an x86_64 host) are classified
    # correctly during configure.
    _mbelib_check_c_source_compiles_static("${_mbelib_x86_32_probe_source}" _mbelib_compiler_targets_x86_32)
    if(_mbelib_compiler_targets_x86_32)
        set(MBELIB_TARGET_X86_32 ON)
    endif()
    unset(_mbelib_compiler_targets_x86_32 CACHE)
    unset(_mbelib_compiler_targets_x86_32)
endif()

# Optional SIMD toggles
if(MBELIB_ENABLE_SIMD)
    target_compile_definitions(mbe-shared PRIVATE MBELIB_ENABLE_SIMD=1)
    target_compile_definitions(mbe-static PRIVATE MBELIB_ENABLE_SIMD=1)
    # Mirror the fork's backend selection, including
    # PFFFT-only backends (AltiVec and SSE1).
    set(_mbelib_pffft_simd_probe_source
        [=[
        #if (defined(__ppc__) || defined(__ppc64__) || defined(__powerpc__) || defined(__powerpc64__)) && (defined(__VEC__) || defined(__ALTIVEC__))
        /* AltiVec */
        #elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64) || defined(_M_ARM64EC) || defined(__wasm_simd128__)
        /* NEON or native WebAssembly SIMD */
        #elif defined(__x86_64__) || defined(__SSE__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 1)
        /* SSE */
        #else
        #error No supported PFFFT SIMD backend
        #endif
        int mbelib_target_pffft_simd;
      ]=]
    )
    _mbelib_check_c_source_compiles_static("${_mbelib_pffft_simd_probe_source}" _mbelib_pffft_simd_available)
    # The fork gates NEON/WASM behind opt-in macros. Its headers additionally
    # check compiler target features; no host-based ISA flags are introduced.
    foreach(_target mbe-shared mbe-static)
        target_compile_definitions(
            ${_target}
            PRIVATE PFFFT_ENABLE_NEON=1 PFFFT_ENABLE_WASM=1
        )
    endforeach()
    if(MBELIB_TARGET_X86_32)
        set(_mbelib_x86_sse2_flag "")
        _mbelib_check_c_source_compiles_static("${_mbelib_x86_32_sse2_probe_source}" _mbelib_x86_targets_sse2_by_default)
        if(_mbelib_x86_targets_sse2_by_default)
            message(
                STATUS
                "MBELIB_ENABLE_SIMD on 32-bit x86: active compiler target already includes SSE2"
            )
        else()
            set(_mbelib_x86_sse2_candidate_flags -msse2 /arch:SSE2)
            if(MSVC)
                set(_mbelib_x86_sse2_candidate_flags /arch:SSE2 -msse2)
            endif()

            foreach(
                _mbelib_candidate_flag
                IN
                LISTS _mbelib_x86_sse2_candidate_flags
            )
                set(_mbelib_saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
                if(CMAKE_REQUIRED_FLAGS)
                    set(CMAKE_REQUIRED_FLAGS
                        "${CMAKE_REQUIRED_FLAGS} ${_mbelib_candidate_flag}"
                    )
                else()
                    set(CMAKE_REQUIRED_FLAGS "${_mbelib_candidate_flag}")
                endif()

                _mbelib_check_c_source_compiles_static(
                  "${_mbelib_x86_32_sse2_probe_source}"
                  _mbelib_x86_sse2_candidate_supported
                )

                set(CMAKE_REQUIRED_FLAGS "${_mbelib_saved_required_flags}")
                set(_mbelib_candidate_flag_supported
                    "${_mbelib_x86_sse2_candidate_supported}"
                )
                unset(_mbelib_x86_sse2_candidate_supported CACHE)
                unset(_mbelib_x86_sse2_candidate_supported)

                if(_mbelib_candidate_flag_supported)
                    set(_mbelib_x86_sse2_flag "${_mbelib_candidate_flag}")
                    break()
                endif()
            endforeach()

            if(NOT _mbelib_x86_sse2_flag)
                message(
                    FATAL_ERROR
                    "MBELIB_ENABLE_SIMD on 32-bit x86 requires the active compiler target to enable SSE2 "
                    "(for example via -msse2 or /arch:SSE2)"
                )
            endif()

            message(
                STATUS
                "MBELIB_ENABLE_SIMD on 32-bit x86: compiling library targets for SSE2 via ${_mbelib_x86_sse2_flag}"
            )
            target_compile_options(mbe-shared PRIVATE ${_mbelib_x86_sse2_flag})
            target_compile_options(mbe-static PRIVATE ${_mbelib_x86_sse2_flag})
        endif()
        unset(_mbelib_x86_sse2_candidate_flags)
        unset(_mbelib_x86_targets_sse2_by_default CACHE)
        unset(_mbelib_x86_targets_sse2_by_default)
        unset(_mbelib_x86_sse2_flag)
        unset(_mbelib_candidate_flag)
        unset(_mbelib_candidate_flag_supported)
        unset(_mbelib_saved_required_flags)
    endif()
    # The x86-32 handling above guarantees SSE2, even if the default target
    # probed without the library's additional SSE2 flag did not provide SSE.
    if(NOT _mbelib_pffft_simd_available AND NOT MBELIB_TARGET_X86_32)
        target_compile_definitions(mbe-shared PRIVATE PFFFT_SIMD_DISABLE=1)
        target_compile_definitions(mbe-static PRIVATE PFFFT_SIMD_DISABLE=1)
        set(_mbelib_pffft_fallback_hint "")
        if(
            CMAKE_SIZEOF_VOID_P EQUAL 4
            AND CMAKE_SYSTEM_PROCESSOR MATCHES "^[Aa][Rr][Mm]"
        )
            set(_mbelib_pffft_fallback_hint
                " On 32-bit ARM, pass -mfpu=neon in CMAKE_C_FLAGS to enable NEON."
            )
        endif()
        message(
            STATUS
            "MBELIB_ENABLE_SIMD: FFT falls back to scalar for this target.${_mbelib_pffft_fallback_hint}"
        )
        unset(_mbelib_pffft_fallback_hint)
    endif()
    unset(_mbelib_pffft_simd_available CACHE)
    unset(_mbelib_pffft_simd_available)
    unset(_mbelib_pffft_simd_probe_source)
else()
    # Keep bundled PFFFT scalar when SIMD is disabled for the project.
    target_compile_definitions(mbe-shared PRIVATE PFFFT_SIMD_DISABLE=1)
    target_compile_definitions(mbe-static PRIVATE PFFFT_SIMD_DISABLE=1)
endif()

# Host CPU tuning (MBELIB_ENABLE_NATIVE). The flags go PRIVATE to the
# libraries, and to the tests that must pick the same compile-time kernels,
# for release-like configurations only: Debug keeps the baseline target its
# golden hashes, sanitizers and coverage expect, and nothing reaches the
# exported targets or pkg-config metadata.
include(cmake/MbeNativeArch.cmake)
mbelib_native_probe_target()
set(_mbelib_native_flags "")
if(MBELIB_ENABLE_NATIVE)
    mbelib_native_arch_flags(_mbelib_native_flags _mbelib_native_reason)
    list(JOIN _mbelib_native_flags " " _mbelib_native_flags_text)
    if(_mbelib_native_flags)
        message(
            STATUS
            "MBELIB_ENABLE_NATIVE: release-like configurations use ${_mbelib_native_flags_text} (${_mbelib_native_reason})"
        )
    else()
        message(
            STATUS
            "MBELIB_ENABLE_NATIVE: not tuning for the build host: ${_mbelib_native_reason}"
        )
    endif()
    unset(_mbelib_native_reason)
    unset(_mbelib_native_flags_text)
else()
    message(
        STATUS
        "MBELIB_ENABLE_NATIVE: off; the libraries keep the compiler's default target"
    )
endif()
# Genex true for the configurations that get the native flags.
set(_mbelib_native_config "0")
if(_mbelib_native_flags)
    set(_mbelib_native_config "${_mbelib_release_like_config}")
endif()
function(_mbelib_apply_native target)
    foreach(_mbelib_native_flag IN LISTS _mbelib_native_flags)
        string(
            REPLACE ","
            "$<COMMA>"
            _mbelib_native_flag
            "${_mbelib_native_flag}"
        )
        target_compile_options(
            ${target}
            PRIVATE "$<${_mbelib_native_config}:${_mbelib_native_flag}>"
        )
    endforeach()
endfunction()
_mbelib_apply_native(mbe-shared)
_mbelib_apply_native(mbe-static)

# Fast-math (optional). clang-cl only accepts GNU-style options behind /clang:.
set(_mbelib_fast_math_options "")
set(_mbelib_ieee_math_options "")
if(
    CMAKE_C_COMPILER_ID STREQUAL "Clang"
    AND CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC"
)
    set(_mbelib_fast_math_options /clang:-ffast-math /clang:-fno-math-errno)
    set(_mbelib_ieee_math_options /clang:-fno-fast-math /clang:-fno-math-errno)
elseif(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    set(_mbelib_fast_math_options -ffast-math -fno-math-errno)
    set(_mbelib_ieee_math_options -fno-fast-math -fno-math-errno)
elseif(MSVC)
    set(_mbelib_fast_math_options /fp:fast)
    set(_mbelib_ieee_math_options /fp:precise)
endif()
if(MBELIB_ENABLE_FAST_MATH AND _mbelib_fast_math_options)
    target_compile_options(mbe-shared PRIVATE ${_mbelib_fast_math_options})
    target_compile_options(mbe-static PRIVATE ${_mbelib_fast_math_options})
endif()

# PCM conversion requires IEEE NaN/Inf handling regardless of the fast-math
# option or externally supplied C flags. Per-source COMPILE_OPTIONS follow the
# target options on Makefile/Ninja command lines; the Visual Studio generator
# maps /fp:precise to the per-file floating-point model instead. The unit also
# refuses to compile under fast-math (see src/core/mbe_pcm_convert.c).
# The SIMD voiced kernels (and their AVX2 copy) seed their oscillators from
# cos and sin of caller-supplied phases of any size; a fast-math math library
# need not reduce huge arguments exactly (MSVC /fp:fast returned values off the
# unit circle), and the kernels themselves are explicit vector code. Both units
# refuse to compile under fast-math.
if(_mbelib_ieee_math_options)
    set_source_files_properties(
        src/core/mbe_pcm_convert.c
        src/core/mbe_voiced.c
        src/core/mbe_voiced_avx2.c
        PROPERTIES COMPILE_OPTIONS "${_mbelib_ieee_math_options}"
    )
endif()
unset(_mbelib_fast_math_options)
unset(_mbelib_ieee_math_options)

# Run-time AVX2/FMA dispatch (src/internal/mbe_cpu.h). A SIMD build for
# baseline x86-64 compiles the kernels whose output depends on the vector
# width, and the integer kernels that gain from AVX2, a second time with
# AVX2/FMA code generation and picks one per call from a cached CPU probe. A
# configuration whose target already has AVX2 (CMAKE_C_FLAGS or a
# CMAKE_C_FLAGS_<CONFIG>) calls its kernels directly: the preprocessor
# (mbe_cpu.h) decides per translation unit, so the copies are built for every
# configuration without native tuning, and those that do not need them leave
# them unused. Natively tuned configurations compile them empty. The flags are
# appended after the IEEE options above.
# 32-bit x86 and universal macOS builds do not dispatch. MinGW and Cygwin GCC
# do not align the stack for 32-byte spills (GCC PR 54412), so they dispatch
# only when the assembler can turn aligned vector moves into unaligned ones.
set(MBELIB_AVX2_DISPATCH_SOURCES
    src/core/mbe_voiced_avx2.c
    src/core/mbe_tone_fit_avx2.c
    src/ecc/ecc_soft_search_avx2.c
    src/tetra/tetra_acelp_avx2.c
)
set(MBELIB_AVX2_DISPATCH_FLAGS "")
set(_mbelib_avx2_reason "")
list(LENGTH CMAKE_OSX_ARCHITECTURES _mbelib_osx_arch_count)
if(NOT MBELIB_ENABLE_SIMD)
    set(_mbelib_avx2_reason "MBELIB_ENABLE_SIMD is off")
elseif(_mbelib_osx_arch_count GREATER 1)
    set(_mbelib_avx2_reason "several CMAKE_OSX_ARCHITECTURES")
else()
    _mbelib_check_c_source_compiles_static(
      [=[
      #if !(defined(__x86_64__) || defined(_M_X64)) || defined(_M_ARM64EC)
      #error not an x86-64 target
      #endif
      int mbelib_target_x86_64;
    ]=]
      MBELIB_TARGET_X86_64
    )
    if(NOT MBELIB_TARGET_X86_64)
        set(_mbelib_avx2_reason "not an x86-64 target")
    else()
        set(_mbelib_avx2_candidate "")
        if(MSVC)
            set(_mbelib_avx2_candidate /arch:AVX2)
        elseif(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
            set(_mbelib_avx2_candidate -mavx2 -mfma)
            if((WIN32 OR CYGWIN) AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
                list(
                    APPEND _mbelib_avx2_candidate
                    -Wa,-muse-unaligned-vector-move
                )
            endif()
        endif()
        if(_mbelib_avx2_candidate)
            list(JOIN _mbelib_avx2_candidate " " _mbelib_avx2_required_flags)
            set(_mbelib_saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
            set(CMAKE_REQUIRED_FLAGS
                "${CMAKE_REQUIRED_FLAGS} ${_mbelib_avx2_required_flags}"
            )
            _mbelib_check_c_source_compiles_static(
              [=[
              #include <immintrin.h>
              int mbelib_avx2_probe(const int* table, const float* x) {
                  __m256i g = _mm256_i32gather_epi32(table, _mm256_set1_epi32(1), 4);
                  __m256 f = _mm256_fmadd_ps(_mm256_loadu_ps(x), _mm256_set1_ps(2.0f), _mm256_set1_ps(1.0f));
                  return _mm256_movemask_epi8(g) + _mm256_movemask_ps(f);
              }
            ]=]
              MBELIB_C_HAS_AVX2_DISPATCH_FLAGS
            )
            set(CMAKE_REQUIRED_FLAGS "${_mbelib_saved_required_flags}")
            unset(_mbelib_saved_required_flags)
            unset(_mbelib_avx2_required_flags)
            if(MBELIB_C_HAS_AVX2_DISPATCH_FLAGS)
                set(MBELIB_AVX2_DISPATCH_FLAGS ${_mbelib_avx2_candidate})
            else()
                set(_mbelib_avx2_reason
                    "the compiler rejects ${_mbelib_avx2_candidate}"
                )
            endif()
        else()
            set(_mbelib_avx2_reason "no AVX2 option for ${CMAKE_C_COMPILER_ID}")
        endif()
        unset(_mbelib_avx2_candidate)
    endif()
    # Probe again on every configure: CMAKE_C_FLAGS may have gained or lost a
    # target ISA in this build directory.
    unset(MBELIB_C_HAS_AVX2_DISPATCH_FLAGS CACHE)
endif()
unset(_mbelib_osx_arch_count)
# Genex true for the configurations that build and call the AVX2 copies.
set(_mbelib_avx2_dispatch_config "0")
if(MBELIB_AVX2_DISPATCH_FLAGS)
    set(_mbelib_avx2_dispatch_config "$<NOT:${_mbelib_native_config}>")
    foreach(_target mbe-shared mbe-static)
        target_compile_definitions(
            ${_target}
            PRIVATE "$<${_mbelib_avx2_dispatch_config}:MBELIB_AVX2_KERNELS=1>"
        )
    endforeach()
    set(_mbelib_avx2_options "")
    foreach(_mbelib_flag IN LISTS MBELIB_AVX2_DISPATCH_FLAGS)
        string(REPLACE "," "$<COMMA>" _mbelib_flag "${_mbelib_flag}")
        list(
            APPEND _mbelib_avx2_options
            "$<${_mbelib_avx2_dispatch_config}:${_mbelib_flag}>"
        )
    endforeach()
    set_property(
        SOURCE ${MBELIB_AVX2_DISPATCH_SOURCES}
        APPEND
        PROPERTY COMPILE_OPTIONS ${_mbelib_avx2_options}
    )
    unset(_mbelib_avx2_options)
    unset(_mbelib_flag)
    list(JOIN MBELIB_AVX2_DISPATCH_FLAGS " " _mbelib_avx2_flags_text)
    if(_mbelib_native_flags)
        message(
            STATUS
            "AVX2 dispatch: ${_mbelib_avx2_flags_text} kernels chosen at run time in configurations without native tuning or an AVX2 target"
        )
    else()
        message(
            STATUS
            "AVX2 dispatch: ${_mbelib_avx2_flags_text} kernels chosen at run time in configurations without an AVX2 target"
        )
    endif()
    unset(_mbelib_avx2_flags_text)
else()
    message(STATUS "AVX2 dispatch: off (${_mbelib_avx2_reason})")
endif()
unset(_mbelib_avx2_reason)

# Sanitizers (Debug configuration only, for single- and multi-config
# generators). Collect the flags once and apply them to the libraries and to
# the executables (tests/examples) so the sanitizer runtimes are linked into
# the final executables.
set(_mbelib_san_cflags "")
set(_mbelib_san_ldflags "")
if(MBELIB_ENABLE_ASAN AND (CMAKE_C_COMPILER_ID MATCHES "Clang|GNU"))
    list(
        APPEND _mbelib_san_cflags
        "$<$<CONFIG:Debug>:-fsanitize=address>"
        "$<$<CONFIG:Debug>:-fno-omit-frame-pointer>"
    )
    # For Clang, enable stronger detection of stack use-after-scope
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        list(
            APPEND _mbelib_san_cflags
            "$<$<CONFIG:Debug>:-fsanitize-address-use-after-scope>"
        )
    endif()
    list(APPEND _mbelib_san_ldflags "$<$<CONFIG:Debug>:-fsanitize=address>")
endif()
if(MBELIB_ENABLE_UBSAN AND (CMAKE_C_COMPILER_ID MATCHES "Clang|GNU"))
    list(
        APPEND _mbelib_san_cflags
        "$<$<CONFIG:Debug>:-fsanitize=undefined>"
        "$<$<CONFIG:Debug>:-fno-sanitize-recover=undefined>"
    )
    list(APPEND _mbelib_san_ldflags "$<$<CONFIG:Debug>:-fsanitize=undefined>")
endif()

if(_mbelib_san_cflags OR _mbelib_san_ldflags)
    # Apply to core libraries first
    foreach(tgt IN ITEMS mbe-shared mbe-static)
        if(_mbelib_san_cflags)
            target_compile_options(${tgt} PRIVATE ${_mbelib_san_cflags})
        endif()
        if(_mbelib_san_ldflags)
            # Add sanitizer link flags before other options for correct ordering
            target_link_options(${tgt} BEFORE PRIVATE ${_mbelib_san_ldflags})
        endif()
    endforeach()
    # Stash the flags for later targets (tests/examples) via directory properties
    set_directory_properties(
        PROPERTIES
            MBELIB_SANITIZER_CFLAGS "${_mbelib_san_cflags}"
            MBELIB_SANITIZER_LDFLAGS "${_mbelib_san_ldflags}"
    )
endif()

include(GNUInstallDirs)

# On Windows, static and import libraries both use the .lib extension.
# Give the static library a distinct name to avoid collisions during parallel builds.
if(WIN32)
    set(_MBELIB_STATIC_OUTPUT_NAME mbe-neo-static)
else()
    set(_MBELIB_STATIC_OUTPUT_NAME mbe-neo)
endif()

set_target_properties(
    mbe-static
    PROPERTIES
        OUTPUT_NAME ${_MBELIB_STATIC_OUTPUT_NAME}
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
)
set_target_properties(
    mbe-shared
    PROPERTIES
        OUTPUT_NAME mbe-neo
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        C_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN YES
        # On macOS, use @rpath so consumers control search paths at runtime
        MACOSX_RPATH ON
        INSTALL_NAME_DIR "@rpath"
)

# Optional LTO/IPO for the Release and RelWithDebInfo configurations, set per
# configuration so multi-config generators get it too.
set(_mbelib_ipo_enabled OFF)
if(MBELIB_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_msg)
    if(_ipo_ok)
        set_target_properties(
            mbe-static
            mbe-shared
            PROPERTIES
                INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE
                INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE
        )
        set(_mbelib_ipo_enabled ON)
    else()
        message(STATUS "IPO/LTO not supported: ${_ipo_msg}")
    endif()
endif()

# Generate pkg-config metadata
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libmbe-neo.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/libmbe-neo.pc"
    @ONLY
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/libmbe-neo.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

install(
    TARGETS mbe-static mbe-shared
    EXPORT mbe-neoTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install public headers (source) and generated headers
install(
    FILES "${PROJECT_SOURCE_DIR}/include/mbelib-neo/mbelib.h"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mbelib-neo
)
install(
    FILES "${MBELIB_GENERATED_INCLUDE_DIR}/mbelib-neo/version.h"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mbelib-neo
)
install(
    FILES
        "${PROJECT_SOURCE_DIR}/LICENSE"
        "${PROJECT_SOURCE_DIR}/LICENSES/GPL-2.0-or-later.txt"
        "${PROJECT_SOURCE_DIR}/LICENSES/ISC.txt"
        "${PROJECT_SOURCE_DIR}/LICENSES/LicenseRef-PFFFT.txt"
        "${PROJECT_SOURCE_DIR}/LICENSES/LicenseRef-ETSI-EN-300-395-2-data.txt"
        "${PROJECT_SOURCE_DIR}/docs/upstream-mbelib-attribution.md"
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
)

install(
    FILES "${MBELIB_PFFFT_SOURCE_DIR}/LICENSE.txt"
    DESTINATION ${CMAKE_INSTALL_DOCDIR}/pffft
)
install(FILES "${PFFFT_NOTICES}" DESTINATION ${CMAKE_INSTALL_DOCDIR}/pffft)

install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/legal/"
    DESTINATION ${CMAKE_INSTALL_DOCDIR}/legal
    FILES_MATCHING
    PATTERN "*.md"
    PATTERN "*.json"
)

# Exported targets and CMake package config
install(
    EXPORT mbe-neoTargets
    NAMESPACE mbe_neo::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mbe-neo
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/mbe-neoConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/mbe-neoConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mbe-neo
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/mbe-neoConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/mbe-neoConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/mbe-neoConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mbe-neo
)

# Choose link target for executables/tests: prefer shared; on Windows use static to avoid DLL path issues
set(MBELIB_EXEC_LINK_TGT mbe-shared)
if(WIN32)
    set(MBELIB_EXEC_LINK_TGT mbe-static)
endif()

# Examples
if(MBELIB_BUILD_EXAMPLES)
    add_executable(example_print_version examples/print_version.c)
    target_link_libraries(example_print_version PRIVATE ${MBELIB_EXEC_LINK_TGT})
    _mbelib_apply_hardening(example_print_version)
    target_include_directories(
        example_print_version
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    # Propagate sanitizer flags to the example executable
    get_directory_property(_SAN_CFLAGS MBELIB_SANITIZER_CFLAGS)
    get_directory_property(_SAN_LDFLAGS MBELIB_SANITIZER_LDFLAGS)
    if(_SAN_CFLAGS)
        target_compile_options(example_print_version PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            example_print_version
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()

    add_executable(example_dstar_encode examples/dstar_encode.c)
    target_link_libraries(example_dstar_encode PRIVATE ${MBELIB_EXEC_LINK_TGT})
    _mbelib_apply_hardening(example_dstar_encode)
    target_include_directories(
        example_dstar_encode
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    set_target_properties(
        example_dstar_encode
        PROPERTIES OUTPUT_NAME dstar_encode
    )
    if(_SAN_CFLAGS)
        target_compile_options(example_dstar_encode PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(example_dstar_encode BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()

    add_executable(example_dstar_decode examples/dstar_decode.c)
    target_link_libraries(example_dstar_decode PRIVATE ${MBELIB_EXEC_LINK_TGT})
    _mbelib_apply_hardening(example_dstar_decode)
    target_include_directories(
        example_dstar_decode
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    set_target_properties(
        example_dstar_decode
        PROPERTIES OUTPUT_NAME dstar_decode
    )
    if(_SAN_CFLAGS)
        target_compile_options(example_dstar_decode PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(example_dstar_decode BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
endif()

# Tests
if(MBELIB_BUILD_TESTS)
    include(CTest)
    enable_testing()

    add_executable(test_assertions tests/test_assertions.c)
    add_test(NAME test_assertions COMMAND test_assertions)

    # Golden data generator tool (not a test): helpful for refreshing hashes
    add_executable(gen_golden tools/gen_golden.c)
    # The kernel-ISA selection is private; link the static library.
    target_link_libraries(gen_golden PRIVATE mbe-static)
    _mbelib_apply_hardening(gen_golden)
    target_include_directories(
        gen_golden
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
            "${PROJECT_SOURCE_DIR}/tests"
    )
    add_executable(test_api tests/test_api.c)
    target_link_libraries(test_api PRIVATE ${MBELIB_EXEC_LINK_TGT})
    target_include_directories(test_api PRIVATE "${PROJECT_SOURCE_DIR}/include")
    # Propagate sanitizer flags to test executables
    get_directory_property(_SAN_CFLAGS MBELIB_SANITIZER_CFLAGS)
    get_directory_property(_SAN_LDFLAGS MBELIB_SANITIZER_LDFLAGS)
    if(_SAN_CFLAGS)
        target_compile_options(test_api PRIVATE ${_SAN_CFLAGS})
        target_compile_options(gen_golden PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_api BEFORE PRIVATE ${_SAN_LDFLAGS})
        target_link_options(gen_golden BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_test(NAME test_api COMMAND test_api)

    # test_api links the static library on Windows; this one calls the TETRA exports through
    # the shared library everywhere. The DLL is built next to the test, so Windows finds it.
    add_executable(test_tetra_shared_api tests/test_tetra_shared_api.c)
    target_link_libraries(test_tetra_shared_api PRIVATE mbe-shared)
    target_include_directories(
        test_tetra_shared_api
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_tetra_shared_api PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            test_tetra_shared_api
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()
    add_test(NAME test_tetra_shared_api COMMAND test_tetra_shared_api)

    add_executable(test_ecc tests/test_ecc.c)
    target_link_libraries(test_ecc PRIVATE ${MBELIB_EXEC_LINK_TGT})
    target_include_directories(
        test_ecc
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_ecc PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_ecc BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_test(NAME test_ecc COMMAND test_ecc)

    add_executable(test_speech_analysis tests/test_speech_analysis.c)
    # The analysis is private; link the static library for its symbols.
    target_link_libraries(test_speech_analysis PRIVATE mbe-static)
    target_include_directories(
        test_speech_analysis
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_speech_analysis PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_speech_analysis BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_test(NAME test_speech_analysis COMMAND test_speech_analysis)

    add_executable(test_ambe2400_encoder tests/test_ambe2400_encoder.c)
    # The private Golay encoder is hidden in the shared library.
    target_link_libraries(test_ambe2400_encoder PRIVATE mbe-static)
    if(NOTONES)
        target_compile_definitions(
            test_ambe2400_encoder
            PRIVATE MBELIB_TEST_NOTONES=1
        )
    endif()
    target_include_directories(
        test_ambe2400_encoder
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_ambe2400_encoder PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            test_ambe2400_encoder
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()
    add_test(NAME test_ambe2400_encoder COMMAND test_ambe2400_encoder)
    # The floor-fit overflow-flag check needs an unoptimized IEEE build of the
    # library and the test, where every floating-point operation runs as
    # written. Elsewhere the flags carry no guarantee: fast-math, optimization
    # without FENV_ACCESS (which GCC does not implement), or MSVC, whose
    # fenv_access is off by default.
    if(NOT MBELIB_ENABLE_FAST_MATH AND NOT MSVC)
        target_compile_definitions(
            test_ambe2400_encoder
            PRIVATE $<$<CONFIG:Debug>:MBELIB_TEST_FENV_OVERFLOW=1>
        )
    endif()
    # The AMBE+2 and IMBE encoder tests drive private quantizer steps and the
    # Hamming encoder, which the shared library does not export.
    foreach(_enc_test test_ambe2450_encoder test_imbe7200_encoder)
        add_executable(${_enc_test} tests/${_enc_test}.c)
        target_link_libraries(${_enc_test} PRIVATE mbe-static)
        if(NOTONES)
            target_compile_definitions(
                ${_enc_test}
                PRIVATE MBELIB_TEST_NOTONES=1
            )
        endif()
        target_include_directories(
            ${_enc_test}
            PRIVATE
                "${PROJECT_SOURCE_DIR}/include"
                "${PROJECT_SOURCE_DIR}/src/internal"
        )
        if(_SAN_CFLAGS)
            target_compile_options(${_enc_test} PRIVATE ${_SAN_CFLAGS})
        endif()
        if(_SAN_LDFLAGS)
            target_link_options(${_enc_test} BEFORE PRIVATE ${_SAN_LDFLAGS})
        endif()
        add_test(NAME ${_enc_test} COMMAND ${_enc_test})
    endforeach()
    # Exercise every partial plan allocation and retry with GNU link wrapping,
    # for each encoder (they share the analysis front end).
    # LTO can resolve allocator calls before the linker can wrap them.
    if(
        CMAKE_SYSTEM_NAME STREQUAL "Linux"
        AND CMAKE_C_COMPILER_ID MATCHES "Clang|GNU"
        AND NOT (
            MBELIB_ENABLE_LTO
            AND (
                _mbelib_multi_config
                OR CMAKE_BUILD_TYPE MATCHES "Release|RelWithDebInfo"
            )
        )
    )
        add_executable(test_ambe2400_encoder_oom tests/test_ambe2400_encoder.c)
        target_link_libraries(test_ambe2400_encoder_oom PRIVATE mbe-static)
        target_include_directories(
            test_ambe2400_encoder_oom
            PRIVATE
                "${PROJECT_SOURCE_DIR}/include"
                "${PROJECT_SOURCE_DIR}/src/internal"
        )
        target_compile_definitions(
            test_ambe2400_encoder_oom
            PRIVATE MBE_ENCODER_TEST_OOM=1
        )
        target_link_options(
            test_ambe2400_encoder_oom
            PRIVATE
                "-Wl,--wrap=pffft_aligned_malloc"
                "-Wl,--wrap=pffft_aligned_free"
                "-Wl,--wrap=calloc"
                "-Wl,--wrap=malloc"
                "-Wl,--wrap=free"
        )
        _mbelib_apply_hardening(test_ambe2400_encoder_oom)
        if(_SAN_CFLAGS)
            target_compile_options(
                test_ambe2400_encoder_oom
                PRIVATE ${_SAN_CFLAGS}
            )
        endif()
        if(_SAN_LDFLAGS)
            target_link_options(
                test_ambe2400_encoder_oom
                BEFORE
                PRIVATE ${_SAN_LDFLAGS}
            )
        endif()
        foreach(allocation RANGE 1 7)
            add_test(
                NAME test_ambe2400_encoder_oom_${allocation}
                COMMAND test_ambe2400_encoder_oom ${allocation}
            )
            add_test(
                NAME test_ambe2450_encoder_oom_${allocation}
                COMMAND test_ambe2400_encoder_oom ${allocation} ambe2450
            )
            add_test(
                NAME test_imbe7200_encoder_oom_${allocation}
                COMMAND test_ambe2400_encoder_oom ${allocation} imbe4400
            )
        endforeach()
    endif()

    add_executable(test_noise_determinism tests/test_noise_determinism.c)
    target_link_libraries(
        test_noise_determinism
        PRIVATE ${MBELIB_EXEC_LINK_TGT}
    )
    target_include_directories(
        test_noise_determinism
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_noise_determinism PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            test_noise_determinism
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()
    add_test(NAME test_noise_determinism COMMAND test_noise_determinism)

    add_executable(test_params tests/test_params.c)
    # The WOLA invariant exercises private DSP symbols without exporting them.
    target_link_libraries(test_params PRIVATE mbe-static)
    target_include_directories(
        test_params
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_params PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_params BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_test(NAME test_params COMMAND test_params)

    # Private-stage tests include direct equation oracles and fixed-point operators.
    foreach(
        _mbelib_private_test
        IN
        ITEMS
            test_fixed_point
            test_tetra_decoder
            test_tetra_simd
            test_dsp_equations
            test_adaptive_equations
    )
        add_executable(${_mbelib_private_test} tests/${_mbelib_private_test}.c)
        target_link_libraries(${_mbelib_private_test} PRIVATE mbe-static)
        target_include_directories(
            ${_mbelib_private_test}
            PRIVATE
                "${PROJECT_SOURCE_DIR}/include"
                "${PROJECT_SOURCE_DIR}/src/internal"
        )
        if(_SAN_CFLAGS)
            target_compile_options(
                ${_mbelib_private_test}
                PRIVATE ${_SAN_CFLAGS}
            )
        endif()
        if(_SAN_LDFLAGS)
            target_link_options(
                ${_mbelib_private_test}
                BEFORE
                PRIVATE ${_SAN_LDFLAGS}
            )
        endif()
        add_test(NAME ${_mbelib_private_test} COMMAND ${_mbelib_private_test})
    endforeach()
    unset(_mbelib_private_test)

    add_executable(
        test_frame_paths
        tests/test_frame_paths.c
        tools/quality/mbe_quality_frames.c
        src/ecc/ecc_const.c
    )
    target_link_libraries(test_frame_paths PRIVATE ${MBELIB_EXEC_LINK_TGT})
    target_include_directories(
        test_frame_paths
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
            "${PROJECT_SOURCE_DIR}/tools/quality"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_frame_paths PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_frame_paths BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_test(NAME test_frame_paths COMMAND test_frame_paths)

    add_executable(
        test_quality_dvsi
        tests/test_quality_dvsi.c
        tools/quality/mbe_quality_dvsi.c
        tools/quality/mbe_quality_frames.c
        src/ecc/ecc_const.c
    )
    target_link_libraries(test_quality_dvsi PRIVATE ${MBELIB_EXEC_LINK_TGT})
    target_include_directories(
        test_quality_dvsi
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
            "${PROJECT_SOURCE_DIR}/tools/quality"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_quality_dvsi PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_quality_dvsi BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_test(NAME test_quality_dvsi COMMAND test_quality_dvsi)

    add_executable(test_input_validation tests/test_input_validation.c)
    target_link_libraries(test_input_validation PRIVATE ${MBELIB_EXEC_LINK_TGT})
    target_include_directories(
        test_input_validation
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_input_validation PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            test_input_validation
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()
    add_test(NAME test_input_validation COMMAND test_input_validation)

    add_executable(test_floattoshort_parity tests/test_floattoshort_parity.c)
    target_link_libraries(
        test_floattoshort_parity
        PRIVATE ${MBELIB_EXEC_LINK_TGT}
    )
    target_include_directories(
        test_floattoshort_parity
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_floattoshort_parity PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            test_floattoshort_parity
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()
    add_test(NAME test_floattoshort_parity COMMAND test_floattoshort_parity)

    add_executable(test_ambe2450_frame_types tests/test_ambe2450_frame_types.c)
    target_link_libraries(
        test_ambe2450_frame_types
        PRIVATE ${MBELIB_EXEC_LINK_TGT}
    )
    target_include_directories(
        test_ambe2450_frame_types
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_ambe2450_frame_types PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            test_ambe2450_frame_types
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()
    add_test(NAME test_ambe2450_frame_types COMMAND test_ambe2450_frame_types)

    add_executable(test_tones tests/test_tones.c)
    target_link_libraries(test_tones PRIVATE ${MBELIB_EXEC_LINK_TGT})
    target_include_directories(
        test_tones
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(NOTONES)
        target_compile_definitions(test_tones PRIVATE MBELIB_TEST_NOTONES=1)
    endif()
    if(_SAN_CFLAGS)
        target_compile_options(test_tones PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_tones BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_test(NAME test_tones COMMAND test_tones)

    add_executable(test_golden_pcm tests/test_golden_pcm.c)
    # The hashes follow the library's kernel ISA, which the private mbe_cpu
    # interface reports and selects; link the static library for it.
    target_link_libraries(test_golden_pcm PRIVATE mbe-static)
    target_include_directories(
        test_golden_pcm
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
    )
    # In Debug, enforce strict float and int16 golden hashes on x86 (single- and multi-config)
    target_compile_definitions(
        test_golden_pcm
        PRIVATE
            $<$<CONFIG:Debug>:MBELIB_TEST_STRICT_FLOAT=1;MBELIB_TEST_STRICT_INT16=1>
    )
    if(_SAN_CFLAGS)
        target_compile_options(test_golden_pcm PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_golden_pcm BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    _mbelib_apply_native(test_golden_pcm)
    # test_golden_pcm runs the CPU probe's choice; the _sse2 and _avx2 variants
    # force either kernel ISA and skip where the build or CPU cannot run it
    # (MBELIB_TEST_REQUIRE_AVX2 turns the AVX2 skip into a failure).
    set(_mbelib_golden_avx2_args avx2)
    if(MBELIB_TEST_REQUIRE_AVX2)
        list(APPEND _mbelib_golden_avx2_args --require)
    endif()
    add_test(NAME test_golden_pcm COMMAND test_golden_pcm auto)
    add_test(NAME test_golden_pcm_sse2 COMMAND test_golden_pcm sse2)
    add_test(
        NAME test_golden_pcm_avx2
        COMMAND test_golden_pcm ${_mbelib_golden_avx2_args}
    )
    set_tests_properties(
        test_golden_pcm
        test_golden_pcm_sse2
        test_golden_pcm_avx2
        PROPERTIES SKIP_RETURN_CODE 77
    )
    unset(_mbelib_golden_avx2_args)

    add_executable(
        test_simd_arch_detection_x64
        tests/test_simd_arch_detection.c
    )
    target_include_directories(
        test_simd_arch_detection_x64
        PRIVATE "${PROJECT_SOURCE_DIR}/src/internal"
    )
    target_compile_definitions(
        test_simd_arch_detection_x64
        PRIVATE MBE_TEST_SCENARIO_X64=1
    )
    if(_SAN_CFLAGS)
        target_compile_options(
            test_simd_arch_detection_x64
            PRIVATE ${_SAN_CFLAGS}
        )
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            test_simd_arch_detection_x64
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()
    add_test(
        NAME test_simd_arch_detection_x64
        COMMAND test_simd_arch_detection_x64
    )

    add_executable(
        test_simd_arch_detection_arm64ec
        tests/test_simd_arch_detection.c
    )
    target_include_directories(
        test_simd_arch_detection_arm64ec
        PRIVATE "${PROJECT_SOURCE_DIR}/src/internal"
    )
    target_compile_definitions(
        test_simd_arch_detection_arm64ec
        PRIVATE MBE_TEST_SCENARIO_ARM64EC=1
    )
    if(_SAN_CFLAGS)
        target_compile_options(
            test_simd_arch_detection_arm64ec
            PRIVATE ${_SAN_CFLAGS}
        )
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            test_simd_arch_detection_arm64ec
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()
    add_test(
        NAME test_simd_arch_detection_arm64ec
        COMMAND test_simd_arch_detection_arm64ec
    )

    add_executable(
        test_simd_arch_detection_x86_32_msvc_sse2
        tests/test_simd_arch_detection.c
    )
    target_include_directories(
        test_simd_arch_detection_x86_32_msvc_sse2
        PRIVATE "${PROJECT_SOURCE_DIR}/src/internal"
    )
    target_compile_definitions(
        test_simd_arch_detection_x86_32_msvc_sse2
        PRIVATE MBE_TEST_SCENARIO_X86_32_MSVC_SSE2=1 _M_IX86_FP=2
    )
    if(_SAN_CFLAGS)
        target_compile_options(
            test_simd_arch_detection_x86_32_msvc_sse2
            PRIVATE ${_SAN_CFLAGS}
        )
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(
            test_simd_arch_detection_x86_32_msvc_sse2
            BEFORE
            PRIVATE ${_SAN_LDFLAGS}
        )
    endif()
    add_test(
        NAME test_simd_arch_detection_x86_32_msvc_sse2
        COMMAND test_simd_arch_detection_x86_32_msvc_sse2
    )

    if(MBELIB_BUILD_EXAMPLES AND NOT CMAKE_CROSSCOMPILING)
        find_package(Python3 3.7 QUIET COMPONENTS Interpreter)
        if(Python3_Interpreter_FOUND)
            add_test(
                NAME test_dstar_examples
                COMMAND
                    ${Python3_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}/tests/test_dstar_examples.py
                    $<TARGET_FILE:example_dstar_encode>
                    $<TARGET_FILE:example_dstar_decode>
            )
        endif()
    endif()

    if(MBELIB_BUILD_TOOLS)
        add_executable(
            test_quality_fs
            tests/test_quality_fs.c
            tools/quality/mbe_quality_fs.c
        )
        target_include_directories(
            test_quality_fs
            PRIVATE "${PROJECT_SOURCE_DIR}/tools/quality"
        )
        if(WIN32)
            if(MSVC)
                target_compile_definitions(
                    test_quality_fs
                    PRIVATE _CRT_SECURE_NO_WARNINGS
                )
            endif()
        else()
            target_compile_definitions(
                test_quality_fs
                PRIVATE _XOPEN_SOURCE=700
            )
            # Darwin hides mkdtemp() behind its full C level under strict X/Open.
            if(APPLE)
                target_compile_definitions(
                    test_quality_fs
                    PRIVATE _DARWIN_C_SOURCE
                )
            endif()
        endif()
        if(_SAN_CFLAGS)
            target_compile_options(test_quality_fs PRIVATE ${_SAN_CFLAGS})
        endif()
        if(_SAN_LDFLAGS)
            target_link_options(test_quality_fs BEFORE PRIVATE ${_SAN_LDFLAGS})
        endif()
        add_test(NAME test_quality_fs COMMAND test_quality_fs)
        if(NOT CMAKE_CROSSCOMPILING)
            find_package(Python3 3.7 QUIET COMPONENTS Interpreter)
            if(Python3_Interpreter_FOUND)
                add_test(
                    NAME test_quality_tools
                    COMMAND
                        ${Python3_EXECUTABLE}
                        ${PROJECT_SOURCE_DIR}/tests/test_quality_tools.py
                        $<TARGET_FILE:mbe_quality_eval>
                        $<TARGET_FILE:mbe_quality_reframe>
                )
                add_test(
                    NAME test_dvsi_scoreboard
                    COMMAND
                        ${Python3_EXECUTABLE}
                        ${PROJECT_SOURCE_DIR}/tests/test_dvsi_scoreboard.py
                )
                # The pitch measurements need NumPy; without it the test reports a skip.
                set_tests_properties(
                    test_dvsi_scoreboard
                    PROPERTIES SKIP_RETURN_CODE 77
                )
                # The A/B runner locks its run directory with fcntl, so it is POSIX-only, and
                # uses assignment expressions, so it needs Python 3.8.
                if(NOT WIN32 AND Python3_VERSION VERSION_GREATER_EQUAL 3.8)
                    add_test(
                        NAME test_quality_ab_finalize
                        COMMAND
                            ${Python3_EXECUTABLE}
                            ${PROJECT_SOURCE_DIR}/tests/test_quality_ab_finalize.py
                            ${PROJECT_SOURCE_DIR}/tools/quality
                    )
                endif()
            endif()
        endif()
    endif()

    find_package(Python3 3.8 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(
            NAME test_ci_scheduling
            COMMAND
                ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tests/test_ci_scheduling.py
        )
        add_test(
            NAME test_legal_release
            COMMAND
                ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tests/test_legal_release.py
        )
        add_test(
            NAME test_release_assets
            COMMAND
                ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tests/test_release_assets.py
        )
        add_test(
            NAME test_pffft_metadata
            COMMAND
                ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tests/test_pffft_metadata.py
        )
        add_test(
            NAME test_tetra_conformance_tools
            COMMAND
                ${Python3_EXECUTABLE} -I
                ${PROJECT_SOURCE_DIR}/tests/test_tetra_conformance_tools.py
        )
        add_test(
            NAME test_tetra_tables
            COMMAND
                ${Python3_EXECUTABLE} -I
                ${PROJECT_SOURCE_DIR}/tests/test_tetra_tables.py
        )
        # The AVX2 units export only *_avx2 symbols (src/internal/mbe_cpu.h).
        if(
            CMAKE_SYSTEM_NAME STREQUAL "Linux"
            AND CMAKE_C_COMPILER_ID MATCHES "Clang|GNU"
            AND CMAKE_NM
            AND MBELIB_AVX2_DISPATCH_FLAGS
        )
            add_test(
                NAME test_avx2_symbols
                COMMAND
                    ${Python3_EXECUTABLE} -I
                    ${PROJECT_SOURCE_DIR}/tests/test_avx2_symbols.py ${CMAKE_NM}
                    $<TARGET_FILE:mbe-static>
                    "--expect-code=$<IF:${_mbelib_avx2_dispatch_config},1,0>"
                    "--lto=$<IF:$<AND:$<BOOL:${_mbelib_ipo_enabled}>,$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>>,1,0>"
            )
            set_tests_properties(
                test_avx2_symbols
                PROPERTIES SKIP_RETURN_CODE 77
            )
        endif()
        # Configure-only checks of the per-configuration SIMD, LTO, native and
        # dispatch settings, with Ninja and Ninja Multi-Config. MSVC needs its
        # developer environment for a Ninja configure, which a Visual Studio
        # generator build does not provide to its tests.
        if(
            NOT CMAKE_CROSSCOMPILING
            AND NOT (MSVC AND NOT CMAKE_GENERATOR MATCHES "Ninja")
        )
            if(CMAKE_GENERATOR MATCHES "Ninja")
                set(_mbelib_test_ninja "${CMAKE_MAKE_PROGRAM}")
            else()
                find_program(MBELIB_NINJA_EXECUTABLE NAMES ninja ninja-build)
                set(_mbelib_test_ninja "")
                if(MBELIB_NINJA_EXECUTABLE)
                    set(_mbelib_test_ninja "${MBELIB_NINJA_EXECUTABLE}")
                endif()
            endif()
            set(_mbelib_test_user_isa "")
            if(MSVC AND MBELIB_NATIVE_TARGET_X86)
                set(_mbelib_test_user_isa "/arch:AVX")
            elseif(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
                # A valid baseline target ISA for each architecture.
                if(MBELIB_NATIVE_TARGET_X86)
                    set(_mbelib_test_user_isa "-march=x86-64")
                elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
                    set(_mbelib_test_user_isa "-march=armv8-a")
                elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^armv7")
                    set(_mbelib_test_user_isa "-march=armv7-a")
                endif()
            endif()
            add_test(
                NAME test_build_configuration
                COMMAND
                    ${Python3_EXECUTABLE} -I
                    ${PROJECT_SOURCE_DIR}/tests/test_build_configuration.py
                    --cmake ${CMAKE_COMMAND} --source ${PROJECT_SOURCE_DIR}
                    --pffft-source ${MBELIB_PFFFT_SOURCE_DIR} --c-compiler
                    ${CMAKE_C_COMPILER}
                    "--c-compiler-arg1=${CMAKE_C_COMPILER_ARG1}"
                    "--ninja=${_mbelib_test_ninja}"
                    "--user-isa=${_mbelib_test_user_isa}"
            )
            set_tests_properties(
                test_build_configuration
                PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900
            )
            unset(_mbelib_test_ninja)
            unset(_mbelib_test_user_isa)
        endif()
        if(MBELIB_BUILD_TOOLS AND NOT CMAKE_CROSSCOMPILING)
            add_test(
                NAME test_tetra_decode_tool
                COMMAND
                    ${Python3_EXECUTABLE} -I
                    ${PROJECT_SOURCE_DIR}/tests/test_tetra_conformance_tools.py
                    $<TARGET_FILE:mbe_tetra_decode>
            )
        endif()
    endif()

    # tools/ci_changed_files.sh decides which pull requests PR fuzzing runs on. It is
    # CI tooling for the Linux runners and needs bash 4.4 or later, so the check
    # registers on Linux hosts that have such a bash, git and Python.
    if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux" AND NOT CMAKE_CROSSCOMPILING)
        find_package(Python3 3.7 QUIET COMPONENTS Interpreter)
        find_package(Git QUIET)
        find_program(MBELIB_BASH_EXECUTABLE bash)
        set(_mbelib_bash_version "")
        if(MBELIB_BASH_EXECUTABLE)
            execute_process(
                COMMAND
                    ${MBELIB_BASH_EXECUTABLE} -c
                    [[echo "${BASH_VERSINFO[0]}.${BASH_VERSINFO[1]}"]]
                OUTPUT_VARIABLE _mbelib_bash_version
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
        endif()
        if(
            Python3_Interpreter_FOUND
            AND GIT_FOUND
            AND _mbelib_bash_version VERSION_GREATER_EQUAL 4.4
        )
            add_test(
                NAME test_ci_changed_files
                COMMAND
                    ${Python3_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}/tests/test_ci_changed_files.py
                    ${MBELIB_BASH_EXECUTABLE} ${GIT_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}/tools/ci_changed_files.sh
            )
            # .githooks/pre-push takes each pushed ref's files from
            # tools/push_changed_files.sh, which has the same requirements.
            add_test(
                NAME test_push_changed_files
                COMMAND
                    ${Python3_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}/tests/test_push_changed_files.py
                    ${MBELIB_BASH_EXECUTABLE} ${GIT_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}/tools/push_changed_files.sh
            )
            add_test(
                NAME test_pre_push_paths
                COMMAND
                    ${Python3_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}/tests/test_pre_push_paths.py
                    ${MBELIB_BASH_EXECUTABLE} ${GIT_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}
            )
        else()
            message(
                STATUS
                "test_ci_changed_files, test_push_changed_files and test_pre_push_paths not registered: needs Python 3, git and bash 4.4 or later"
            )
        endif()
        # tools/shell_lint.sh must fail when either of its tools does. The shell/workflow
        # lint job also runs this test, where both tools are installed.
        find_program(MBELIB_SHELLCHECK_EXECUTABLE shellcheck)
        find_program(MBELIB_SHFMT_EXECUTABLE shfmt)
        if(
            Python3_Interpreter_FOUND
            AND GIT_FOUND
            AND _mbelib_bash_version VERSION_GREATER_EQUAL 4.4
            AND MBELIB_SHELLCHECK_EXECUTABLE
            AND MBELIB_SHFMT_EXECUTABLE
        )
            add_test(
                NAME test_shell_lint
                COMMAND
                    ${Python3_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}/tests/test_shell_lint.py
                    ${PROJECT_SOURCE_DIR}/tools/shell_lint.sh
                    ${MBELIB_BASH_EXECUTABLE} ${GIT_EXECUTABLE}
                    ${MBELIB_SHELLCHECK_EXECUTABLE} ${MBELIB_SHFMT_EXECUTABLE}
            )
        endif()
        unset(_mbelib_bash_version)
    endif()

    add_executable(test_simd_numerics tests/test_simd_numerics.c)
    target_link_libraries(test_simd_numerics PRIVATE mbe-static)
    target_include_directories(
        test_simd_numerics
        PRIVATE "${PROJECT_SOURCE_DIR}/src/internal"
    )
    # Select the same private kernels as the library; a target the test sees
    # without vectors (e.g. 32-bit x86 without SSE2 flags) skips their checks.
    if(MBELIB_ENABLE_SIMD)
        target_compile_definitions(
            test_simd_numerics
            PRIVATE MBELIB_ENABLE_SIMD=1
        )
    endif()
    if(_SAN_CFLAGS)
        target_compile_options(test_simd_numerics PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_simd_numerics BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    _mbelib_apply_native(test_simd_numerics)
    add_test(NAME test_simd_numerics COMMAND test_simd_numerics)

    add_executable(test_cpu_dispatch tests/test_cpu_dispatch.c)
    # The kernel-ISA selection is private; link the static library for it.
    target_link_libraries(test_cpu_dispatch PRIVATE mbe-static)
    target_include_directories(
        test_cpu_dispatch
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
    )
    # Declare the voiced kernels as the library's SIMD build does.
    if(MBELIB_ENABLE_SIMD)
        target_compile_definitions(
            test_cpu_dispatch
            PRIVATE MBELIB_ENABLE_SIMD=1
        )
    endif()
    if(_SAN_CFLAGS)
        target_compile_options(test_cpu_dispatch PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(test_cpu_dispatch BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_test(NAME test_cpu_dispatch COMMAND test_cpu_dispatch)

    set(_MBELIB_TEST_TARGETS
        test_cpu_dispatch
        test_dsp_equations
        test_adaptive_equations
        test_assertions
        test_api
        test_tetra_shared_api
        test_ecc
        test_ambe2400_encoder
        test_ambe2450_encoder
        test_imbe7200_encoder
        test_speech_analysis
        test_simd_numerics
        test_noise_determinism
        test_params
        test_fixed_point
        test_tetra_decoder
        test_tetra_simd
        test_frame_paths
        test_quality_dvsi
        test_input_validation
        test_floattoshort_parity
        test_ambe2450_frame_types
        test_tones
        test_golden_pcm
        test_simd_arch_detection_x64
        test_simd_arch_detection_arm64ec
        test_simd_arch_detection_x86_32_msvc_sse2
    )
    if(MBELIB_BUILD_TOOLS)
        list(APPEND _MBELIB_TEST_TARGETS test_quality_fs)
    endif()
    foreach(_mbelib_test_target IN LISTS _MBELIB_TEST_TARGETS)
        _mbelib_apply_hardening(${_mbelib_test_target})
        if(MSVC)
            # /UNDEBUG conflicts with the release configuration's /DNDEBUG
            # and emits D9025. Undefine it before the test source is read.
            target_compile_options(
                ${_mbelib_test_target}
                PRIVATE "/FI${PROJECT_SOURCE_DIR}/tests/enable_assertions.h"
            )
        elseif(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
            target_compile_options(${_mbelib_test_target} PRIVATE -UNDEBUG)
        endif()
    endforeach()
    unset(_MBELIB_TEST_TARGETS)
    unset(_mbelib_test_target)
endif()

# Benchmarks (opt-in)
if(MBELIB_BUILD_BENCHMARKS)
    add_executable(bench_synth bench/bench_synth.c)
    target_link_libraries(bench_synth PRIVATE ${MBELIB_EXEC_LINK_TGT})
    _mbelib_apply_hardening(bench_synth)
    target_include_directories(
        bench_synth
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    # Propagate sanitizer flags if present
    get_directory_property(_SAN_CFLAGS MBELIB_SANITIZER_CFLAGS)
    get_directory_property(_SAN_LDFLAGS MBELIB_SANITIZER_LDFLAGS)
    if(_SAN_CFLAGS)
        target_compile_options(bench_synth PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(bench_synth BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_executable(bench_encode bench/bench_encode.c)
    target_link_libraries(bench_encode PRIVATE ${MBELIB_EXEC_LINK_TGT})
    _mbelib_apply_hardening(bench_encode)
    target_include_directories(
        bench_encode
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(bench_encode PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(bench_encode BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_executable(bench_replay bench/bench_replay.c)
    target_link_libraries(bench_replay PRIVATE ${MBELIB_EXEC_LINK_TGT})
    _mbelib_apply_hardening(bench_replay)
    target_include_directories(
        bench_replay
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(bench_replay PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(bench_replay BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_executable(bench_tetra bench/bench_tetra.c)
    target_link_libraries(bench_tetra PRIVATE ${MBELIB_EXEC_LINK_TGT})
    _mbelib_apply_hardening(bench_tetra)
    target_include_directories(
        bench_tetra
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(bench_tetra PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(bench_tetra BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_executable(bench_unvoiced bench/bench_unvoiced.c)
    target_link_libraries(bench_unvoiced PRIVATE ${MBELIB_EXEC_LINK_TGT})
    _mbelib_apply_hardening(bench_unvoiced)
    target_include_directories(
        bench_unvoiced
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(bench_unvoiced PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(bench_unvoiced BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    add_executable(bench_convert bench/bench_convert.c)
    target_link_libraries(bench_convert PRIVATE ${MBELIB_EXEC_LINK_TGT})
    _mbelib_apply_hardening(bench_convert)
    target_include_directories(
        bench_convert
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(_SAN_CFLAGS)
        target_compile_options(bench_convert PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(bench_convert BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
endif()

# Quality analysis and private frame fixtures (opt-in)
if(MBELIB_BUILD_TOOLS)
    add_executable(
        mbe_quality_eval
        tools/quality/mbe_quality_eval.c
        tools/quality/mbe_quality_fs.c
    )
    # Enable the POSIX realpath/lstat interfaces used for output alias checks.
    # Win32 uses GetFinalPathNameByHandle/GetFileInformationByHandle instead.
    if(NOT WIN32)
        target_compile_definitions(mbe_quality_eval PRIVATE _XOPEN_SOURCE=700)
    endif()
    # Always measure the shared library selected at runtime. The usual Windows
    # executable target is static, which would bake in one implementation and
    # prevent measuring a different DLL. Linux A/B runs verify the selected
    # library via loader tracing; Windows supports local single measurements.
    target_link_libraries(mbe_quality_eval PRIVATE mbe-shared)
    _mbelib_apply_hardening(mbe_quality_eval)
    if(MSVC)
        target_compile_definitions(
            mbe_quality_eval
            PRIVATE _CRT_SECURE_NO_WARNINGS
        )
    endif()
    target_include_directories(
        mbe_quality_eval
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    add_executable(
        mbe_quality_reframe
        tools/quality/mbe_quality_reframe.c
        tools/quality/mbe_quality_fs.c
        tools/quality/mbe_quality_frames.c
        tools/quality/mbe_quality_dvsi.c
    )
    if(NOT WIN32)
        target_compile_definitions(
            mbe_quality_reframe
            PRIVATE _XOPEN_SOURCE=700
        )
    endif()
    target_link_libraries(mbe_quality_reframe PRIVATE mbe-static)
    _mbelib_apply_hardening(mbe_quality_reframe)
    target_include_directories(
        mbe_quality_reframe
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
    )
    if(MSVC)
        target_compile_definitions(
            mbe_quality_reframe
            PRIVATE _CRT_SECURE_NO_WARNINGS
        )
    endif()
    add_executable(
        mbe_quality_encode
        tools/quality/mbe_quality_encode.c
        tools/quality/mbe_quality_fs.c
    )
    if(NOT WIN32)
        target_compile_definitions(mbe_quality_encode PRIVATE _XOPEN_SOURCE=700)
    endif()
    # Static: each binary measures the encoder it was built with.
    target_link_libraries(mbe_quality_encode PRIVATE mbe-static)
    _mbelib_apply_hardening(mbe_quality_encode)
    target_include_directories(
        mbe_quality_encode
        PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/src/internal"
    )
    if(MSVC)
        target_compile_definitions(
            mbe_quality_encode
            PRIVATE _CRT_SECURE_NO_WARNINGS
        )
    endif()
    # TETRA speech frame file (.CHD layout) to PCM through the public API.
    add_executable(
        mbe_tetra_decode
        tools/quality/mbe_tetra_decode.c
        tools/quality/mbe_quality_fs.c
    )
    if(NOT WIN32)
        target_compile_definitions(mbe_tetra_decode PRIVATE _XOPEN_SOURCE=700)
    endif()
    target_link_libraries(mbe_tetra_decode PRIVATE mbe-static)
    _mbelib_apply_hardening(mbe_tetra_decode)
    target_include_directories(
        mbe_tetra_decode
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
    )
    if(MSVC)
        target_compile_definitions(
            mbe_tetra_decode
            PRIVATE _CRT_SECURE_NO_WARNINGS
        )
    endif()
    get_directory_property(_SAN_CFLAGS MBELIB_SANITIZER_CFLAGS)
    get_directory_property(_SAN_LDFLAGS MBELIB_SANITIZER_LDFLAGS)
    if(_SAN_CFLAGS)
        target_compile_options(mbe_quality_eval PRIVATE ${_SAN_CFLAGS})
        target_compile_options(mbe_quality_reframe PRIVATE ${_SAN_CFLAGS})
        target_compile_options(mbe_quality_encode PRIVATE ${_SAN_CFLAGS})
        target_compile_options(mbe_tetra_decode PRIVATE ${_SAN_CFLAGS})
    endif()
    if(_SAN_LDFLAGS)
        target_link_options(mbe_quality_eval BEFORE PRIVATE ${_SAN_LDFLAGS})
        target_link_options(mbe_quality_reframe BEFORE PRIVATE ${_SAN_LDFLAGS})
        target_link_options(mbe_quality_encode BEFORE PRIVATE ${_SAN_LDFLAGS})
        target_link_options(mbe_tetra_decode BEFORE PRIVATE ${_SAN_LDFLAGS})
    endif()
    # Let the CodeQL extraction target compile the quality tools too.
    target_sources(
        mbe-codeql-first-party
        PRIVATE
            tools/quality/mbe_quality_eval.c
            tools/quality/mbe_quality_reframe.c
            tools/quality/mbe_quality_encode.c
            tools/quality/mbe_tetra_decode.c
            tools/quality/mbe_quality_frames.c
            tools/quality/mbe_quality_dvsi.c
            tools/quality/mbe_quality_fs.c
    )
    target_include_directories(
        mbe-codeql-first-party
        PRIVATE "${PROJECT_SOURCE_DIR}/tools/quality"
    )
    if(NOT WIN32)
        target_compile_definitions(
            mbe-codeql-first-party
            PRIVATE _XOPEN_SOURCE=700
        )
    endif()
endif()

# Uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE
    @ONLY
)

add_custom_target(
    uninstall
    COMMAND
        ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)
