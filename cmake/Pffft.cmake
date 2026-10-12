# SPDX-License-Identifier: GPL-2.0-or-later
# Compile the pinned fork into mbelib; do not configure its C++ project.
include(FetchContent)

set(_pffft_pins "${CMAKE_CURRENT_LIST_DIR}/../tools/ci-dependency-pins.env")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_pffft_pins}")
foreach(_pin PFFFT_SHA PFFFT_ARCHIVE_SHA256)
    file(STRINGS "${_pffft_pins}" _lines REGEX "^${_pin}=")
    list(LENGTH _lines _count)
    if(NOT _count EQUAL 1)
        message(FATAL_ERROR "Expected one ${_pin} in ${_pffft_pins}")
    endif()
    string(REPLACE "${_pin}=" "" ${_pin} "${_lines}")
    string(LENGTH "${${_pin}}" _length)
    if(_pin STREQUAL "PFFFT_SHA")
        set(_expected 40)
    else()
        set(_expected 64)
    endif()
    if(NOT _length EQUAL _expected OR NOT "${${_pin}}" MATCHES "^[0-9a-f]+$")
        message(FATAL_ERROR "Invalid ${_pin} in ${_pffft_pins}")
    endif()
endforeach()

if(FETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT)
    if(
        NOT IS_ABSOLUTE "${FETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT}"
        OR NOT IS_DIRECTORY "${FETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT}"
    )
        message(
            FATAL_ERROR
            "FETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT must name an existing absolute PFFFT source directory"
        )
    endif()
endif()
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()
FetchContent_Declare(
    mbelib_pffft
    URL "https://codeload.github.com/arancormonk/pffft/tar.gz/${PFFFT_SHA}"
    URL_HASH "SHA256=${PFFFT_ARCHIVE_SHA256}"
    TLS_VERIFY TRUE
    # Supported since CMake 3.18: a directory without CMakeLists.txt populates
    # sources without add_subdirectory(), preserving our C-only build.
    SOURCE_SUBDIR include
)
FetchContent_MakeAvailable(mbelib_pffft)

set(MBELIB_PFFFT_SOURCE_DIR
    "${mbelib_pffft_SOURCE_DIR}"
    CACHE INTERNAL
    "PFFFT source directory for analyzers"
    FORCE
)
set(PFFFT_SOURCES
    "${MBELIB_PFFFT_SOURCE_DIR}/src/pffft.c"
    "${MBELIB_PFFFT_SOURCE_DIR}/src/pffft_common.c"
)
set(PFFFT_INCLUDE_DIR "${MBELIB_PFFFT_SOURCE_DIR}/include/pffft")
set(_pffft_notice_sources
    src/pffft.c
    src/pffft_priv_impl.h
    include/pffft/pffft.h
    include/pffft/pffft_double.h
    src/simd/pf_float.h
    src/simd/pf_scalar_float.h
    src/simd/pf_sse1_float.h
    src/simd/pf_fma_x86.h
    src/simd/pf_neon_float.h
    src/simd/pf_wasm_float.h
    src/simd/pf_altivec_float.h
)
foreach(_file IN ITEMS src/pffft_common.c LICENSE.txt ${_pffft_notice_sources})
    if(NOT EXISTS "${MBELIB_PFFFT_SOURCE_DIR}/${_file}")
        message(
            FATAL_ERROR
            "Incompatible PFFFT source directory ${MBELIB_PFFFT_SOURCE_DIR}: missing ${_file}. Use the pinned fork; for offline builds set FETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT."
        )
    endif()
endforeach()
file(READ "${MBELIB_PFFFT_SOURCE_DIR}/src/simd/pf_fma_x86.h" _pffft_fma)
if(NOT _pffft_fma MATCHES "PFFFT_DISABLE_FMA")
    message(
        FATAL_ERROR
        "Incompatible PFFFT revision: PFFFT_DISABLE_FMA support is required. Use the pinned fork."
    )
endif()

# Include the license headers of every compiled backend, not just LICENSE.txt.
set(PFFFT_NOTICES "${PROJECT_BINARY_DIR}/PFFFT-NOTICES.txt")
file(
    WRITE "${PFFFT_NOTICES}"
    "PFFFT: https://github.com/arancormonk/pffft\nPinned revision: ${PFFFT_SHA}\n\n"
)
foreach(_file IN LISTS _pffft_notice_sources)
    file(READ "${MBELIB_PFFFT_SOURCE_DIR}/${_file}" _text)
    string(FIND "${_text}" "*/" _end)
    math(EXPR _end "${_end} + 2")
    string(SUBSTRING "${_text}" 0 ${_end} _notice)
    file(APPEND "${PFFFT_NOTICES}" "${_file}\n${_notice}\n\n")
endforeach()

if(MBELIB_ENABLE_WARNINGS AND CMAKE_C_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
    # Upstream's private FFT helpers have external linkage and optional stack
    # scratch uses VLAs. These exceptions apply only to dependency sources.
    set_source_files_properties(
        ${PFFFT_SOURCES}
        PROPERTIES COMPILE_OPTIONS "-Wno-missing-prototypes;-Wno-vla"
    )
endif()
