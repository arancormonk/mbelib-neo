# SPDX-License-Identifier: GPL-2.0-or-later
# Embed the pinned pffft-neo implementation in both library artifacts.
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
set(_pffft_url
    "https://codeload.github.com/arancormonk/pffft-neo/tar.gz/${PFFFT_SHA}"
)
set(_pffft_bundle
    "${PROJECT_SOURCE_DIR}/third_party/pffft-neo-${PFFFT_SHA}.tar.gz"
)
set(_pffft_disconnected "${FETCHCONTENT_FULLY_DISCONNECTED}")
if(NOT FETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT AND EXISTS "${_pffft_bundle}")
    file(SHA256 "${_pffft_bundle}" _pffft_bundle_hash)
    if(NOT _pffft_bundle_hash STREQUAL PFFFT_ARCHIVE_SHA256)
        message(
            FATAL_ERROR
            "Bundled PFFFT SHA256 does not match the dependency pin"
        )
    endif()
    set(_pffft_url "${_pffft_bundle}")
    # A fresh offline bundle still needs local extraction. No network URL is
    # involved, even when the caller requests fully disconnected operation.
    set(FETCHCONTENT_FULLY_DISCONNECTED OFF)
endif()
FetchContent_Declare(
    mbelib_pffft
    URL "${_pffft_url}"
    URL_HASH "SHA256=${PFFFT_ARCHIVE_SHA256}"
    TLS_VERIFY TRUE
    # Supported since CMake 3.18: a directory without CMakeLists.txt populates
    # sources without add_subdirectory(), preserving our C-only build.
    SOURCE_SUBDIR include
)
FetchContent_MakeAvailable(mbelib_pffft)
set(FETCHCONTENT_FULLY_DISCONNECTED "${_pffft_disconnected}")

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
foreach(
    _file
    IN
    ITEMS
        ${PFFFT_SOURCES}
        "${PFFFT_INCLUDE_DIR}/pffft.h"
        "${MBELIB_PFFFT_SOURCE_DIR}/README.md"
        "${MBELIB_PFFFT_SOURCE_DIR}/LICENSE"
        "${MBELIB_PFFFT_SOURCE_DIR}/third_party/pffft/LICENSE.txt"
        "${MBELIB_PFFFT_SOURCE_DIR}/third_party/pffft/LicenseRef-PFFFT.txt"
)
    if(NOT EXISTS "${_file}")
        message(
            FATAL_ERROR
            "Incompatible PFFFT source directory ${MBELIB_PFFFT_SOURCE_DIR}: missing ${_file}. Use the pinned pffft-neo checkout."
        )
    endif()
endforeach()
