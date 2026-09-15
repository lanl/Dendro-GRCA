# ----------
# BASE CMAKE CONFIGURATION
#
# This file in particular contains all of the options that apply to all aspects
# of the program, including some Dendro options that are used when including the
# library.
#
# Note that all options and compiler flags set here will propagate to other
# directories such as src/ and dendro/

cmake_minimum_required(VERSION 3.15)
project(DendroGR C CXX Fortran)

if(POLICY CMP0169)
  cmake_policy(SET CMP0169 OLD)
endif()

# set C++ standard to 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# OpenMP, MPI, and GSL are required packages for all aspects of the project
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(GSL REQUIRED)

# =============================
# == CPU ARCH OPTIONS ==
# =============================
set(CPU_ARCH "native" CACHE STRING "Target CPU Architecture, see list")

# this list of strings are the valid list of options we have available for compiling
set(CPU_ARCH_STRINGS
    "native"
    "znver1" "znver2" "znver3" "znver4"
    "haswell" "broadwell" "skylake-avx512" "cascadelake" "icelake-server"
    "sapphirerapids"
    "generic_avx2" "x86-64-v3"
)
set_property(CACHE CPU_ARCH PROPERTY STRINGS ${CPU_ARCH_STRINGS})

set (ARCH_FLAGS "")

if(CPU_ARCH STREQUAL "generic_avx2")
  # this is like a custom alias for generic avx2
  list(APPEND ARCH_FLAGS "-mavx2" "-mfma")
elseif(CPU_ARCH IN_LIST CPU_ARCH_STRINGS)
  # we pass the name directly to "march" if it's in the list of strings
  list(APPEND ARCH_FLAGS "-march=${CPU_ARCH}")
else()
  message(FATAL_ERROR "Unknown CPU_ARCH: ${CPU_ARCH}")
endif()

# A named -march is only understood by new enough compilers (e.g. sapphirerapids needs
# GCC >= 11 / Clang >= 12 / ICX). Fail here with the remedy rather than deep in a TU.
include(CheckCXXCompilerFlag)
foreach(_f IN LISTS ARCH_FLAGS)
  string(MAKE_C_IDENTIFIER "HAVE_ARCH${_f}" _have)
  check_cxx_compiler_flag("${_f}" ${_have})
  if(NOT ${_have})
    message(FATAL_ERROR
      "CPU_ARCH=${CPU_ARCH} produces '${_f}', which ${CMAKE_CXX_COMPILER_ID} "
      "${CMAKE_CXX_COMPILER_VERSION} does not accept. Load a newer compiler, or pick a "
      "CPU_ARCH it supports (e.g. icelake-server for an AVX-512 host). Do NOT silently "
      "fall back to a different arch: it changes the ISA under a benchmark.")
  endif()
endforeach()

# by setting compile options, we're forcing things to be set properly
add_compile_options(${ARCH_FLAGS})

message(STATUS "Build Target Architecture: ${CPU_ARCH}")
message(STATUS "Compiler Flags: ${ARCH_FLAGS}")

# set up fsanatize flags to check for memory leaks (Debug only)
# if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
#   set(SAN_FLAGS "-fsanitize=address,undefined" "-fno-sanitize-recover=all")
#   foreach(FLAG ${SAN_FLAGS})
#     add_compile_options($<$<CONFIG:Debug>:${FLAG}>)
#     add_link_options($<$<CONFIG:Debug>:${FLAG}>)
#   endforeach()
#   message(STATUS "Sanitizers configured for Debug builds (Address + Undefined)")
# endif()


# then we need to enable CUDA based on the user options
option(WITH_CUDA "Build dendro with CUDA" OFF)

if(WITH_CUDA)
  enable_language(CUDA)
endif()

# =============================
# Project Configuration and Dependencies
# =============================
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

if(MANUAL_BLAS_LAPACK)
  if("$ENV{BLAS}" STREQUAL "" OR "$ENV{LAPACK}" STREQUAL "")
    message(FATAL_ERROR "MANUAL_BLAS_LAPACK is ON, but BLAS/LAPACK env vars are empty.")
  endif()

  message(STATUS "Configuring BLAS/LAPACK Manually...")

  set(LAPACKE_DIR "$ENV{LAPACK}/LAPACKE")
  set(BLAS_LIBS "$ENV{BLAS}/lib")
  set(LAPACK_LIBS "$ENV{LAPACK}/lib")

  set(MANUAL_LIBS
        "${LAPACK_LIBS}/liblapacke.a"
        "${LAPACK_LIBS}/liblapack.a"
        "${BLAS_LIBS}/libblas.a"
        "-static" "libgfortran.a" "libquadmath.a"
    )

  if(NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${MANUAL_LIBS}")
  endif()

  if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${MANUAL_LIBS}")
  endif()

  # Prevent dendrolib from trying to find_package and overwriting this
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
endif()


set(CHECK_GIT_INCLUDE_PATH ${CMAKE_SOURCE_DIR}/cmake-modules)
include(${CHECK_GIT_INCLUDE_PATH}/CheckGit.cmake)
message(
  STATUS
    "Setting up adding Git Information directly to compile options, header file should be in ${CHECK_GIT_INCLUDE_PATH}"
)
checkgitsetup()

include(FetchContent)
# we require toml11
FetchContent_Declare(
    toml11
    GIT_REPOSITORY "https://github.com/ToruNiina/toml11"
    GIT_TAG        "v4.4.0"
)
FetchContent_MakeAvailable(toml11)

# settable option for the user on where dendrolib is available
set(DENDRO_dendrolib_DIR "" CACHE PATH "Path to the local dendrolib source directory (optional)")
# Default to the dfvk fork's experimental branch: it carries bh_history.h
# (BHHistory tracker) and the compact-FD work this branch requires. Upstream
# paralab/Dendro-5.01 master lacks these and fails to build. Override either
# cache var, or set DENDRO_dendrolib_DIR, to use a different dendrolib.
set(DENDRO_dendrolib_GIT_REPOSITORY "https://github.com/dfvankomen/Dendro-5.01" CACHE STRING "Git repository URL for dendrolib")
set(DENDRO_dendrolib_GIT_TAG "experimental" CACHE STRING "The Git tag, branch, or commit hash for dendrolib")


# now include the dendro repository for compiling the library
if (IS_DIRECTORY ${DENDRO_dendrolib_DIR})
  message(STATUS "Using user-provided dendrolib from: ${DENDRO_dendrolib_DIR}")
  add_subdirectory(${DENDRO_dendrolib_DIR} dendrolib)
  set(DENDRO_FETCHED FALSE)
else()
  message(STATUS "Dendrolib not found locally, fetching from repository...")
  FetchContent_Declare(
    dendrolib
    GIT_REPOSITORY ${DENDRO_dendrolib_GIT_REPOSITORY}
    GIT_TAG ${DENDRO_dendrolib_GIT_TAG}
  )

  FetchContent_MakeAvailable(dendrolib)

  set(DENDRO_dendrolib_DIR ${dendrolib_SOURCE_DIR} CACHE INTERNAL "Path to the fetched dendrolib source")

endif()

add_subdirectory(BSSN_GR)
add_subdirectory(Fluid)
