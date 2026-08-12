cmake_minimum_required(VERSION 3.9)
project(pydegensac)

add_definitions(-fPIC)

FIND_PACKAGE(LAPACK REQUIRED)
if (LAPACK_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LAPACK_CXX_FLAGS}")
endif (LAPACK_FOUND)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Link-time optimisation. The estimators' hot calls cross translation units --
# the error functions go through a function pointer, and inlidxs/cov_mat live
# in their own files -- so there is real cross-TU inlining available. Measured
# on M1/Clang 2026-08-12: F +3.5% (342 -> 354 estimator calls per 10 s), H
# within noise, and bit-identical outputs.
#
# Guarded three ways, because this is the wheel build for five platforms:
# check_ipo_supported() skips toolchains that cannot do it, PYDEGENSAC_NO_LTO=1
# turns it off for bisecting, and -fno-strict-aliasing rides along -- this C
# core predates the strict-aliasing rule (the C++ half already silences the
# warning), and LTO is precisely what would start acting on it.
#
# Also measured and NOT adopted: -funroll-loops/-ftree-vectorize on the C core
# (noise), -fno-math-errno (Clang already defaults to it here), and
# -mcpu=apple-m1 (no gain, so portable wheels cost nothing on this ISA).
include(CheckIPOSupported)
check_ipo_supported(RESULT PYDEGENSAC_IPO_OK OUTPUT PYDEGENSAC_IPO_MSG)
if(PYDEGENSAC_IPO_OK AND NOT DEFINED ENV{PYDEGENSAC_NO_LTO})
  message(STATUS "pydegensac: enabling LTO")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  if(NOT MSVC)
    add_compile_options(-fno-strict-aliasing)
  endif()
else()
  message(STATUS "pydegensac: LTO off (${PYDEGENSAC_IPO_MSG})")
endif()

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-aliasing")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-write-strings")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory(src/pydegensac/matutls)
add_subdirectory(lib/pybind11)


####### DEGENSAC
include_directories(
src/pydegensac
src/pydegensac/degensac)
set(degensac_srcs
        src/pydegensac/degensac/bsd_random.c
        src/pydegensac/degensac/DegUtils.c
        src/pydegensac/degensac/exp_ranF.c
        src/pydegensac/degensac/exp_ranH.c
        src/pydegensac/degensac/Ftools.c
        src/pydegensac/degensac/hash.c
        src/pydegensac/degensac/Htools.c
        src/pydegensac/degensac/ranH.c
        src/pydegensac/degensac/rtools.c
        src/pydegensac/degensac/utools.c
        src/pydegensac/degensac/lapwrap.c
)

add_library(pydegensac_support ${degensac_srcs})
# pydegensac_support calls into matutls (mattr, mattr's neighbours in ccmath.h),
# so it has to come first on the link line: a static archive only yields the
# members needed by what the linker has already seen. Stating the dependency
# here rather than relying on the order in target_link_libraries below lets
# CMake order and repeat the archives itself.
#
# Without LTO the flat order happened to work. With it, GCC's linker plugin
# resolves archive members differently and the reference to `mattr` from
# DegUtils.c/exp_ranH.c arrived after matutls had already been passed, so every
# LTO build on Linux/GCC linked an extension that failed to import with
# "undefined symbol: mattr" — both manylinux images and local gcc 13.3.
target_link_libraries(pydegensac_support PUBLIC matutls)
pybind11_add_module(pydegensac "src/pydegensac/bindings.cpp")
target_link_libraries(pydegensac PRIVATE ${LAPACK_LIBRARIES} pydegensac_support)
