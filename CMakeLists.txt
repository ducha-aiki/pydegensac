cmake_minimum_required(VERSION 3.15...3.29)
project(pydegensac LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

add_definitions(-fPIC)

find_package(LAPACK REQUIRED)

add_subdirectory(src/pydegensac/matutls)
add_subdirectory(lib/pybind11)

include_directories(
  src/pydegensac
  src/pydegensac/degensac
)

set(degensac_srcs
  src/pydegensac/degensac/DegUtils.c
  src/pydegensac/degensac/exp_ranF.c
  src/pydegensac/degensac/exp_ranH.c
  src/pydegensac/degensac/Ftools.c
  src/pydegensac/degensac/hash.c
  src/pydegensac/degensac/Htools.c
  src/pydegensac/degensac/ranF.c
  src/pydegensac/degensac/ranH2el.c
  src/pydegensac/degensac/ranH.c
  src/pydegensac/degensac/rtools.c
  src/pydegensac/degensac/utools.c
  src/pydegensac/degensac/lapwrap.c
)

add_library(pydegensac_support ${degensac_srcs})

# Core speed flags for BOTH C and C++
target_compile_options(pydegensac_support PRIVATE
  $<$<COMPILE_LANGUAGE:C>:
    -O3 -fno-math-errno -ffp-contract=fast
  >
  $<$<COMPILE_LANGUAGE:CXX>:
    -O3 -fno-math-errno -ffp-contract=fast
  >
)

# Optional local-only (don't use for distributing wheels):
# target_compile_options(pydegensac_support PRIVATE -march=native)

pybind11_add_module(pydegensac "src/pydegensac/bindings.cpp")
target_compile_options(pydegensac PRIVATE -O3 -fno-math-errno -ffp-contract=fast)

target_link_libraries(pydegensac_support PRIVATE ${LAPACK_LIBRARIES} matutls)
target_link_libraries(pydegensac PRIVATE ${LAPACK_LIBRARIES} matutls pydegensac_support)
