cmake_minimum_required(VERSION 2.8.12)
project(pydegensac)

add_definitions(-DA64)
add_definitions(-std=c++11)
add_definitions(-fPIC)



FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif (OPENMP_FOUND)


enable_language (Fortran)
FIND_PACKAGE(LAPACK REQUIRED)
if (LAPACK_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LAPACK_CXX_FLAGS}")
endif (LAPACK_FOUND)

SET(CMAKE_BUILD_TYPE "RELEASE")

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-aliasing")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-write-strings")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -ansi")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory(src/pydegensac/matutls)
add_subdirectory(lib/pybind11)


####### DEGENSAC
include_directories(
src/pydegensac
src/pydegensac/degensac)
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


pybind11_add_module(pydegensac "src/pydegensac/bindings.cpp" ${degensac_srcs})
target_link_libraries(pydegensac PRIVATE ${LAPACK_LIBRARIES} matutls)
