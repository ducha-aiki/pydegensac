CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
PROJECT(MatUtls C)

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-aliasing")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-write-strings")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -ansi")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

include_directories(.)
set(matutls_srcs
ldumat.c
ldvmat.c
mattr.c
minv.c
mmul.c
qrbdv.c
rmmult.c
svduv.c
trnm.c
matconsts.c
ccmath_globals.c)
add_library(matutls ${matutls_srcs})
