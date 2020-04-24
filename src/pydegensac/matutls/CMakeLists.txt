CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(MatUtls)

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-aliasing")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-write-strings")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -ansi")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

include_directories(.)
set(matutls_srcs
atou1.c   cmprt.c    hevmax.c
mcopy.c   qrbdv.c   solvru.c
trncm.c   atovm.c   csolv.c
hmgen.c   minv.c    qrecvc.c
solvtd.c  trnm.c    chouse.c  cvmul.c    house.c
mmul.c    qreval.c  sv2u1v.c  unitary.c
chousv.c  eigen.c    housev.c  ortho.c   qrevec.c  sv2uv.c   utrncm.c
cmattr.c  eigval.c   ldumat.c  otrma.c
qrlsq.c   sv2val.c  utrnhm.c
cmcpy.c   evmax.c    ldvmat.c  otrsm.c   rmmult.c  vmul.c
cminv.c   hconj.c    lsqsv.c   psinv.c   ruinv.c   svdu1v.c
cmmul.c   heigval.c  matprt.c  qrbdi.c   smgen.c   svduv.c
cmmult.c  heigvec.c  mattr.c   qrbdu1.c  solvps.c  svdval.c
)
add_library(matutls ${matutls_srcs})
