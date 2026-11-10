#ifndef __DEGENSAC_SIMD_H__
#define __DEGENSAC_SIMD_H__

/* Runtime-dispatched AVX kernels for the per-correspondence loops.

   Wheels are built for the baseline ISA, so wider kernels are selected at
   run time from CPUID. The AVX kernels use only separate multiplies, adds,
   subtractions, divisions and max in the scalar evaluation order (no FMA,
   no reassociation), so both paths produce bit-identical results. This
   assumes that the scalar path is not contracted into FMAs either. That
   holds for the baseline-ISA wheels, but not necessarily for a source build
   with -march=native or another flag that enables FMA.

   PYDEGENSAC_SIMD=0 (or "scalar") in the environment forces the scalar
   path; it is read once, at the first dispatch. */

#include "rtools.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 1 when the AVX kernels are in use, 0 for scalar. */
int degensac_use_avx(void);
/* Name of the active backend: "avx" or "scalar". */
const char *degensac_simd_backend(void);

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#define DEGENSAC_HAVE_AVX_KERNELS 1
void FDs_soa_avx(const double *s, const double *F, double *p, int len);
void FDsSym_soa_avx(const double *s, const double *F, double *p, int len);
Score inlidxs_avx(const double *err, int len, double th, int *inl);
#endif

#ifdef __cplusplus
}
#endif

#endif /* __DEGENSAC_SIMD_H__ */
