/* AVX kernels and CPU dispatch; see simd.h.
   Each kernel mirrors its scalar twin operation by operation: the comments
   give the scalar expression with its implicit left-to-right association. */
#include <stdlib.h>
#include <string.h>
#include "simd.h"
#include "utools.h"

static int simd_level = -1;

#ifdef DEGENSAC_HAVE_AVX_KERNELS

#include <immintrin.h>
/* MSVC proper accepts AVX intrinsics anywhere; GCC and Clang (clang-cl too)
   need the target attribute. Never add "fma" to it: contracted multiply-adds
   would break bit-equality with the scalar path. */
#if defined(__GNUC__) || defined(__clang__)
#define DEGENSAC_AVX_FN __attribute__((target("avx")))
#else
#define DEGENSAC_AVX_FN
#endif
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#define DEGENSAC_MSVC_CPUID 1
#else
#include <cpuid.h>
#endif

static int cpu_has_avx(void)
{
    unsigned long long xcr0;
#ifdef DEGENSAC_MSVC_CPUID
    int r[4];
    __cpuid(r, 1);
    if (!((r[2] >> 27) & 1) || !((r[2] >> 28) & 1)) /* OSXSAVE, AVX */
        return 0;
    xcr0 = _xgetbv(0);
#else
    unsigned a, b, c, d, lo, hi;
    if (!__get_cpuid(1, &a, &b, &c, &d))
        return 0;
    if (!(c & bit_OSXSAVE) || !(c & bit_AVX))
        return 0;
    __asm__ volatile ("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));
    xcr0 = ((unsigned long long)hi << 32) | lo;
#endif
    return (xcr0 & 6) == 6; /* OS saves XMM and YMM state */
}

#else
static int cpu_has_avx(void) { return 0; }
#endif /* DEGENSAC_HAVE_AVX_KERNELS */

static int detect_level(void)
{
    const char *env = getenv("PYDEGENSAC_SIMD");
    if (env != NULL && (strcmp(env, "0") == 0 || strcmp(env, "scalar") == 0))
        return 0;
    return cpu_has_avx();
}

/* Racing first calls compute and store the same value. */
int degensac_use_avx(void)
{
    int level = simd_level;
    if (level < 0) {
        level = detect_level();
        simd_level = level;
    }
    return level;
}

const char *degensac_simd_backend(void)
{
    return degensac_use_avx() ? "avx" : "scalar";
}

#ifdef DEGENSAC_HAVE_AVX_KERNELS

#define LD(p) _mm256_loadu_pd(p)
#define ADD _mm256_add_pd
#define SUB _mm256_sub_pd
#define MUL _mm256_mul_pd
#define DIV _mm256_div_pd

/* Sampson error, FDs_soa. */
DEGENSAC_AVX_FN
void FDs_soa_avx(const double *s, const double *F, double *p, int len)
{
    const double *x1 = s, *y1 = s + len, *x2 = s + 2*len, *y2 = s + 3*len;
    const __m256d f1 = _mm256_set1_pd(_f1), f2 = _mm256_set1_pd(_f2),
        f3 = _mm256_set1_pd(_f3), f4 = _mm256_set1_pd(_f4),
        f5 = _mm256_set1_pd(_f5), f6 = _mm256_set1_pd(_f6),
        f7 = _mm256_set1_pd(_f7), f8 = _mm256_set1_pd(_f8),
        f9 = _mm256_set1_pd(_f9);
    int i, n4 = len & ~3;

    for (i = 0; i < n4; i += 4) {
        __m256d a1 = LD(x1+i), a2 = LD(y1+i), b1 = LD(x2+i), b2 = LD(y2+i);
        /* rxc = _f1*b1 + _f4*b2 + _f7, and likewise ryc, rwc */
        __m256d rxc = ADD(ADD(MUL(f1, b1), MUL(f4, b2)), f7);
        __m256d ryc = ADD(ADD(MUL(f2, b1), MUL(f5, b2)), f8);
        __m256d rwc = ADD(ADD(MUL(f3, b1), MUL(f6, b2)), f9);
        /* r = a1*rxc + a2*ryc + rwc */
        __m256d r = ADD(ADD(MUL(a1, rxc), MUL(a2, ryc)), rwc);
        /* rx = _f1*a1 + _f2*a2 + _f3; ry = _f4*a1 + _f5*a2 + _f6 */
        __m256d rx = ADD(ADD(MUL(f1, a1), MUL(f2, a2)), f3);
        __m256d ry = ADD(ADD(MUL(f4, a1), MUL(f5, a2)), f6);
        /* r*r / (rxc*rxc + ryc*ryc + rx*rx + ry*ry) */
        __m256d den = ADD(ADD(ADD(MUL(rxc, rxc), MUL(ryc, ryc)), MUL(rx, rx)),
                          MUL(ry, ry));
        _mm256_storeu_pd(p + i, DIV(MUL(r, r), den));
    }
    _mm256_zeroupper();
    for (i = n4; i < len; i++) {
        double a1 = x1[i], a2 = y1[i], b1 = x2[i], b2 = y2[i];
        double rxc = _f1 * b1 + _f4 * b2 + _f7;
        double ryc = _f2 * b1 + _f5 * b2 + _f8;
        double rwc = _f3 * b1 + _f6 * b2 + _f9;
        double r = (a1 * rxc + a2 * ryc + rwc);
        double rx = _f1 * a1 + _f2 * a2 + _f3;
        double ry = _f4 * a1 + _f5 * a2 + _f6;
        p[i] = r*r / (rxc*rxc + ryc*ryc + rx*rx + ry*ry);
    }
}

/* Symmetric epipolar distance, FDsSym_soa. */
DEGENSAC_AVX_FN
void FDsSym_soa_avx(const double *s, const double *F, double *p, int len)
{
    const double *x1 = s, *y1 = s + len, *x2 = s + 2*len, *y2 = s + 3*len;
    const __m256d f1 = _mm256_set1_pd(_f1), f2 = _mm256_set1_pd(_f2),
        f3 = _mm256_set1_pd(_f3), f4 = _mm256_set1_pd(_f4),
        f5 = _mm256_set1_pd(_f5), f6 = _mm256_set1_pd(_f6),
        f7 = _mm256_set1_pd(_f7), f8 = _mm256_set1_pd(_f8),
        f9 = _mm256_set1_pd(_f9);
    int i, n4 = len & ~3;

    for (i = 0; i < n4; i += 4) {
        __m256d a1 = LD(x1+i), a2 = LD(y1+i), b1 = LD(x2+i), b2 = LD(y2+i);
        __m256d rxc = ADD(ADD(MUL(f1, b1), MUL(f4, b2)), f7);
        __m256d ryc = ADD(ADD(MUL(f2, b1), MUL(f5, b2)), f8);
        __m256d rwc = ADD(ADD(MUL(f3, b1), MUL(f6, b2)), f9);
        __m256d r = ADD(ADD(MUL(a1, rxc), MUL(a2, ryc)), rwc);
        __m256d rx = ADD(ADD(MUL(f1, a1), MUL(f2, a2)), f3);
        __m256d ry = ADD(ADD(MUL(f4, a1), MUL(f5, a2)), f6);
        /* a = rxc*rxc + ryc*ryc; b = rx*rx + ry*ry; r*r*(a+b) / (a*b) */
        __m256d a = ADD(MUL(rxc, rxc), MUL(ryc, ryc));
        __m256d b = ADD(MUL(rx, rx), MUL(ry, ry));
        _mm256_storeu_pd(p + i, DIV(MUL(MUL(r, r), ADD(a, b)), MUL(a, b)));
    }
    _mm256_zeroupper();
    for (i = n4; i < len; i++) {
        double a1 = x1[i], a2 = y1[i], b1 = x2[i], b2 = y2[i];
        double rxc = _f1 * b1 + _f4 * b2 + _f7;
        double ryc = _f2 * b1 + _f5 * b2 + _f8;
        double rwc = _f3 * b1 + _f6 * b2 + _f9;
        double r = (a1 * rxc + a2 * ryc + rwc);
        double rx = _f1 * a1 + _f2 * a2 + _f3;
        double ry = _f4 * a1 + _f5 * a2 + _f6;
        double a = rxc*rxc + ryc*ryc;
        double b = rx*rx + ry*ry;
        p[i] = r*r* (a+b)/(a*b);
    }
}

/* inlidxs with th != 0. Its four partial gain sums j0..j3 are exactly the
   four lanes here, so the final (j0 + j1) + (j2 + j3) is unchanged.
   max_pd(g, 0) returns its second operand unless g > 0, which is
   truncQuadInv's g > 0.0 ? g : 0.0, NaN included. */
DEGENSAC_AVX_FN
Score inlidxs_avx(const double *err, int len, double th, int *inl)
{
    Score s = {0,0,0,0};
    const double inv = 1/(th*9/4);
    const __m256d vinv = _mm256_set1_pd(inv), one = _mm256_set1_pd(1.0),
        zero = _mm256_setzero_pd(), vth = _mm256_set1_pd(th);
    __m256d j = _mm256_setzero_pd();
    double jl[4];
    int i, n4 = len & ~3;

    for (i = 0; i < n4; i += 4) {
        __m256d e = LD(err + i);
        int m;
        j = ADD(j, _mm256_max_pd(SUB(one, MUL(e, vinv)), zero));
        m = _mm256_movemask_pd(_mm256_cmp_pd(e, vth, _CMP_LE_OQ));
        inl[s.I] = i;     s.I += (m & 1);
        inl[s.I] = i+1;   s.I += (m >> 1) & 1;
        inl[s.I] = i+2;   s.I += (m >> 2) & 1;
        inl[s.I] = i+3;   s.I += (m >> 3) & 1;
    }
    _mm256_storeu_pd(jl, j);
    _mm256_zeroupper();
    for (i = n4; i < len; ++i) {
        jl[0] += truncQuadInv(err[i], inv);
        inl[s.I] = i;
        s.I += (err[i] <= th);
    }
    s.J = (jl[0] + jl[1]) + (jl[2] + jl[3]);
    return s;
}

#endif /* DEGENSAC_HAVE_AVX_KERNELS */
