/* Numerical regression tests for the fixed-size symmetric 9x9 eigensolver. */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "lapwrap.h"

#define N 9

static uint64_t state = UINT64_C(0x69bf52738a4d1c0e);

static double random_double(void)
{
    state = state * UINT64_C(6364136223846793005) + 1;
    return ((state >> 11) / 9007199254740992.0) * 2.0 - 1.0;
}

static void set_upper(double *a, int i, int j, double value)
{
    if (i > j) { int t = i; i = j; j = t; }
    a[i + N*j] = value;
}

static double upper(const double *a, int i, int j)
{
    if (i > j) { int t = i; i = j; j = t; }
    return a[i + N*j];
}

static void fill_lower_poison(double *a)
{
    int i, j;
    for (j = 0; j < N; ++j)
        for (i = j + 1; i < N; ++i)
            a[i + N*j] = NAN;
}

static void dense(double *a, double scale)
{
    int i, j;
    for (j = 0; j < N; ++j)
        for (i = 0; i <= j; ++i)
            a[i + N*j] = scale * random_double();
    fill_lower_poison(a);
}

static void psd(double *a, int rank, double scale)
{
    double x[N][N], value;
    int i, j, k;
    for (i = 0; i < N; ++i)
        for (k = 0; k < N; ++k)
            x[i][k] = k < rank ? random_double() : 0.0;
    for (j = 0; j < N; ++j)
        for (i = 0; i <= j; ++i) {
            value = 0.0;
            for (k = 0; k < N; ++k)
                value += x[i][k] * x[j][k];
            a[i + N*j] = scale * value;
        }
    fill_lower_poison(a);
}

static void check_success(const char *name, double *input)
{
    double original[N*N], ev[N], norm = 0.0;
    int i, j, k;
    memcpy(original, input, sizeof(original));
    assert(eig9_symmetric(input, ev) == 0);
    for (i = 0; i < N; ++i) {
        assert(isfinite(ev[i]));
        if (i) assert(ev[i - 1] <= ev[i]);
        for (j = 0; j < N; ++j) {
            double dot = 0.0;
            for (k = 0; k < N; ++k) dot += input[k + N*i] * input[k + N*j];
            assert(fabs(dot - (i == j)) < 1e-11);
        }
    }
    for (j = 0; j < N; ++j)
        for (i = 0; i < N; ++i)
            if (fabs(upper(original, i, j)) > norm) norm = fabs(upper(original, i, j));
    if (norm != 0.0) {
        for (j = 0; j < N; ++j)
            for (i = 0; i < N; ++i) {
                double av = 0.0, lv = (ev[j] / norm) * input[i + N*j];
                for (k = 0; k < N; ++k)
                    av += (upper(original, i, k) / norm) * input[k + N*j];
                if (fabs(av - lv) > 1e-11) {
                    fprintf(stderr, "%s: residual too large at (%d,%d): %.17g\n",
                            name, i, j, fabs(av - lv));
                    assert(0);
                }
            }
    }
}

static void check_failure(const char *name, double *a)
{
    double before[N*N], ev[N], ev_before[N];
    int i;
    for (i = 0; i < N; ++i) ev[i] = 1000.0 + i;
    memcpy(before, a, sizeof(before));
    memcpy(ev_before, ev, sizeof(ev_before));
    if (eig9_symmetric(a, ev) == 0) {
        fprintf(stderr, "%s: unexpectedly succeeded\n", name);
        assert(0);
    }
    assert(memcmp(a, before, sizeof(before)) == 0);
    assert(memcmp(ev, ev_before, sizeof(ev)) == 0);
}

int main(void)
{
    double a[N*N];
    int i, j, trial;

    for (trial = 0; trial < 300; ++trial) {
        dense(a, 1.0); check_success("dense", a);
        psd(a, trial % (N + 1), 1.0); check_success("PSD/rank-deficient", a);
    }
    dense(a, 1e200); check_success("large finite scale", a);
    dense(a, 1e-200); check_success("small normal scale", a);
    for (j = 0; j < N; ++j)
        for (i = 0; i <= j; ++i)
            a[i + N*j] = i == j ? (i < 3 ? 4.0 : 1.0) : 0.0;
    fill_lower_poison(a);
    check_success("repeated eigenvalues", a);
    memset(a, 0, sizeof(a));
    fill_lower_poison(a);
    check_success("zero", a);

    dense(a, 1.0); set_upper(a, 0, 4, NAN); check_failure("NaN", a);
    dense(a, 1.0); set_upper(a, 2, 7, INFINITY); check_failure("infinity", a);
    memset(a, 0, sizeof(a));
    set_upper(a, 0, 0, nextafter(0.0, 1.0));
    fill_lower_poison(a);
    check_failure("subnormal", a);
    memset(a, 0, sizeof(a));
    set_upper(a, 0, 0, DBL_MAX); set_upper(a, 0, 1, DBL_MAX);
    set_upper(a, 1, 1, DBL_MAX); fill_lower_poison(a);
    check_failure("rescaling overflow", a);
    puts("eig9 symmetric regression passed");
    return 0;
}
