/* Regression tests for the guarded PSD 9x9 minimum-eigenvector solver. */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 9

extern int eig9_psd_smallest(const double *a, double *out);

static void householder(double q[N][N])
{
    double norm = 0.0;
    int i, j;
    for (i = 0; i < N; ++i) norm += (i + 1) * (i + 1);
    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j)
            q[i][j] = (i == j ? 1.0 : 0.0) - 2.0 * (i + 1) * (j + 1) / norm;
}

static void upper_gapped(double a[N*N], double scale)
{
    static const double eigenvalues[N] = {0.05, 0.3, 0.8, 1.4, 2.0, 3.0, 4.0, 6.0, 9.0};
    double q[N][N];
    int i, j, k;
    householder(q);
    for (j = 0; j < N; ++j)
        for (i = 0; i <= j; ++i) {
            double value = 0.0;
            for (k = 0; k < N; ++k) value += q[i][k] * eigenvalues[k] * q[j][k];
            a[i + N*j] = scale * value;
        }
    for (j = 0; j < N; ++j)
        for (i = j + 1; i < N; ++i) a[i + N*j] = NAN;
}

static double upper_at(const double *a, int i, int j)
{
    return i <= j ? a[i + N*j] : a[j + N*i];
}

static void check_gapped(double scale)
{
    double a[N*N], before[N*N], out[N], q[N][N], norm = 0.0, dot = 0.0,
           rayleigh = 0.0, max_residual = 0.0;
    int i, j;
    upper_gapped(a, scale); memcpy(before, a, sizeof(a)); householder(q);
    assert(eig9_psd_smallest(a, out) == 0);
    assert(memcmp(a, before, sizeof(a)) == 0);
    for (i = 0; i < N; ++i) {
        double av = 0.0;
        norm += out[i] * out[i]; dot += out[i] * q[i][0];
        for (j = 0; j < N; ++j) av += upper_at(a, i, j) / scale * out[j];
        rayleigh += out[i] * av;
    }
    assert(fabs(norm - 1.0) < 2e-10);
    assert(fabs(fabs(dot) - 1.0) < 2e-9);
    assert(fabs(rayleigh - 0.05) < 2e-11);
    for (i = 0; i < N; ++i) {
        double residual = 0.0;
        for (j = 0; j < N; ++j) residual += upper_at(a, i, j) / scale * out[j];
        residual -= rayleigh * out[i];
        max_residual = fmax(max_residual, fabs(residual));
    }
    assert(max_residual < 5e-11);
}

static void check_fallback(const char *name, double a[N*N])
{
    double before[N*N], out[N], out_before[N];
    int i;
    for (i = 0; i < N; ++i) out[i] = 100.0 + i;
    memcpy(before, a, sizeof(before)); memcpy(out_before, out, sizeof(out));
    if (eig9_psd_smallest(a, out) == 0) {
        fprintf(stderr, "%s unexpectedly succeeded\n", name); assert(0);
    }
    assert(memcmp(a, before, sizeof(before)) == 0);
    assert(memcmp(out, out_before, sizeof(out)) == 0);
}

int main(void)
{
    double a[N*N], out[N];
    int i, j;
    check_gapped(1.0); check_gapped(1e-200); check_gapped(1e200);
    memset(a, 0, sizeof(a));
    for (j = 0; j < N; ++j) for (i = j + 1; i < N; ++i) a[i + N*j] = NAN;
    assert(eig9_psd_smallest(a, out) == 0);
    assert(fabs(out[0] - 1.0) == 0.0);
    for (i = 1; i < N; ++i) assert(out[i] == 0.0);

    upper_gapped(a, 1.0); a[0 + N*4] = NAN; check_fallback("NaN", a);
    upper_gapped(a, 1.0); a[2 + N*7] = INFINITY; check_fallback("infinity", a);
    memset(a, 0, sizeof(a)); a[0] = nextafter(0.0, 1.0); check_fallback("subnormal", a);
    memset(a, 0, sizeof(a));
    for (i = 0; i < N; ++i) a[i + N*i] = i == 0 ? -0.1 : (double)i;
    check_fallback("indefinite", a);
    memset(a, 0, sizeof(a));
    for (i = 0; i < N; ++i) a[i + N*i] = i < 2 ? 0.0 : (double)i;
    check_fallback("repeated minimum", a);
    puts("PSD eig9 minimum-eigenvector regression passed");
    return 0;
}
