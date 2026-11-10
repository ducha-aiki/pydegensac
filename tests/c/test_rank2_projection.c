/* Regression tests for the guarded 3x3 rank-two projection. */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

extern int rank2_project3(double *f);

#define N 3

static void multiply(double out[N][N], const double a[N][N],
                     const double b[N][N])
{
    int i, j, k;
    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j) {
            double sum = 0.0;
            for (k = 0; k < N; ++k) sum += a[i][k] * b[k][j];
            out[i][j] = sum;
        }
}

static void transpose(double out[N][N], const double a[N][N])
{
    int i, j;
    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j) out[i][j] = a[j][i];
}

/* Build U diag(s) V' using fixed, non-axis-aligned orthonormal rotations. */
static void spectrum(double out[9], const double s[3], double scale)
{
    const double U[N][N] = {{0.8, -0.6, 0.0},
                            {0.36, 0.48, -0.8},
                            {0.48, 0.64, 0.6}};
    const double V[N][N] = {{0.6, 0.0, 0.8},
                            {0.64, 0.6, -0.48},
                            {-0.48, 0.8, 0.36}};
    double d[N][N] = {{0}}, ud[N][N], vt[N][N], result[N][N];
    int i, j;

    for (i = 0; i < N; ++i) d[i][i] = s[i] * scale;
    multiply(ud, U, d);
    transpose(vt, V);
    multiply(result, ud, vt);
    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j) out[3*i + j] = result[i][j];
}

static double frobenius(const double a[9])
{
    double sum = 0.0;
    int i;
    for (i = 0; i < 9; ++i) sum += a[i] * a[i];
    return sqrt(sum);
}

static double determinant_normalized(const double a[9], double scale)
{
    double b[9];
    int i;
    for (i = 0; i < 9; ++i) b[i] = a[i] / scale;
    return b[0] * (b[4]*b[8] - b[5]*b[7])
         - b[1] * (b[3]*b[8] - b[5]*b[6])
         + b[2] * (b[3]*b[7] - b[4]*b[6]);
}

static void check_projection(const double singulars[3], double scale)
{
    double input[9], projected[9], difference[9];
    int i;

    spectrum(input, singulars, scale);
    memcpy(projected, input, sizeof(input));
    assert(rank2_project3(projected) == 0);
    /* Normalize before squaring so the finite extreme-scale cases stay
       meaningful rather than overflowing the test's Frobenius accumulator. */
    for (i = 0; i < 9; ++i)
        difference[i] = (input[i] - projected[i]) / scale;
    assert(fabs(frobenius(difference) - singulars[2]) < 1e-12);
    assert(fabs(determinant_normalized(projected, scale)) < 1e-12);
}

static void check_fallback(const char *name, double a[9])
{
    double before[9];
    memcpy(before, a, sizeof(before));
    if (rank2_project3(a) == 0) {
        fprintf(stderr, "%s unexpectedly used the fast path\n", name);
        assert(0);
    }
    assert(memcmp(a, before, sizeof(before)) == 0);
}

int main(void)
{
    const double well_conditioned[3] = {0.9, 0.5, 0.05};
    const double near_rank_one[3] = {1.0, 1e-6, 1e-10};
    const double tied_smallest[3] = {1.0, 0.5, 0.5};
    const double below_guard[3] = {1.0, 1e-4, 1e-9};
    const double above_guard[3] = {1.0, 1e-3, 1e-5};
    double a[9];
    int i;

    check_projection(well_conditioned, 1.0);
    check_projection(above_guard, 1.0);
    check_projection(well_conditioned, 1e200);
    check_projection(well_conditioned, 1e-200);

    spectrum(a, below_guard, 1.0); check_fallback("below guard", a);
    spectrum(a, near_rank_one, 1.0); check_fallback("near rank one", a);
    spectrum(a, tied_smallest, 1.0); check_fallback("tied smallest", a);
    memset(a, 0, sizeof(a));
    assert(rank2_project3(a) == 0);
    for (i = 0; i < 9; ++i) assert(a[i] == 0.0);

    spectrum(a, well_conditioned, 1.0);
    a[4] = NAN; check_fallback("NaN", a);
    spectrum(a, well_conditioned, 1.0);
    a[7] = INFINITY; check_fallback("infinity", a);
    puts("rank-two projection regression passed");
    return 0;
}
