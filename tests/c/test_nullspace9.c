/* Exact regression gate for the fixed-size 9x9 nullspace path used by H. */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utools.h"

#define N 9
#define SENTINEL 0x5a

/* This harness never reaches cov_mat's BLAS path, but utools.c needs it. */
void dsyrk_(char *uplo, char *trans, ptrdiff_t *n, ptrdiff_t *k,
            double *alpha, double *a, ptrdiff_t *lda, double *beta,
            double *c, ptrdiff_t *ldc)
{
    (void)uplo; (void)trans; (void)n; (void)k; (void)alpha; (void)a;
    (void)lda; (void)beta; (void)c; (void)ldc;
    abort();
}

typedef struct {
    unsigned char before[16];
    double values[N*N];
    unsigned char after[16];
} GuardedBasis;

static uint64_t state = UINT64_C(0x168db26c4f7a9301);

static double random_double(void)
{
    state = state * UINT64_C(6364136223846793005) + 1;
    return ((state >> 11) * 0x1.0p-53) * 2.0 - 1.0;
}

static void check_case(const char *name, const double *source)
{
    double generic_matrix[N*N], fixed_matrix[N*N];
    GuardedBasis generic, fixed;
    int generic_buffer[2*N], fixed_buffer[2*N];
    int generic_size, fixed_size, i;

    memcpy(generic_matrix, source, sizeof(generic_matrix));
    memcpy(fixed_matrix, source, sizeof(fixed_matrix));
    memset(&generic, SENTINEL, sizeof(generic));
    memset(&fixed, SENTINEL, sizeof(fixed));
    memset(generic_buffer, SENTINEL, sizeof(generic_buffer));
    memset(fixed_buffer, SENTINEL, sizeof(fixed_buffer));
    generic_size = nullspace(generic_matrix, generic.values, N, generic_buffer);
    fixed_size = nullspace9(fixed_matrix, fixed.values, fixed_buffer);
    if (generic_size != fixed_size ||
        memcmp(generic.values, fixed.values, sizeof(generic.values)) != 0 ||
        memcmp(generic_buffer, fixed_buffer, sizeof(generic_buffer)) != 0) {
        fprintf(stderr, "%s: basis, rank, or pivot buffer differs (%d vs %d)\n",
                name, generic_size, fixed_size);
        assert(0);
    }
    for (i = 0; i < (int)sizeof(generic.before); ++i) {
        assert(generic.before[i] == SENTINEL && fixed.before[i] == SENTINEL);
        assert(generic.after[i] == SENTINEL && fixed.after[i] == SENTINEL);
    }
}

static void random_matrix(double *matrix, double scale)
{
    int i;
    for (i = 0; i < N*N; ++i)
        matrix[i] = scale * random_double();
}

/* Exact layout produced by exp_ranH.c's sampleH4: four correspondences and
   the ninth zero row that makes the 8x9 DLT system square. */
static void h_matrix(double *matrix, double noise)
{
    int i;
    for (i = 0; i < 4; ++i) {
        double x1 = random_double(), y1 = random_double();
        double x2 = x1 + 0.3 + noise * random_double();
        double y2 = y1 - 0.2 + noise * random_double();
        double *a = matrix + 18*i, *b = a + N;
        a[0]=x2; a[1]=0; a[2]=-x1*x2; a[3]=y2; a[4]=0; a[5]=-x1*y2;
        a[6]=1; a[7]=0; a[8]=-x1;
        b[0]=0; b[1]=x2; b[2]=-y1*x2; b[3]=0; b[4]=y2; b[5]=-y1*y2;
        b[6]=0; b[7]=1; b[8]=-y1;
    }
    memset(matrix + 8*N, 0, N * sizeof(*matrix));
}

int main(void)
{
    double matrix[N*N];
    int i, trial;

    for (trial = 0; trial < 5000; ++trial) {
        random_matrix(matrix, 1.0);
        check_case("full random", matrix);
        h_matrix(matrix, 0.02);
        check_case("H DLT", matrix);
        h_matrix(matrix, 1e-14);
        check_case("nearly degenerate H DLT", matrix);
    }
    for (trial = 0; trial < 5000; ++trial) {
        for (i = 0; i < N*N; ++i)
            matrix[i] = ldexp(random_double(), (int)(random_double() * 1000));
        check_case("mixed exponents", matrix);
    }
    memset(matrix, 0, sizeof(matrix));
    check_case("zero", matrix);
    for (trial = 0; trial <= N; ++trial) {
        memset(matrix, 0, sizeof(matrix));
        for (i = 0; i < trial; ++i) matrix[N*i+i] = 1.0;
        check_case("finite rank", matrix);
    }
    memset(matrix, 0, sizeof(matrix));
    matrix[0] = nextafter(1e-12, 0.0);
    matrix[N+1] = 1e-12;
    matrix[2*N+2] = nextafter(1e-12, INFINITY);
    check_case("tolerance boundary", matrix);
    memset(matrix, 0, sizeof(matrix));
    for (i = 0; i < N; ++i) matrix[N*i+i] = 1.0;
    matrix[0] = 1e-12; matrix[1] = 1e308;
    check_case("finite overflow during normalisation", matrix);
    random_matrix(matrix, 1.0); matrix[2*N+3] = NAN;
    check_case("nan", matrix);
    random_matrix(matrix, 1.0); matrix[4*N+5] = INFINITY;
    check_case("infinity", matrix);
    random_matrix(matrix, 1.0); matrix[0] = -0.0; matrix[1] = 0.0;
    check_case("signed zero", matrix);
    puts("nullspace9 matches generic 9x9 nullspace");
    return 0;
}
