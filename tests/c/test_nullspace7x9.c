/* Exact regression gate for the seven-correspondence padded nullspace path. */
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

/* utools.c contains cov_mat as well as nullspace. This harness never reaches
   the BLAS branch, but supplies the symbol so it can link without a BLAS
   dependency. */
void dsyrk_(char *uplo, char *trans, ptrdiff_t *n, ptrdiff_t *k,
            double *alpha, double *a, ptrdiff_t *lda, double *beta,
            double *c, ptrdiff_t *ldc)
{
    (void)uplo; (void)trans; (void)n; (void)k; (void)alpha; (void)a;
    (void)lda; (void)beta; (void)c; (void)ldc;
    abort();
}

static uint64_t state = UINT64_C(0x53a9f481d2c7b695);

typedef struct {
    unsigned char before[16];
    double values[N*N];
    unsigned char after[16];
} GuardedBasis;

static double random_double(void)
{
    state = state * UINT64_C(6364136223846793005) + 1;
    return ((state >> 11) * 0x1.0p-53) * 2.0 - 1.0;
}

static void zero_padding(double *matrix)
{
    int i;
    for (i = 7 * N; i < N * N; ++i)
        matrix[i] = 0.0;
}

static void random_matrix(double *matrix, double scale)
{
    int i;
    for (i = 0; i < 7 * N; ++i)
        matrix[i] = scale * random_double();
    zero_padding(matrix);
}

/* One row per correspondence, in the same ordering as sampleF7. */
static void correspondence_matrix(double *matrix, double noise)
{
    int row;
    for (row = 0; row < 7; ++row) {
        double x1 = random_double();
        double y1 = random_double();
        double x2 = x1 + 0.3 + noise * random_double();
        double y2 = y1 - 0.2 + noise * random_double();
        double values[N] = {x2*x1, x2*y1, x2, y2*x1, y2*y1, y2, x1, y1, 1.0};
        memcpy(matrix + row * N, values, sizeof(values));
    }
    zero_padding(matrix);
}

static void check_case(const char *name, const double *source)
{
    double generic_matrix[N*N], fast_matrix[N*N];
    double generic_basis[N*N], fast_basis[N*N];
    GuardedBasis generic_guard, fast_guard;
    int generic_buffer[2*N], fast_buffer[2*N];
    int generic_rank, fast_rank, i;

    memcpy(generic_matrix, source, sizeof(generic_matrix));
    memcpy(fast_matrix, source, sizeof(fast_matrix));
    memset(&generic_guard, SENTINEL, sizeof(generic_guard));
    memset(&fast_guard, SENTINEL, sizeof(fast_guard));
    memset(generic_buffer, SENTINEL, sizeof(generic_buffer));
    memset(fast_buffer, SENTINEL, sizeof(fast_buffer));

    generic_rank = nullspace(generic_matrix, generic_guard.values, N,
                             generic_buffer);
    fast_rank = nullspace7x9(fast_matrix, fast_guard.values, fast_buffer);
    memcpy(generic_basis, generic_guard.values, sizeof(generic_basis));
    memcpy(fast_basis, fast_guard.values, sizeof(fast_basis));

    if (generic_rank != fast_rank ||
        memcmp(generic_basis, fast_basis, sizeof(generic_basis)) != 0 ||
        memcmp(generic_buffer, fast_buffer, sizeof(generic_buffer)) != 0) {
        fprintf(stderr, "%s: rank/basis/pivot buffer differs (%d vs %d)\n",
                name, generic_rank, fast_rank);
        assert(0);
    }
    for (i = 0; i < (int)sizeof(generic_guard.before); ++i) {
        assert(generic_guard.before[i] == SENTINEL);
        assert(fast_guard.before[i] == SENTINEL);
    }
    for (i = 0; i < (int)sizeof(generic_guard.after); ++i) {
        assert(generic_guard.after[i] == SENTINEL);
        assert(fast_guard.after[i] == SENTINEL);
    }
}

int main(void)
{
    double matrix[N*N];
    int i, trial;

    for (trial = 0; trial < 4000; ++trial) {
        random_matrix(matrix, 1.0);
        check_case("random", matrix);
        correspondence_matrix(matrix, 0.02);
        check_case("geometric", matrix);
        correspondence_matrix(matrix, 1e-14);
        check_case("nearly planar", matrix);
    }
    for (trial = 0; trial < 4000; ++trial) {
        for (i = 0; i < 7 * N; ++i)
            matrix[i] = ldexp(random_double(),
                              (int)(random_double() * 1000.0));
        zero_padding(matrix);
        check_case("mixed exponents", matrix);
    }

    correspondence_matrix(matrix, 0.0);
    check_case("exact planar", matrix);
    for (i = 0; i < 7; ++i)
        memcpy(matrix + i*N, matrix, N * sizeof(double));
    check_case("repeated rows", matrix);
    memset(matrix, 0, sizeof(matrix));
    check_case("zero", matrix);
    random_matrix(matrix, 1e-150);
    check_case("tiny", matrix);
    random_matrix(matrix, 1e150);
    check_case("large", matrix);

    memset(matrix, 0, sizeof(matrix));
    for (i = 0; i < 7; ++i)
        matrix[i*N + i] = 1.0;
    matrix[0] = 1e-12;
    matrix[1] = 1e308;
    check_case("finite overflow during normalisation", matrix);
    random_matrix(matrix, 1.0);
    matrix[2*N + 3] = NAN;
    check_case("nan", matrix);
    random_matrix(matrix, 1.0);
    matrix[4*N + 5] = INFINITY;
    check_case("infinity", matrix);
    random_matrix(matrix, 1.0);
    matrix[0] = -0.0;
    matrix[1] = 0.0;
    check_case("signed zero", matrix);

    /* Exercise every finite rank and the exact tolerance boundary. */
    for (trial = 0; trial <= 7; ++trial) {
        memset(matrix, 0, sizeof(matrix));
        for (i = 0; i < trial; ++i)
            matrix[i*N + i] = 1.0;
        check_case("finite rank", matrix);
    }
    memset(matrix, 0, sizeof(matrix));
    matrix[0] = nextafter(1e-12, 0.0);
    matrix[N + 1] = 1e-12;
    matrix[2*N + 2] = nextafter(1e-12, INFINITY);
    check_case("tolerance boundary", matrix);
    puts("nullspace7x9 matches generic padded nullspace");
    return 0;
}
