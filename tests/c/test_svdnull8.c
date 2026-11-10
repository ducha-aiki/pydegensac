/* Regression gate for the eight-correspondence nullspace shortcut.

   F's eight-point LO builds a 9 by 8 matrix, and consumes only column 8 of
   the 9 by 9 left factor returned by svduv.  svdnull8 must produce that
   column exactly, including signed zero and NaN payloads, while avoiding the
   work needed only for the other singular vectors. */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "matutl.h"

#define ROWS 9
#define COLS 8

static uint64_t state = UINT64_C(0x9e3779b97f4a7c15);

static double uniform(void)
{
    state = state * UINT64_C(6364136223846793005) + 1;
    return ((state >> 11) * 0x1.0p-53) * 2.0 - 1.0;
}

static void random_matrix(double *a, double scale)
{
    int i;
    for (i = 0; i < ROWS * COLS; ++i)
        a[i] = scale * uniform();
}

/* Match lin_fm's 9-by-8 layout: each column is one correspondence. */
static void correspondence_matrix(double *a, double noise, double scale)
{
    int i;
    for (i = 0; i < COLS; ++i) {
        double x1 = uniform();
        double y1 = uniform();
        double x2 = x1 + 0.3 + noise * uniform();
        double y2 = y1 - 0.2 + noise * uniform();
        double row[ROWS] = {x2*x1, x2*y1, x2, y2*x1, y2*y1, y2, x1, y1, 1.0};
        int j;
        for (j = 0; j < ROWS; ++j)
            a[j * COLS + i] = scale * row[j];
    }
}

static void repeated_column_matrix(double *a)
{
    int i, j;
    for (j = 0; j < ROWS; ++j)
        for (i = 0; i < COLS; ++i)
            a[j * COLS + i] = (i < 4) ? (double)(j - 4) : (double)(2 * j - 3);
}

static void compare_against_full(const char *name, const double *input,
                                 int check_residual)
{
    double full_a[ROWS * COLS], fast_a[ROWS * COLS];
    double singular[COLS], full_u[ROWS * ROWS], full_v[COLS * COLS];
    double fast_u[ROWS * ROWS];
    int i;

    memcpy(full_a, input, sizeof(full_a));
    memcpy(fast_a, input, sizeof(fast_a));
    assert(svduv(singular, full_a, full_u, ROWS, full_v, COLS) == 0);
    assert(svdnull8(fast_a, fast_u) == 0);
    for (i = 0; i < ROWS; ++i) {
        if (memcmp(&full_u[i * ROWS + COLS], &fast_u[i * ROWS + COLS],
                   sizeof(double)) != 0) {
            fprintf(stderr, "%s: null vector differs at %d\n", name, i);
            assert(0);
        }
    }

    if (check_residual) {
        double residual2 = 0.0, matrix2 = 0.0, vector2 = 0.0;
        int column;
        for (column = 0; column < COLS; ++column) {
            double dot = 0.0;
            for (i = 0; i < ROWS; ++i) {
                double x = input[i * COLS + column];
                dot += x * fast_u[i * ROWS + COLS];
                matrix2 += x * x;
            }
            residual2 += dot * dot;
        }
        for (i = 0; i < ROWS; ++i)
            vector2 += fast_u[i * ROWS + COLS] * fast_u[i * ROWS + COLS];
        assert(sqrt(residual2) <= 1e-10 * sqrt(matrix2 * vector2));
    }
}

int main(void)
{
    double a[ROWS * COLS];
    int i, trial;

    for (trial = 0; trial < 80; ++trial) {
        random_matrix(a, 1.0);
        compare_against_full("random", a, 1);
    }
    for (trial = 0; trial < 40; ++trial) {
        correspondence_matrix(a, 0.02, 1.0);
        compare_against_full("correspondence", a, 1);
        for (i = 0; i < COLS; ++i) {
            double weight = pow(10.0, 6.0 * uniform());
            int row;
            for (row = 0; row < ROWS; ++row)
                a[row * COLS + i] *= weight;
        }
        compare_against_full("weighted correspondence", a, 1);
        correspondence_matrix(a, 0.0, 1.0);
        compare_against_full("planar", a, 0);
        correspondence_matrix(a, 1e-12, 1.0);
        compare_against_full("nearly planar", a, 1);
    }
    repeated_column_matrix(a);
    compare_against_full("repeated columns", a, 0);
    memset(a, 0, sizeof(a));
    compare_against_full("zero", a, 0);
    for (i = 0; i < ROWS * COLS; ++i)
        a[i] = (i & 1) ? -0.0 : 0.0;
    compare_against_full("signed zero", a, 0);
    random_matrix(a, 1e-150);
    compare_against_full("tiny", a, 0);
    random_matrix(a, 1e150);
    compare_against_full("large", a, 0);
    random_matrix(a, 1.0);
    a[0] = NAN;
    compare_against_full("nan", a, 0);
    random_matrix(a, 1.0);
    a[ROWS * COLS - 1] = INFINITY;
    compare_against_full("infinity", a, 0);
    puts("svdnull8 matches svduv's consumed nullspace column");
    return 0;
}
