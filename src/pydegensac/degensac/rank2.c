/* Rank-2 projection for a row-major 3x3 F.
   Return one on any unsafe input so a caller can retain its LAPACK fallback. */
#include "lapwrap.h"
#include <float.h>
#include <math.h>
#include <string.h>

static int jacobi_smallest3(const double in[9], double vmin[3])
{
    double b[3][3], v[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
    int i, j, k, step, order[3];

    for (i = 0; i < 3; ++i)
        for (j = i; j < 3; ++j) {
            double value = 0.0;
            for (k = 0; k < 3; ++k) value += in[k * 3 + i] * in[k * 3 + j];
            b[i][j] = b[j][i] = value;
        }

    for (step = 0; step < 32; ++step) {
        int p = 0, q = 1;
        double off, diagonal;
        if (fabs(b[0][2]) > fabs(b[p][q])) { p = 0; q = 2; }
        if (fabs(b[1][2]) > fabs(b[p][q])) { p = 1; q = 2; }
        off = fabs(b[p][q]);
        diagonal = fabs(b[0][0]) + fabs(b[1][1]) + fabs(b[2][2]);
        if (off <= 2.0 * DBL_EPSILON * (diagonal + 1.0)) break;
        {
            double app = b[p][p], aqq = b[q][q], apq = b[p][q];
            double tau = (aqq - app) / (2.0 * apq);
            double t = copysign(1.0, tau) / (fabs(tau) + sqrt(1.0 + tau * tau));
            double c = 1.0 / sqrt(1.0 + t * t), s = t * c;
            for (k = 0; k < 3; ++k) {
                double bkp, bkq;
                if (k == p || k == q) continue;
                bkp = b[k][p]; bkq = b[k][q];
                b[k][p] = b[p][k] = c * bkp - s * bkq;
                b[k][q] = b[q][k] = s * bkp + c * bkq;
            }
            b[p][p] = app - t * apq;
            b[q][q] = aqq + t * apq;
            b[p][q] = b[q][p] = 0.0;
            for (k = 0; k < 3; ++k) {
                double vkp = v[k][p], vkq = v[k][q];
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    if (step == 32) return 1;
    order[0] = 0; order[1] = 1; order[2] = 2;
    for (i = 0; i < 2; ++i)
        for (j = i + 1; j < 3; ++j)
            if (b[order[j]][order[j]] > b[order[i]][order[i]]) {
                int t = order[i]; order[i] = order[j]; order[j] = t;
            }
    /* Forming F'F squares the condition number.  The smallest right vector
       is not trustworthy when the second singular value is below 4*eps^1/4 of
       the first, or when the two smallest singular directions are inseparable. */
    if (!(b[order[1]][order[1]] > 16.0 * sqrt(DBL_EPSILON) * b[order[0]][order[0]] &&
          b[order[0]][order[0]] > 0.0)) return 1;
    if (!((b[order[1]][order[1]] - b[order[2]][order[2]]) >
          16.0 * sqrt(DBL_EPSILON) * b[order[0]][order[0]])) return 1;
    for (i = 0; i < 3; ++i) {
        if (!isfinite(v[i][order[2]])) return 1;
        vmin[i] = v[i][order[2]];
    }
    return 0;
}

/* On success, overwrite F with its closest Frobenius-norm rank-2 matrix. */
int rank2_project3(double *f)
{
    double a[9], out[9], v[3], fv[3], scale = 0.0;
    int i, j;
    for (i = 0; i < 9; ++i) {
        double x = fabs(f[i]);
        if (!isfinite(x)) return 1;
        if (x > scale) scale = x;
    }
    if (scale == 0.0) return 0;
    for (i = 0; i < 9; ++i) a[i] = f[i] / scale;
    if (jacobi_smallest3(a, v) != 0) return 1;
    for (i = 0; i < 3; ++i) {
        fv[i] = 0.0;
        for (j = 0; j < 3; ++j) fv[i] += a[3 * i + j] * v[j];
    }
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) {
            out[3 * i + j] = (a[3 * i + j] - fv[i] * v[j]) * scale;
            if (!isfinite(out[3 * i + j])) return 1;
        }
    memcpy(f, out, sizeof(out));
    return 0;
}
