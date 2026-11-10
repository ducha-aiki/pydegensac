/*
 * Fixed-size symmetric eigensolver for the 9x9 normal matrices.
 *
 * Adapted from mpmath 1.3.0, matrices/eigen_symmetric.py:
 * r_sy_tridiag and tridiag_eigen, Copyright 2013 Timo Hartmann.
 * https://github.com/mpmath/mpmath/blob/1.3.0/mpmath/matrices/eigen_symmetric.py
 * Those routines derive from EISPACK tred2.f and imtql2.f.
 *
 * Copyright (c) 2005-2021 Fredrik Johansson and mpmath contributors
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *   a. Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 *   b. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   c. Neither the name of the copyright holder nor the names of its
 *      contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 */
#include "lapwrap.h"
#include <float.h>
#include <math.h>
#include <string.h>
#define N 9
#define MAXIT 96

static void swap_col(double v[N][N], int a, int b)
{
    int r;
    for (r = 0; r < N; ++r) {
        double t = v[r][a]; v[r][a] = v[r][b]; v[r][b] = t;
    }
}

/* Householder tridiagonalization, adapted from mpmath r_sy_tridiag. */
static void tred2(double v[N][N], double d[N], double e[N])
{
    int i, j, k;
    for (i = N - 1; i > 0; --i) {
        double scale = 0.0, h = 0.0;
        for (k = 0; k < i; ++k) scale += fabs(v[k][i]);
        if (i == 1 || scale == 0.0 || !isfinite(1.0 / scale)) {
            e[i] = v[i - 1][i];
            d[i] = 0.0;
        }
        else {
            double f, g, hh;
            for (k = 0; k < i; ++k) { v[k][i] /= scale; h += v[k][i] * v[k][i]; }
            f = v[i - 1][i]; g = sqrt(h); if (f > 0.0) g = -g;
            e[i] = scale * g; h -= f * g; v[i - 1][i] = f - g;
            f = 0.0;
            for (j = 0; j < i; ++j) {
                v[i][j] = v[j][i] / h;
                g = 0.0;
                for (k = 0; k <= j; ++k) g += v[k][j] * v[k][i];
                for (k = j + 1; k < i; ++k) g += v[j][k] * v[k][i];
                e[j] = g / h;
                f += e[j] * v[j][i];
            }
            hh = f / (2.0 * h);
            for (j = 0; j < i; ++j) {
                f = v[j][i]; g = e[j] - hh * f; e[j] = g;
                for (k = 0; k <= j; ++k)
                    v[k][j] -= f * e[k] + g * v[k][i];
            }
            d[i] = h;
        }
    }
    d[0] = e[0] = 0.0;
    for (i = 0; i < N; ++i) {
        if (d[i] != 0.0) for (j = 0; j < i; ++j) {
            double g = 0.0;
            for (k = 0; k < i; ++k) g += v[i][k] * v[k][j];
            for (k = 0; k < i; ++k) v[k][j] -= g * v[k][i];
        }
        d[i] = v[i][i]; v[i][i] = 1.0;
        for (j = 0; j < i; ++j) v[j][i] = v[i][j] = 0.0;
    }
}

/* Between 2^-500 and 2^500, each nonzero square stays normal and their
   sum cannot overflow. Other ranges retain hypot's scaling safeguards. */
static double hypot_normal(double x, double y)
{
    double ax = fabs(x), ay = fabs(y);
    if ((ax == 0.0 || (ax >= 3.054936363499605e-151 && ax <= 3.273390607896142e150)) &&
        (ay == 0.0 || (ay >= 3.054936363499605e-151 && ay <= 3.273390607896142e150)))
        return sqrt(x*x + y*y);
    return hypot(x, y);
}

/* Implicit shifted QL.  A cap permits safe fallback instead of partial roots. */
static int tqli(double v[N][N], double d[N], double e[N])
{
    int i, k, l, m;
    for (i = 1; i < N; ++i) e[i - 1] = e[i];
    e[N - 1] = 0.0;
    for (l = 0; l < N; ++l) {
        int its = 0;
        do {
            for (m = l; m < N - 1; ++m)
                if (fabs(e[m]) <= DBL_EPSILON * (fabs(d[m]) + fabs(d[m + 1]))) break;
            if (m != l) {
                double b, c, f, g, p, r, s;
                if (++its > MAXIT) return 1;
                g = (d[l + 1] - d[l]) / (2.0 * e[l]); r = hypot_normal(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + copysign(r, g)); s = c = 1.0; p = 0.0;
                for (i = m - 1; i >= l; --i) {
                    f = s * e[i]; b = c * e[i]; r = hypot_normal(f, g); e[i + 1] = r;
                    if (r == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
                    s = f / r; c = g / r; g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b; p = s * r;
                    d[i + 1] = g + p; g = c * r - b;
                    for (k = 0; k < N; ++k) {
                        f = v[k][i + 1]; v[k][i + 1] = s * v[k][i] + c * f;
                        v[k][i] = c * v[k][i] - s * f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p; e[l] = g; e[m] = 0.0;
            }
        } while (m != l);
    }
    return 0;
}

/* lap_eig-compatible success ABI: upper triangle of column-major a on input;
   column-major eigenvectors and ascending ev on output.  Failure preserves a. */
int eig9_symmetric(double *a, double *ev)
{
    double v[N][N], d[N], e[N], scale = 0.0;
    int i, j;
    for (j = 0; j < N; ++j) for (i = 0; i <= j; ++i) {
        double x = a[i + N * j];
        if (!isfinite(x)) return 1;
        if (fabs(x) > scale) scale = fabs(x);
    }
    if (scale == 0.0) {
        memset(a, 0, N * N * sizeof(*a));
        for (i = 0; i < N; ++i) { a[i + N * i] = 1.0; ev[i] = 0.0; }
        return 0;
    }
    if (scale < DBL_MIN) return 1;
    for (i = 0; i < N; ++i) for (j = i; j < N; ++j)
        v[i][j] = v[j][i] = a[i + N * j] / scale;
    tred2(v, d, e);
    if (tqli(v, d, e)) return 1;
    for (i = 0; i < N - 1; ++i) for (j = i + 1; j < N; ++j)
        if (d[j] < d[i]) { double t = d[i]; d[i] = d[j]; d[j] = t; swap_col(v, i, j); }
    /* Complete all failure checks before touching the caller's output. */
    for (i = 0; i < N; ++i) {
        if (!isfinite(d[i] * scale)) return 1;
        for (j = 0; j < N; ++j)
            if (!isfinite(v[j][i])) return 1;
    }
    for (i = 0; i < N; ++i) {
        ev[i] = d[i] * scale;
        for (j = 0; j < N; ++j) a[j + N * i] = v[j][i];
    }
    return 0;
}

static int psd_ldlt(const double d[N], const double e[N], double shift,
                    double norm, double ld[N], double l[N])
{
    int i;
    ld[0] = d[0] - shift;
    if (!(ld[0] > 32.0 * DBL_EPSILON * norm) || !isfinite(ld[0])) return 1;
    l[0] = 0.0;
    for (i = 1; i < N; ++i) {
        l[i] = e[i] / ld[i - 1];
        ld[i] = d[i] - shift - l[i] * e[i];
        if (!(ld[i] > 32.0 * DBL_EPSILON * norm) || !isfinite(l[i]) || !isfinite(ld[i])) return 1;
    }
    return 0;
}

static void psd_solve(const double ld[N], const double l[N], double x[N])
{
    int i;
    for (i = 1; i < N; ++i) x[i] -= l[i] * x[i - 1];
    for (i = 0; i < N; ++i) x[i] /= ld[i];
    for (i = N - 2; i >= 0; --i) x[i] -= l[i + 1] * x[i + 1];
}

/* Negative-eigenvalue count for a tridiagonal T-shift I. A near-zero
   recurrence pivot makes the count unreliable, so force the full solver. */
static int sturm_count(const double d[N], const double e[N], double shift, double norm)
{
    double p = d[0] - shift;
    int i, count = 0;
    if (!isfinite(p) || fabs(p) <= 32.0 * DBL_EPSILON * norm) return -1;
    if (p < 0.0) ++count;
    for (i = 1; i < N; ++i) {
        p = d[i] - shift - e[i] * e[i] / p;
        if (!isfinite(p) || fabs(p) <= 32.0 * DBL_EPSILON * norm) return -1;
        if (p < 0.0) ++count;
    }
    return count;
}

/* For PSD normal matrices, recover only the needed minimum eigenvector.
   Inverse iteration avoids the full QL decomposition. The Sturm counts require
   exactly one eigenvalue near the Rayleigh value and none below it, so a
   converged higher eigenvector cannot be accepted. Ambiguous or slowly
   converging cases leave a/out untouched for the general solver. */
int eig9_psd_smallest(const double *a, double *out)
{
    double q[N][N], d[N], e[N], ld[N], l[N], x[N], y[N], local[N];
    double scale = 0.0, norm, shift, lambda, normx, residual = 0.0, delta;
    int i, j, iter, below, above, negative;
    for (j = 0; j < N; ++j) for (i = 0; i <= j; ++i) {
        double z = a[i + N*j];
        if (!isfinite(z)) return 1;
        if (fabs(z) > scale) scale = fabs(z);
    }
    if (scale == 0.0) { memset(out, 0, N * sizeof(*out)); out[0] = 1.0; return 0; }
    if (scale < DBL_MIN) return 1;
    for (i = 0; i < N; ++i) for (j = i; j < N; ++j)
        q[i][j] = q[j][i] = a[i + N*j] / scale;
    tred2(q, d, e);
    norm = 0.0;
    for (i = 0; i < N; ++i) {
        double row = fabs(d[i]);
        if (i) row += fabs(e[i]);
        if (i + 1 < N) row += fabs(e[i + 1]);
        if (row > norm) norm = row;
    }
    if (!(norm > 0.0) || !isfinite(norm)) return 1;
    delta = 16.0 * sqrt(DBL_EPSILON) * norm;
    /* Reject inputs that are not PSD apart from roundoff. */
    negative = sturm_count(d, e, -delta, norm);
    if (negative != 0) return 1;
    shift = -64.0 * DBL_EPSILON * norm;
    if (psd_ldlt(d, e, shift, norm, ld, l)) return 1;
    for (i = 0; i < N; ++i) x[i] = 1.0 + 0.125 * i;
    for (iter = 0; iter < 20; ++iter) {
        double tr_residual = 0.0;
        psd_solve(ld, l, x);
        normx = 0.0; for (i = 0; i < N; ++i) normx += x[i] * x[i]; normx = sqrt(normx);
        if (!(normx > 0.0) || !isfinite(normx)) return 1;
        normx = 1.0 / normx;
        for (i = 0; i < N; ++i) x[i] *= normx;
        lambda = 0.0;
        for (i = 0; i < N; ++i) {
            double tx = d[i] * x[i];
            if (i) tx += e[i] * x[i - 1];
            if (i + 1 < N) tx += e[i + 1] * x[i + 1];
            lambda += x[i] * tx;
        }
        for (i = 0; i < N; ++i) {
            double tx = d[i] * x[i];
            if (i) tx += e[i] * x[i - 1];
            if (i + 1 < N) tx += e[i + 1] * x[i + 1];
            tx -= lambda * x[i]; tr_residual += tx * tx;
        }
        if (sqrt(tr_residual) <= 32.0 * DBL_EPSILON * norm) break;
    }
    if (iter == 20) return 1;
    if (!(lambda >= -64.0 * DBL_EPSILON * norm) || !isfinite(lambda)) return 1;
    below = sturm_count(d, e, lambda - delta, norm);
    above = sturm_count(d, e, lambda + delta, norm);
    if (below != 0 || above != 1) return 1;
    for (i = 0; i < N; ++i) { y[i] = 0.0; for (j = 0; j < N; ++j) y[i] += q[i][j] * x[j]; }
    normx = 0.0; for (i = 0; i < N; ++i) normx += y[i] * y[i]; normx = sqrt(normx);
    if (!(normx > 0.0) || !isfinite(normx)) return 1;
    for (i = 0; i < N; ++i) local[i] = y[i] / normx;
    for (i = 0; i < N; ++i) {
        double r = 0.0;
        for (j = 0; j < N; ++j) r += (i <= j ? a[i + N*j] : a[j + N*i]) / scale * local[j];
        r -= lambda * local[i]; residual += r * r;
    }
    if (!(sqrt(residual) <= 64.0 * DBL_EPSILON * norm)) return 1;
    memcpy(out, local, sizeof(local));
    return 0;
}
