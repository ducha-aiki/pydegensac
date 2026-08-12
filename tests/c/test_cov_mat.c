/* Checks cov_mat against a naive reference over many shapes.

   cov_mat is not reachable from Python, so this is a plain C harness rather
   than a pytest case. It exists because the fast form accumulates in a
   different order than the obvious one, and "different order" must stay
   "same answer to within rounding", not "same answer usually".

   Build and run:
     mkdir -p tests/bin && cc -O2 -I src/pydegensac/degensac \
       tests/c/test_cov_mat.c src/pydegensac/degensac/utools.c \
       -o tests/bin/test_cov_mat && ./tests/bin/test_cov_mat            */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cov_mat(double *Cv, const double * Z, int len, int siz);

/* The pre-2026-08-12 implementation, kept as the reference. */
static void reference(double *Cv, const double *Z, int len, int siz) {
    int i, j, k;
    for (i = 0; i < siz; i++)
        for (j = 0; j <= i; j++) {
            double val = 0;
            for (k = 0; k < len*siz; k += siz) val += Z[k+i] * Z[k+j];
            Cv[siz*i + j] = val;
            Cv[i + siz*j] = val;
        }
}

int main(void) {
    /* siz 9 is the common case; nullsize (exp_ranH.c:750) makes it variable,
       and len runs from a handful (the per-iteration MCE solve) to the full
       inlier count. */
    const int sizes[] = {2, 3, 7, 8, 9, 10};
    const int lens[] = {1, 2, 8, 10, 33, 64, 257};
    int si, li, t, i;
    double worst = 0.0;
    srand(1);
    for (si = 0; si < 6; si++)
        for (li = 0; li < 7; li++)
            for (t = 0; t < 5; t++) {
                int siz = sizes[si], len = lens[li];
                double *Z = (double *) malloc(sizeof(double)*len*siz);
                double *A = (double *) malloc(sizeof(double)*siz*siz);
                double *B = (double *) malloc(sizeof(double)*siz*siz);
                for (i = 0; i < len*siz; i++)
                    Z[i] = (double)rand()/RAND_MAX*2.0 - 1.0;
                reference(A, Z, len, siz);
                cov_mat(B, Z, len, siz);
                for (i = 0; i < siz*siz; i++) {
                    double denom = fabs(A[i]) > 1.0 ? fabs(A[i]) : 1.0;
                    double rel = fabs(A[i] - B[i]) / denom;
                    if (rel > worst) worst = rel;
                }
                free(Z); free(A); free(B);
            }
    printf("worst relative difference vs reference: %.3g\n", worst);
    if (worst > 1e-12) { printf("FAIL\n"); return 1; }
    printf("PASS\n");
    return 0;
}
