#ifndef __UTOOLS_H_
#define __UTOOLS_H_

#define _f1 (*F)
#define _f2 (*(F+1))
#define _f3 (*(F+2))
#define _f4 (*(F+3))
#define _f5 (*(F+4))
#define _f6 (*(F+5))
#define _f7 (*(F+6))
#define _f8 (*(F+7))
#define _f9 (*(F+8))

#define crossprod(out,a,b) crossprod_st(out,a,b,1)
#define crossprodT(out,a,b) crossprod_st(out,a,b,3)

void normu (const double *u, const int * inl, int len,
            double *A1, double *A2);

void denormF (double *F, double *A1, double *A2);

void denormH (double *F, double *A1, double *A2);

void scalmul (double *data, double m, int len, int step) ;

int nullspace(double *matrix, double *nullspace, int n, int * buffer);

/* Fixed 9x9 row-major elimination. matrix is overwritten; basis needs up
   to 81 doubles and buffer needs 18 ints, as in nullspace(..., 9, ...). */
int nullspace9(double *matrix, double *basis, int *buffer);

/* Seven-point specialization. matrix is 9x9 row-major scratch with its last
   two rows zero-padded. Basis/rank and the 18-entry pivot buffer match
   nullspace(..., 9, ...); basis needs room for up to 81 doubles. */
int nullspace7x9(double *matrix, double *basis, int *buffer);

void cov_mat(double *Cv, const double * Z, int len, int siz);

void crossprod_st(double *out, const double *a, const double *b, int st);

double det3 (const double *A);

#endif /* __UTOOLS_H_ */
