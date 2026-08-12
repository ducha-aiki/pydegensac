#undef __STRICT_ANSI__
#include <math.h>
//#include <stdio.h>

#include "utools.h"
#include "lapwrap.h"

/* Length from which cov_mat hands Z^T Z to BLAS dsyrk instead of doing it
   itself. In isolation dsyrk already wins at len=10 (M1/Accelerate, siz=9, ns
   per call: len 8 -- 213 by hand vs 216 dsyrk; 10 -- 257 vs 164; 64 -- 1543 vs
   244; 256 -- 6519 vs 516), but in the estimators a threshold that low costs
   F 3% while gaining H 25%, because F's small-len calls pay the call overhead
   more often than they save. Measured end to end, estimator calls per 10 s:

     threshold      F      H
     no dsyrk     297   3880
     >= 10        288   4860
     >= 32        306   5394
     >= 64        306   5346

   32 is the only setting that wins on both. */
#define COV_BLAS_MIN 32

void normu (const double *u, const int * inl, int len, 
           double *A1, double *A2)
{
  int i,j;
  double a,b;
  const double *p = u;

  for (j = 0; j < 3; j++)
    {
      A1[j] = 0; A2[j] = 0;
    }

  for (j = 0; j < len; j++)
    {
      u = p+ 6*inl[j];
      A1[1] += u[0]; A1[2] += u[1]; 
      A2[1] += u[3]; A2[2] += u[4]; 
    }

  if (len > 0)
    for (i = 1; i < 3; i++)
      {
        A1[i] /= len; A2[i] /= len;
      }

  /* Two partial sums per image rather than one: each mean-distance
     accumulator was a serial chain with a square root in it, so the loop ran
     at sqrt-plus-add latency on work that is independent per correspondence.
     Reassociation changes rounding. */
  {
    double s1a = 0, s1b = 0, s2a = 0, s2b = 0;
    int n2 = len & ~1;
    for (j = 0; j < n2; j += 2)
      {
        const double *q = p + 6*inl[j];
        const double *w = p + 6*inl[j+1];
        double qa = q[0] - A1[1], qb = q[1] - A1[2];
        double wa = w[0] - A1[1], wb = w[1] - A1[2];
        s1a += sqrt(qa*qa + qb*qb);
        s1b += sqrt(wa*wa + wb*wb);
        qa = q[3] - A2[1]; qb = q[4] - A2[2];
        wa = w[3] - A2[1]; wb = w[4] - A2[2];
        s2a += sqrt(qa*qa + qb*qb);
        s2b += sqrt(wa*wa + wb*wb);
      }
    A1[0] = s1a + s1b;
    A2[0] = s2a + s2b;
    j = n2;
  }

  for (; j < len; j++)
    {
      u = p+ 6*inl[j];
      a = u[0] - A1[1];
      b = u[1] - A1[2];
      A1[0] += sqrt(a*a + b*b);

      a = u[3] - A2[1];
      b = u[4] - A2[2];
      A2[0] += sqrt(a*a + b*b);
    }

  if (A1[0] != 0)
      A1[0] = len * sqrt(2) / A1[0];
  if (A2[0] != 0)
      A2[0] = len * sqrt(2) / A2[0];

   A1[1] *= -A1[0]; A1[2] *= -A1[0];
   A2[1] *= -A2[0]; A2[2] *= -A2[0];
}

void denormF (double *F, double *A1, double *A2)
{
  double r, x, y;
  
  r = A2[0]; x = A2[1]; y = A2[2];
  _f7 += x * _f1 + y*_f4;
  _f8 += x * _f2 + y*_f5;
  _f9 += x * _f3 + y*_f6;
  _f1 *= r; _f2 *= r; _f3 *= r;
  _f4 *= r; _f5 *= r; _f6 *= r;

  r = A1[0]; x = A1[1]; y = A1[2];
  _f3 += x * _f1 + y*_f2;
  _f6 += x * _f4 + y*_f5;
  _f9 += x * _f7 + y*_f8;
  _f1 *= r; _f4 *= r; _f7 *= r;
  _f2 *= r; _f5 *= r; _f8 *= r;
}

void denormH (double *F, double *A1, double *A2)
{
  double r, x, y;
  int i;  

  r = A2[0]; x = A2[1]; y = A2[2];
  _f7 += x * _f1 + y*_f4;
  _f8 += x * _f2 + y*_f5;
  _f9 += x * _f3 + y*_f6;
  _f1 *= r; _f2 *= r; _f3 *= r;
  _f4 *= r; _f5 *= r; _f6 *= r;

  r = 1/A1[0]; x = -A1[1] * r; y = -A1[2] * r;

  for (i = 0; i < 9; i+=3)
    {
      F[i]   = r * F[i]   + x * F[i+2];
      F[i+1] = r * F[i+1] + y * F[i+2];
    }

}

void scalmul (double *data, double m, int len, int step) 
{int i; for (i =0; i < len; i++, data += step) *data *= m;}

int nullspace(double *matrix, double *nullspace, int n, int * buffer) /* Expects matrix to be stored row-wise */
     /* buffer size 2*n*sizeof(int) */
{
   int *pnopivot=buffer, nonpivot=0;
   int *ppivot=buffer+n;
   int i, j, k, l, ptr, max;
   double pivot, t;
   double tol=1e-12;
   
   ptr = 0;
   i = 0;
   for (j=0;j<n;j++)
   {
      /* find pivot, start with diagonal element */
      pivot = fabs(matrix[n*i+j]); max = i;
      for (k=i+1; k<n; k++)
      {
         t = fabs(matrix[n*k+j]);
         if (pivot<t) { pivot=t; max=k; }
      }
      if (pivot<tol)
      {
         *(pnopivot++) = j; nonpivot++;
         /* negligible column, zero out */
         for (k=i;k<n;k++) matrix[n*k+j]=0;
      } else {
         *(ppivot++) = j;
         /* swap rows i <-> max */
         for (k=j; k<n; k++)
         {
            t = matrix[i*n+k]; 
            matrix[i*n+k] = matrix[max*n+k];
            matrix[max*n+k]=t;
         }
         pivot = matrix[i*n+j];
         /* divide the pivot row by the pivot element. */
         for (k=j; k<n; k++)
            matrix[i*n+k] /= pivot;

         /* Subtract multiples of the pivot row from all the other rows. */
         for (k=0; k<i; k++)
         {
            pivot = -matrix[k*n+j];
            for (l=j; l<n; l++)
               matrix[k*n+l] += pivot*matrix[i*n+l];
         }
         
         for (k=i+1; k<n; k++)
         {
            pivot = matrix[k*n+j];
            for (l=j; l<n; l++)
               matrix[k*n+l] -= pivot*matrix[i*n+l];
         }
         i++;
      }
   }
   
   /* initialize null space vectors */
   for (k=0;k<nonpivot;k++)
   {      
      j=buffer[k];
      /* copy nonpivot -column above diagonal */
      for (l=0;l<n-nonpivot;l++)
         nullspace[k*n+buffer[n+l]]=-matrix[l*n+j];
      
      for (l=0;l<nonpivot;l++)
         nullspace[k*n+buffer[l]]=(j==buffer[l])?1:0;
   }
   /* number of nullspace vectors */
   return nonpivot;
}


/* Cv = Z^T Z for Z of shape len x siz, row-major; both triangles filled.

   One pass over the points, accumulating every unique entry at once. The
   previous form ran one pass per entry -- 45 of them at siz=9 -- and each was
   a single serial FP accumulation chain, so it was latency-bound rather than
   throughput-bound. The hot caller is exp_ranH's per-iteration MCE solve
   (len=10), which runs this on every RANSAC iteration: ~45 dependent chains
   against ~45 independent accumulators.

   Summation over points still runs in ascending order; what changes is that
   the partial sums live in a different order of operations, so rounding can
   differ in the last bits. Deliberate -- see docs/superpowers/specs/
   2026-08-12-covmat-inlidxs-perf-design.md. */
void cov_mat(double *Cv, const double * Z, int len, int siz)
{
   int i, j, k, lenM = len * siz;

   for (i=0; i<siz*siz; i++)
      Cv[i] = 0;

   if (len >= COV_BLAS_MIN)
   {
      /* Fortran reads the row-major len x siz array Z as a column-major
         siz x len matrix A, so A*A^T ("N", no transpose) is the Z^T Z we
         want. The result is symmetric, so the row-major/column-major
         distinction does not matter for Cv -- but dsyrk writes only one
         triangle, hence the mirror below. Integer widths follow the rest of
         this library: lapack_int is ptrdiff_t and the vendor BLAS reads the
         low half, which is correct for these small positive values on any
         little-endian target. */
      lapack_int n = siz, kk = len, lda = siz, ldc = siz;
      double alpha = 1.0, beta = 0.0;
      dsyrk_("L", "N", &n, &kk, &alpha, (double *) Z, &lda, &beta, Cv, &ldc);
      for (i=0; i<siz; i++)
         for (j=0; j<i; j++)
            Cv[siz*i + j] = Cv[i + siz*j];
      return;
   }

   for (k=0; k<lenM; k+=siz)
      for (i=0; i<siz; i++)
      {
         const double zi = Z[k+i];
         for (j=0; j<=i; j++)
            Cv[siz*i + j] += zi * Z[k+j];
      }

   for (i=0; i<siz; i++)
      for (j=0; j<i; j++)
         Cv[i + siz*j] = Cv[siz*i + j];
}


void crossprod_st(double *out, const double *a, const double *b, int st)
{
   int st2 = 2 * st;
   *out   = a[st]*b[st2] - a[st2]*b[st];
   out[1] = a[st2]*b[0]  - a[0]*b[st2];
   out[2] = a[0]*b[st]   - a[st]*b[0];
}


double det3 (const double *A)
{
   double r;
   r = (A[0]*A[4]*A[8] + A[2]*A[3]*A[7] + A[1]*A[5]*A[6]);
   r -=(A[2]*A[4]*A[6] + A[0]*A[5]*A[7] + A[1]*A[3]*A[8]);
   return(r);
}

