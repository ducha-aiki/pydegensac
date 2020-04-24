#undef __STRICT_ANSI__
#include <math.h>
//#include <stdio.h>

#include "utools.h"

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

  for (j = 0; j < len; j++)
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


void cov_mat(double *Cv, const double * Z, int len, int siz)
{
   int i, j, k, lenM = len * siz;
   double val;

   for (i=0; i<siz; i++)
      for (j=0; j<=i; j++)
      {
         val = 0;
         for (k=0; k< lenM; k+=siz)
            val += Z[k+i] * Z[k+j];
         Cv[siz*i + j] = val;
         Cv[i + siz*j] = val;
      }
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

