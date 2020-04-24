#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stddef.h>

#include "../matutls/matutl.h"
#include "utools.h"
#include "lapwrap.h"

#include "Ftools.h"
#define SYMMETRIC_ERROR
#define SYMMETRIC_ERROR_CHECK
#define pit 1.0471975511965967
void lin_fm(const double *u, double *p, const int* inl, const int len)
{
    /* linearizes corresp. with respect to entries of fundamental matrix,
      so that u' F u -> A f */

    const double  *s;
    int i,k,l,pos;

    for (i = 0; i < len; i++)
    {
        s = u + 6*inl[i];
        pos = 0;
        for (k = 0; k < 3; k++)
        {
            for (l = 0; l < 3; l++)
            {
                *(p+pos) = *(s+k+3) * (*(s+l));
                pos += len;
            }
        }
        p++;
    }
}

void slcm(double *A, double *B, double *p)
{
    /* calculates polynomial p in x, so that det(xA + (1-x)B) = 0
      where A,B are [3][3] and p is [4] arrays
      ** CHANGES B to A-B ***
      so finally det(A + (x-1) B) = 0 */

    int i;

    *p = -(b13*b22*b31) + b12*b23*b31 + b13*b21*b32 -
            b11*b23*b32 - b12*b21*b33 + b11*b22*b33;

    *(p+1) = -(a33*b12*b21) + a32*b13*b21 + a33*b11*b22 -
            a31*b13*b22 - a32*b11*b23 + a31*b12*b23 +
            a23*b12*b31 - a22*b13*b31 - a13*b22*b31 +
            3*b13*b22*b31 + a12*b23*b31 - 3*b12*b23*b31 -
            a23*b11*b32 + a21*b13*b32 + a13*b21*b32 -
            3*b13*b21*b32 - a11*b23*b32 + 3*b11*b23*b32 +
            (a22*b11 - a21*b12 - a12*b21 + 3*b12*b21 + a11*b22 -
             3*b11*b22)*b33;

    *(p+2) = -(a21*a33*b12) + a21*a32*b13 +
            a13*a32*b21 - a12*a33*b21 + 2*a33*b12*b21 -
            2*a32*b13*b21 - a13*a31*b22 + a11*a33*b22 -
            2*a33*b11*b22 + 2*a31*b13*b22 + a12*a31*b23 -
            a11*a32*b23 + 2*a32*b11*b23 - 2*a31*b12*b23 +
            2*a13*b22*b31 - 3*b13*b22*b31 - 2*a12*b23*b31 +
            3*b12*b23*b31 + a13*a21*b32 - 2*a21*b13*b32 -
            2*a13*b21*b32 + 3*b13*b21*b32 + 2*a11*b23*b32 -
            3*b11*b23*b32 + a23*
            (-(a32*b11) + a31*b12 + a12*b31 - 2*b12*b31 -
             a11*b32 + 2*b11*b32) +
            (-(a12*a21) + 2*a21*b12 + 2*a12*b21 - 3*b12*b21 -
             2*a11*b22 + 3*b11*b22)*b33 +
            a22*(a33*b11 - a31*b13 - a13*b31 + 2*b13*b31 +
                 a11*b33 - 2*b11*b33);

    for (i=0; i < 9; i++)
        B[i] = A[i] - B[i];

    *(p+3) =-(b13*b22*b31) + b12*b23*b31 + b13*b21*b32 -
            b11*b23*b32 - b12*b21*b33 + b11*b22*b33;
}

void FDs (const double *u, const double *F, double *p, int len)
{
    double rx, ry, rwc, ryc, rxc, r;
    double a,b;
    int i;

    for (i=1; i<=len; i++)
    {
        rxc = _f1 * u4 + _f4 * u5 + _f7;
        ryc = _f2 * u4 + _f5 * u5 + _f8;
        rwc = _f3 * u4 + _f6 * u5 + _f9;
        r =(u1 * rxc + u2 * ryc + rwc);
        rx = _f1 * u1 + _f2 * u2 + _f3;
        ry = _f4 * u1 + _f5 * u2 + _f6;
        *p = r*r / (rxc*rxc + ryc*ryc + rx*rx + ry*ry); //original, Sampson`s error
        p ++;
        u += 6;
    }
}

void FDsidx (const double *mu, const double *F, double *p, int len,  int *idx, int siz)
{
    double rx, ry, rwc, ryc, rxc, r;
    const double *u;
    int i, mi;

    for (mi=0; mi<siz; mi++)
    {
        i = idx[mi];
        u = mu + 6*i;

        rxc = _f1 * u4 + _f4 * u5 + _f7;
        ryc = _f2 * u4 + _f5 * u5 + _f8;
        rwc = _f3 * u4 + _f6 * u5 + _f9;
        r =(u1 * rxc + u2 * ryc + rwc);
        rx = _f1 * u1 + _f2 * u2 + _f3;
        ry = _f4 * u1 + _f5 * u2 + _f6;
        p[i] = r*r / (rxc*rxc + ryc*ryc + rx*rx + ry*ry); //original, Sampson`s error

    }
}
void exFDs (const double *u, const double *F, double *p, double *w, int len)
{
    double rx, ry, rwc, ryc, rxc, r,a,b;
    int i;

    for (i=1; i<=len; i++)
    {
        rxc = _f1 * u4 + _f4 * u5 + _f7;
        ryc = _f2 * u4 + _f5 * u5 + _f8;
        rwc = _f3 * u4 + _f6 * u5 + _f9;
        r =(u1 * rxc + u2 * ryc + rwc);
        rx = _f1 * u1 + _f2 * u2 + _f3;
        ry = _f4 * u1 + _f5 * u2 + _f6;

        *w = rxc*rxc + ryc*ryc + rx*rx + ry*ry;
        *p = r*r / *w;
        *w = 1 / sqrt(*w);//original, Sampson`s error

        p ++;
        w ++;
        u += 6;
    }
}
void FDsSym (const double *u, const double *F, double *p, int len)
{
    double rx, ry, rwc, ryc, rxc, r;
    double a,b;
    int i;

    for (i=1; i<=len; i++)
    {
        rxc = _f1 * u4 + _f4 * u5 + _f7;
        ryc = _f2 * u4 + _f5 * u5 + _f8;
        rwc = _f3 * u4 + _f6 * u5 + _f9;
        r =(u1 * rxc + u2 * ryc + rwc);
        rx = _f1 * u1 + _f2 * u2 + _f3;
        ry = _f4 * u1 + _f5 * u2 + _f6;
        a =  rxc*rxc + ryc*ryc;
        b = rx*rx + ry*ry;
        *p = r*r* (a+b)/(a*b); //Mishkin.  Symmetric epipolar distance

        p ++;
        u += 6;
    }
}

void FDsSymidx (const double *mu, const double *F, double *p, int len,  int *idx, int siz)
{
    double rx, ry, rwc, ryc, rxc, r;
    const double *u;
    int i, mi;
    double a,b;

    for (mi=0; mi<siz; mi++)
    {
        i = idx[mi];
        u = mu + 6*i;



        rxc = _f1 * u4 + _f4 * u5 + _f7;
        ryc = _f2 * u4 + _f5 * u5 + _f8;
        rwc = _f3 * u4 + _f6 * u5 + _f9;
        r =(u1 * rxc + u2 * ryc + rwc);
        rx = _f1 * u1 + _f2 * u2 + _f3;
        ry = _f4 * u1 + _f5 * u2 + _f6;
        a =  rxc*rxc + ryc*ryc;
        b = rx*rx + ry*ry;
        p[i] = r*r* (a+b)/(a*b); //Mishkin.  Symmetric epipolar distance


    }
}


void FDsfull (const double *u, const double *F, double *p, int len)
{
    double rx, ry, rwc, ryc, rxc, r;
    double a,b;
    int i;

    for (i=1; i<=len; i++)
    {
        rxc = _f1 * u4 + _f4 * u5 + _f7;
        ryc = _f2 * u4 + _f5 * u5 + _f8;
        rwc = _f3 * u4 + _f6 * u5 + _f9;
        r =(u1 * rxc + u2 * ryc + rwc);
        rx = _f1 * u1 + _f2 * u2 + _f3;
        ry = _f4 * u1 + _f5 * u2 + _f6;
#ifdef SYMMETRIC_ERROR_CHECK
        a =  rxc*rxc + ryc*ryc;
        b = rx*rx + ry*ry;
        *p = r*r* (a+b)/(a*b); //Mishkin. Slower, but more precise. Symmetric epipolar distance
#else
        *p = r*r / (rxc*rxc + ryc*ryc + rx*rx + ry*ry); //original, Sampson`s error
#endif

        p ++;
        u += 6;
    }
}



void exFDsSym (const double *u, const double *F, double *p, double *w, int len)
{
    double rx, ry, rwc, ryc, rxc, r,a,b;
    int i;

    for (i=1; i<=len; i++)
    {
        rxc = _f1 * u4 + _f4 * u5 + _f7;
        ryc = _f2 * u4 + _f5 * u5 + _f8;
        rwc = _f3 * u4 + _f6 * u5 + _f9;
        r =(u1 * rxc + u2 * ryc + rwc);
        rx = _f1 * u1 + _f2 * u2 + _f3;
        ry = _f4 * u1 + _f5 * u2 + _f6;
        a =  rxc*rxc + ryc*ryc;
        b = rx*rx + ry*ry;
        *w = (a*b)/(a+b);
        *p = r*r / *w; //Mishkin. Slower, but more precise. Symmetric epipolar distance

        p ++;
        w ++;
        u += 6;
    }
}
int rroots3 (double *po, double *r)
{
    /* real roots of the polynomial of degree 3 */

    double b,c, b2, bt, v,  e;
    double p, q, D, A, cosphi, phit, R, _2R;
    b = *(po + 1) / rr_a;
    c = *(po + 2) / rr_a;
    b2 = b*b;
    bt = b/3;

    p = (3*c - b2)/ 9;
    q = ((2 * b2 * b)/27 - b*c/3 + rr_d/rr_a) / 2;

    D = q*q + p*p*p;

    if (D > 0)
    {
        A = sqrt(D) - q;
        if (A > 0)
        {
            v = pow(A,1.0/3);
            *r = v - p/v - bt;
        } else
        {
            v = pow(-A,1.0/3);
            *r = p/v - v - bt;
        }

        return 1;
    } else
    {
        if (q > 0) e = 1; else e = -1;
        R = e * sqrt(-p);
        _2R = R *2;
        cosphi = q / (R*R*R);
        if (cosphi > 1) cosphi = 1; else
            if (cosphi < -1) cosphi = -1;
        phit = acos(cosphi) /3;
       // pit = 3.14159265358979/3;

        r[0] = -_2R * cos(phit) -bt;
        r[1] =  _2R * cos(pit - phit) -bt;
        r[2] =  _2R * cos(pit + phit) -bt;

        return 3;
    }
}

void lin_fmN(const double *u, double *p, const int *inl, int len,
             double *A1, double *A2)
{
    /* linearizes corresp. with respect to entries of fundamental matrix,
      so that u' F u -> A f */

    const double  *s;
    double a[3], b[3];
    int i,k,l;

    a[2] = 1; b[2] = 1;

    s = u;
    for (i = 0; i < len; i++)
    {
        s = u + 6*inl[i];

        a[0] = *(s) * A1[0] + A1[1];
        a[1] = *(s+1) * A1[0] + A1[2];
        b[0] = *(s+3) * A2[0] + A2[1];
        b[1] = *(s+4) * A2[0] + A2[2];
        for (k = 0; k < 3; k++)
            for (l = 0; l < 3; l++)
            {
                *p = a[l] * b[k];
                p++;
            }
    }
}

void singulF(double *F)
{
    double S[3], UT[9], V[9], D[9] = {1,0,0,0,1,0,0,0,1}, VD[9];
    trnm(F,3);

    if( lap_SVD (S, F, UT, 3, V, 3) != 0 ) {
        memcpy(F, D, 9*sizeof(double));
        return;
    }

    D[0] = S[0];
    D[4] = S[1];
    D[8] = 0;

    mmul(VD,V,D,3); /*F = U.D.V^T = (V.D^T.U^T)^T = (V.D.U^T)^T*/
    mmul(F,VD,UT,3);
    trnm(F,3);
}


void u2f(const double *u, const int *inl, int len,
         double *F, double *buffer)
{
    double A1[3], A2[3];
    double *Z, V[9*9], U[8*8], D[9], *p;
    int i, j;

    if (buffer == NULL)
        Z = (double *) malloc(sizeof(double) * 9 * len);
    else
        Z = buffer;

    if (len > 8)
    {
        normu (u, inl, len, A1, A2);
        lin_fmN(u, Z, inl, len, A1, A2);

        cov_mat(V, Z, len, 9);
        lap_eig(V,D,9);
        trnm(V,9); /* lapack stores column-wise */
    } else
    {
        lin_fm(u, Z, inl, len);
        svduv(D,Z,V,9,U,8);
    }

    if (len > 8)
    {
        j = 0;
        for (i = 1; i<9; i++)
            if (D[i] < D[j]) j = i;
        p = V + j;
    } else
        p = V + 8;

    for (i = 0; i<9; i++)
    {
        F[i] = *p;
        p += 9;
    }

    singulF(F);

    if (len > 8)
        denormF(F, A1, A2);

    if (buffer == NULL)
        free (Z);
}

void u2fw(const double *u, const int *inl, const double * w,
          int len, double *F, double *buffer)
{
    double A1[3], A2[3];
    double *Z, V[9*9], U[8*8], D[9], *p;
    int i, j;

    if (buffer == NULL)
        Z = (double *) malloc(sizeof(double) * 9 * len);
    else
        Z = buffer;

    if (len > 8)
    {
        normu (u, inl, len, A1, A2);
        lin_fmN(u, Z, inl, len, A1, A2);
        for (i=0; i<len; i++)
        {
            j = inl[i];
            scalmul(Z + 9*i, w[j], 9, 1);
        }

        cov_mat(V, Z, len, 9);
        lap_eig(V,D,9);
        trnm(V,9); /* lapack stores column-wise */
    } else
    {
        lin_fm(u, Z, inl, len);
        for (i=0; i<len; i++)
        {
            j = inl[i];
            scalmul(Z+i, w[j], 9, 9);
        }
        svduv(D,Z,V,9,U,8);
    }

    if (len > 8)
    {
        j = 0;
        for (i = 1; i<9; i++)
            if (D[i] < D[j]) j = i;
        p = V + j;
    } else
        p = V + 8;

    for (i = 0; i<9; i++)
    {
        F[i] = *p;
        p += 9;
    }

    singulF(F);

    if (len > 8)
        denormF(F, A1, A2);

    if (buffer == NULL)
        free (Z);
}

/************** oriented constraints ******************/
#define xeps 1.9984e-15

void epipole(double *ec, const double *F)
{
    int i;
    crossprod(ec,F,F+6);
    for(i =0; i<3; i++)
        if ((ec[i] > xeps) || (ec[i] < -xeps)) return;
    crossprod(ec,F+3,F+6);
}

double getorisig(double *F, double *ec, double *u)
{
    double s1, s2;

    s1 = F[0]*u[3] + F[3]*u[4] + F[6]*u[5];
    s2 = ec[1]*u[2] -ec[2]*u[1];
    return(s1 * s2);
}

int all_ori_valid(double *F, double *us, int *idx, int N)
{
    double sig, sig1, ec[3], *u;
    int i;
    epipole(ec, F);
    sig1 = getorisig(F, ec, us+6*idx[0]);
    for(i = 1; i < N; i++)
    {
        u = us+6*idx[i];
        sig = getorisig(F, ec, u);
        if (sig1 * sig < 0) return 0;
    }
    return 1;
}

/***********    oriented error    ***********/

int exFDso (const double *u, const double *F, double *p, double *w, int len,
            double th, int * inl_buff, int **inls)
{
    double rx, ry, rwc, ryc, rxc, r;
    int i, po = 0, ne = 0;
    double ec[3];
    double sx, sy, sgn;

    epipole(ec, F);

    for (i=1; i<=len; i++)
    {
        rxc = _f1 * u4 + _f4 * u5 + _f7 * u6;
        ryc = _f2 * u4 + _f5 * u5 + _f8 * u6;
        rwc = _f3 * u4 + _f6 * u5 + _f9 * u6;
        r =(u1 * rxc + u2 * ryc + u3 * rwc);
        rx = _f1 * u1 + _f2 * u2 + _f3 * u3;
        ry = _f4 * u1 + _f5 * u2 + _f6 * u3;

        *w = rxc*rxc + ryc*ryc + rx*rx + ry*ry;
        *p = r*r / *w;
        if (*p < th)
        {
            sx = ec[1]*u[2] -ec[2]*u[1];
            sy = ec[2]*u[0] -ec[0]*u[2];
            sgn = sx * rxc + sy * ryc;
            if (sgn > 0)
            {
                inl_buff[po] = i;
                po ++;
            } else
            {
                ne ++;
                inl_buff[len-ne] = i;
            }
            *w = 1 / sqrt(*w);
        }
        p ++;
        w ++;
        u += 6;
    }

    if (po >= ne)
    {
        *inls = inl_buff;
        return po;
    }
    *inls = inl_buff + len - ne;
    return ne;
}

void FDso (const double *u, const double *F, double *p, double *sgn, int len)
{
    double rx, ry, rwc, ryc, rxc, r;
    int i;
    double ec[3];
    double sx, sy;
    epipole(ec, F);

    for (i=1; i<=len; i++)
    {
        rxc = _f1 * u4 + _f4 * u5 + _f7 * u6;
        ryc = _f2 * u4 + _f5 * u5 + _f8 * u6;
        rwc = _f3 * u4 + _f6 * u5 + _f9 * u6;
        r =(u1 * rxc + u2 * ryc + u3 * rwc);
        rx = _f1 * u1 + _f2 * u2 + _f3 * u3;
        ry = _f4 * u1 + _f5 * u2 + _f6 * u3;

        *p = r*r / (rxc*rxc + ryc*ryc + rx*rx + ry*ry);

        sx = ec[1]*u[2] -ec[2]*u[1];
        sy = ec[2]*u[0] -ec[0]*u[2];
        *sgn = sx * rxc + sy * ryc;

        sgn ++;
        p ++;
        u += 6;
    }
}

/*
  Insipired by code of Frederik Schaffalitzky.

  Suppose an underdetermined linear equation Ax=0, first decompose matrix A = Q R.
  We are looking for solution of
  
  A x = 0 ...  QR x = 0

  We know that Q if orthonormal QR x = 0 <-> R x = 0

  We know that there are k = cols(A) - rank(A) solutions x_1,..,x_k, with arbitrary LI tail.

  For our special case of matrix 7x9, we choose base x_1 = (*,*,*,*,*,*,*,0,1), x_2 = (*,*,*,*,*,*,*,1,0) and
  do backsubstition of R x.

*/
int nullspace_qr7x9(const double *A, double *N)
{
    const lapack_int rows=7;
    const lapack_int cols=9;
    int i,j;
    // allocate workspaces
    // change row->column organization for Fortran
#ifndef _MSC_VER
    double T[rows*cols];
    double tau[cols];
    double work[3*cols+1];
    lapack_int p[cols];
#else
    double T[7*9];
    double tau[9];
    double work[3*9+1];
    lapack_int p[9];
#endif

    lapack_int work_size = 3*cols+1;
    lapack_int info;
    // assume underdetermined system with full possible rank...
    int null_size = cols - rows;
    lapack_int k,r,c;
    double *sol = N;
    double a;

    for (i=0; i<rows; i++)
        for (j=0; j<cols; j++)
            T[i + rows*j] = A[cols*i + j];

    // prepare permutation vector
    for (j=0; j<cols; j++) p[j] = 0;

    r = rows; c = cols;
    // call Fortran LAPACK function
#ifdef _WIN32
    dgeqp3_(&r, &c, T, &r, p, tau, work, &work_size, &info);
#endif

#ifdef __linux__
    dgeqp3_(&r, &c, T, &r, p, tau, work, &work_size, &info);
#endif
    if (info!=0)
        return -1;

    // correct permutation offset
    for (j=0; j<cols; j++)
        p[j]--;

    // do backsubstitution, resulting T is column organized rows x cols
    // matrix, only elements on and above diagonal are valid and permuted
    // with permutation in p
    for (k=1;k<=null_size;k++)
    {
        // setup arbitrary part of solution vector
        for (c=rows;c<cols; c++) sol[p[c]]=0;
        sol[p[cols-k]]=1;

        // do backsubstitution
        for (r=rows-1; r>=0; r--)
        {
            a=0;
            if (T[r*rows+r]==0.0)
                return -1;
            for (c=r+1;c<cols;c++)
                a += T[c*rows+r]*sol[p[c]];
            // newvalue = -a/diagonal element
            sol[p[r]]=-a/T[r*rows+r];
        }
        sol+=cols;
    }
    return 0;
}

