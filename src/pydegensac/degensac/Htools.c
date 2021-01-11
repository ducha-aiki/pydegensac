#include <stdlib.h>

#include <memory.h>

#include <stdio.h>

#include "lapwrap.h"
#include "../matutls/matutl.h"
#include "utools.h"
#include "rtools.h"

#include "Htools.h"
#include "math.h"

#ifndef MAX
#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#endif

//#define USE_MAX_ERROR
void lin_hg(const double *u, double *dst, const int* inl, int len) /* output matrix stored column-wise */
{
    /* linearizes corresp. with respect to entries of homography matrix,
      so that u' = H u -> A h */

    const double  *s;
    double *p;
    int i,j,len6 = 6*len;

    for (i = 0; i < len; i++)
    {
        s = u + 6*inl[i];
        p = dst + 2*i;

        for (j = 3; j < 6; j++, p+=len6)
            *p = s[j];

        p -= 16*len;
        for (j =0; j<3; j++, p+=len6)
            *p = 0;

        p -= 16*len;
        for (j =3; j<6; j++, p+=len6)
            *p = -s[0] * s[j];

        p = dst + 2*i + 1;
        for (j =0; j<3; j++, p+=len6)
            *p = 0;

        p -= 16*len;
        for (j = 3; j < 6; j++, p+=len6)
            *p = s[j];

        p -= 16*len;
        for (j =3; j<6; j++, p+=len6)
            *p = -s[1] * s[j];

    }
}

void lin_hgN(const double *u, double *p, const int* inl, int len,
             double *A1, double *A2) /* output matrix stored row-wise */
{
    /* linearizes corresp. with respect to entries of homography matrix,
      so that u' = H u -> A h */

    const double  *s;
    double a[3], b[3];
    int i,j;

    a[2] = 1; b[2] = 1;

    for (i = 0; i < len; i++)
    {
        s = u + 6*inl[i];
        a[0] = *(s) * A1[0] + A1[1];
        a[1] = *(s+1) * A1[0] + A1[2];
        b[0] = *(s+3) * A2[0] + A2[1];
        b[1] = *(s+4) * A2[0] + A2[2];

        for (j = 0; j < 3; j++, p+=3)
            *p = b[j];
        p -= 8;
        for (j =0; j<3; j++, p+=3)
            *p = 0;
        p -= 8;
        for (j =0; j<3; j++, p+=3)
            *p = -a[0] * b[j];
        p -= 2;
        for (j =0; j<3; j++, p+=3)
            *p = 0;
        p -= 8;
        for (j = 0; j < 3; j++, p+=3)
            *p = b[j];
        p -= 8;
        for (j =0; j<3; j++, p+=3)
            *p = -a[1] * b[j];
        p -= 2;
    }
}

void u2h(const double *u, const int *inl, int len, double *H, double *buffer) {
    double A1[3], A2[3];
    double *Z, V[9*9], D[9];
    int i, nullspace_buff[2*9];

    if (len < 4) { /* Nothing to do */
        return;
    } else if (len == 4) { /* Nullspace */
        double Z2[9*9];
        lin_hg(u, Z2, inl, len); /* Stored column-wise... */
        trnm(Z2, 9); /* ...but we (nullspace) need it row-wise. */
        for (i = 9*8; i < 9*9; ++i) { /* Fill with zeros to square */
            Z2[i] = 0.0;
        }
        nullspace(Z2, V, 9, nullspace_buff);
        memcpy(H, V, 3*3 * sizeof(double));
    } else { /* Least Squares */
        if (!buffer) {
            Z = (double *) malloc(sizeof(double) * 9 * len * 2);
        } else {
            Z = buffer;
        }
        normu (u, inl, len, A1, A2);
        lin_hgN(u, Z, inl, len, A1, A2);
        cov_mat(V, Z, 2*len, 9);
        lap_eig(V, D, 9);
        memcpy(H, V, 3*3 * sizeof(double));
        denormH(H, A1, A2);
        if (!buffer) {
            free(Z);
        }
    }
}

void pinvJ (double a, double b, double c, double d, double e, double *pJ)
{
    double a2=a*a, b2=b*b, c2=c*c, d2=d*d, e2=e*e;
    double c2pd2 = c2+d2, ab = a*b, de = d*e;
    double Q = c * (c2pd2 + e2);
    double N;
    int i;

    pJ[0] = -b * de + a * (c2 + e2);
    pJ[1] = b * c2pd2 - a * de;
    pJ[2] = Q;
    pJ[3] = -c * (a*d + b*e);

    pJ[4] = d * (b2 + c2) - ab * e;
    pJ[5] = -ab * d + e * (a2 + c2);
    pJ[6] = pJ[3];
    pJ[7] = c * (a2 + b2 + c2);

    /*  N = b2 * c2pd2 + a2 * (c2 + e2) + c * Q +
       a2 * (c2 + e2) - 2 * ab * de; */
    N = a * pJ[0] + b * pJ[1] + c * pJ[2];

    for (i=0; i < 8; i++)
        pJ[i] /= N;
}

void HDs(const double *lin, const double * u, 
         const double *H, double *p, int len)
{
    int i, j, shift = 2*len;
    const double *l;
    double pJ[8];
    double r1, r2, a, b, c, d, e;

    for (i=0; i<len; i++)
    {
        r1 = 0;
        r2 = 0;
        l = lin + 2*i;
        for (j = 0; j < 9; j++)
        {
            r1 += H[j] * *l;
            r2 += H[j] * l[1];
            l += shift;
        }

        a = H[0] - H[2] * u[0];
        b = H[3] - H[5] * u[0];
        c = -H[8] - H[2] * u[3] - H[5] * u[4];
        d = H[1] - H[2] * u[1];
        e = H[4] - H[5] * u[1];

        pinvJ(a,b,c,d,e,pJ);


        *p = 0;
        for (j = 0; j < 4; j++)
        {
            a = pJ[j] * r1 + pJ[j+4] * r2;
            *p += a * a;
        }
        p++;
        u += 6;
    }
}


void HDsSymSumSq(const double *lin, const double * u,
            const double *H, double *p, int len) //Mishkin
{
    int i;
    double Hinv[9],H1[9];

    double den,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (i=0; i<len; i++)
    {
        den = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/den;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/den;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= xdiff*xdiff+ydiff*ydiff;

        den = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;
        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/den;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/den;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= xdiff*xdiff+ydiff*ydiff;
        *p = d1+d2;
        p++;
        u += 6;
    }
}

void HDsSymSum(const double *lin, const double * u,
            const double *H, double *p, int len) //Mishkin
{
    int i;
    double Hinv[9],H1[9];

    double den,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (i=0; i<len; i++)
    {
        den = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/den;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/den;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= sqrt(xdiff*xdiff+ydiff*ydiff);

        den = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;
        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/den;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/den;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= sqrt(xdiff*xdiff+ydiff*ydiff);
        *p = d1+d2;
        p++;
        u += 6;
    }
}
void HDsSymMax(const double *lin, const double * u,
               const double *H, double *p, int len) //Mishkin
{
    int i;
    double Hinv[9],H1[9];

    double a,b,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (i=0; i<len; i++)
    {
        a = H1[6]*u[0]+H1[7]*u[1]+H1[8];
        b = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8];

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/a;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/a;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= xdiff*xdiff+ydiff*ydiff;

        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/b;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/b;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= xdiff*xdiff+ydiff*ydiff;
        *p = sqrt(MAX(d1,d2));
        p++;
        u += 6;
    }
}
void HDsSymMaxSq(const double *lin, const double * u,
               const double *H, double *p, int len) //Mishkin
{
    int i;
    double Hinv[9],H1[9];

    double a,b,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (i=0; i<len; i++)
    {
        a = H1[6]*u[0]+H1[7]*u[1]+H1[8];
        b = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8];

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/a;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/a;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= xdiff*xdiff+ydiff*ydiff;

        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/b;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/b;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= xdiff*xdiff+ydiff*ydiff;
        *p = MAX(d1,d2);
        p++;
        u += 6;
    }
}
/* Sampson error for homography and point correspondences, computed only on a subset pts */
void HDsi(const double *lin, const double * u6,
          const double *H, double *p, int len, int *pts, int ni)
{
    int i, j, shift = 2*len;
    const double *l;
    double pJ[8];
    double r1, r2, a, b, c, d, e;
    const double *u;

    for (i=0; i<ni; i++)
    {
        u = u6 + 6*pts[i];
        r1 = 0;
        r2 = 0;
        l = lin + 2*pts[i];
        for (j = 0; j < 9; j++)
        {
            r1 += H[j] * *l;
            r2 += H[j] * l[1];
            l += shift;
        }

        a = H[0] - H[2] * u[0];
        b = H[3] - H[5] * u[0];
        c = -H[8] - H[2] * u[3] - H[5] * u[4];
        d = H[1] - H[2] * u[1];
        e = H[4] - H[5] * u[1];

        pinvJ(a,b,c,d,e,pJ);


        *p = 0;
        for (j = 0; j < 4; j++)
        {
            a = pJ[j] * r1 + pJ[j+4] * r2;
            *p += a * a;
        }
        p++;
    }
}

void invertH( const double *H, double *Hinv){
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];
    minv(Hinv,3);
}



void HDsiSymSumSq(const double *lin, const double * u6,
                  const double *H, double *p, int len, int *pts, int ni)
{
    int i;
    double Hinv[9],H1[9];
    const double *u;

    double den,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (i=0; i<ni; i++)
    {
        u = u6 + 6*pts[i];
        den = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/den;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/den;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= xdiff*xdiff+ydiff*ydiff;

        den = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;
        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/den;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/den;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= xdiff*xdiff+ydiff*ydiff;
        *p = d1+d2;
        p++;
        u += 6;
    }
}
void HDsiSymSum(const double *lin, const double * u6,
                const double *H, double *p, int len, int *pts, int ni)
{
    int i;
    double Hinv[9],H1[9];
    const double *u;

    double den,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (i=0; i<ni; i++)
    {
        u = u6 + 6*pts[i];
        den = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/den;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/den;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= sqrt(xdiff*xdiff+ydiff*ydiff);

        den = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;
        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/den;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/den;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= sqrt(xdiff*xdiff+ydiff*ydiff);
        *p = d1+d2;
        p++;
        u += 6;
    }
}
void HDsiSymMaxSq(const double *lin, const double * u6,
                const double *H, double *p, int len, int *pts, int ni)
{
    int i;
    double Hinv[9],H1[9];
    const double *u;

    double a,b,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (i=0; i<ni; i++)
    {
        u = u6 + 6*pts[i];
        a = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;
        b = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/a;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/a;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= xdiff*xdiff+ydiff*ydiff;

        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/b;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/b;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= xdiff*xdiff+ydiff*ydiff;
        *p = MAX(d1,d2);

        p++;
        u += 6;
    }
}
void HDsiSymMax(const double *lin, const double * u6,
                const double *H, double *p, int len, int *pts, int ni)
{
    int i;
    double Hinv[9],H1[9];
    const double *u;

    double a,b,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (i=0; i<ni; i++)
    {
        u = u6 + 6*pts[i];
        a = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;
        b = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/a;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/a;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= xdiff*xdiff+ydiff*ydiff;

        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/b;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/b;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= xdiff*xdiff+ydiff*ydiff;
        *p = sqrt(MAX(d1,d2));

        p++;
        u += 6;
    }
}

void HDsidx(const double *lin, const double * mu, const double *H,
            double *p, int len, int *idx, int siz)
{
    int mi, i, j, shift = 2*len;
    const double *l;
    double pJ[8];
    const double *u;
    double r1, r2, a, b, c, d, e;

    for (mi=0; mi<siz; mi++)
    {
        i = idx[mi];
        r1 = 0;
        r2 = 0;
        l = lin + 2*i;
        u = mu + 6*i;

        for (j = 0; j < 9; j++)
        {
            r1 += H[j] * *l;
            r2 += H[j] * l[1];
            l += shift;
        }

        a = H[0] - H[2] * u[0];
        b = H[3] - H[5] * u[0];
        c = -H[8] - H[2] * u[3] - H[5] * u[4];
        d = H[1] - H[2] * u[1];
        e = H[4] - H[5] * u[1];

        pinvJ(a,b,c,d,e,pJ);


        p[i] = 0;
        for (j = 0; j < 4; j++)
        {
            a = pJ[j] * r1 + pJ[j+4] * r2;
            p[i] += a * a;
        }
    }
}

void HDsSymSumSqidx(const double *lin, const double * mu,
               const double *H, double *p, int len, int *idx, int siz) //Mishkin
{
    int i,mi;
    double Hinv[9],H1[9];
    const double *u;
    double den,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (mi=0; mi<siz; mi++)
    {
        i = idx[mi];
        u = mu + 6*i;
        den = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/den;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/den;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1 = xdiff*xdiff+ydiff*ydiff;

        den = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;
        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/den;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/den;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2 = xdiff*xdiff+ydiff*ydiff;
        p[i] = d1+d2;
    }
}

void HDsSymSumidx(const double *lin, const double * mu,
               const double *H, double *p, int len, int *idx, int siz) //Mishkin
{
    int i,mi;
    double Hinv[9],H1[9];
    const double *u;
    double den,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (mi=0; mi<siz; mi++)
    {
        i = idx[mi];
        u = mu + 6*i;
        den = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/den;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/den;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1 = sqrt(xdiff*xdiff+ydiff*ydiff);

        den = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;
        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/den;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/den;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2 = sqrt(xdiff*xdiff+ydiff*ydiff);
        p[i] = d1+d2;
    }
}
void HDsSymMaxidx(const double *lin, const double * mu,
                  const double *H, double *p, int len, int *idx, int siz) //Mishkin
{
    int i,mi;
    double Hinv[9],H1[9];
    const double *u;
    double a,b,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (mi=0; mi<siz; mi++)
    {
        i = idx[mi];
        u = mu + 6*i;
        a = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;
        b = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/a;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/a;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= xdiff*xdiff+ydiff*ydiff;
        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/b;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/b;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= xdiff*xdiff+ydiff*ydiff;
        p[i] = sqrt(MAX(d1,d2));
    }
}

void HDsSymMaxSqidx(const double *lin, const double * mu,
                  const double *H, double *p, int len, int *idx, int siz) //Mishkin
{
    int i,mi;
    double Hinv[9],H1[9];
    const double *u;
    double a,b,xa,ya,d1,d2,xdiff,ydiff;
    Hinv[0] = H[0];
    Hinv[1] = H[3];
    Hinv[2] = H[6];
    Hinv[3] = H[1];
    Hinv[4] = H[4];
    Hinv[5] = H[7];
    Hinv[6] = H[2];
    Hinv[7] = H[5];
    Hinv[8] = H[8];

    for (i=0; i<9; i++)
        H1[i] = Hinv[i];
    minv(H1,3);

    for (mi=0; mi<siz; mi++)
    {
        i = idx[mi];
        u = mu + 6*i;
        a = H1[6]*u[0]+H1[7]*u[1]+H1[8] + 1e-10;
        b = Hinv[6]*u[3]+Hinv[7]*u[4]+Hinv[8] + 1e-10;

        xa = (H1[0]*u[0]+H1[1]*u[1]+H1[2])/a;
        ya = (H1[3]*u[0]+H1[4]*u[1]+H1[5])/a;
        xdiff = u[3]-xa;
        ydiff = u[4]-ya;
        d1= xdiff*xdiff+ydiff*ydiff;
        xa = (Hinv[0]*u[3]+Hinv[1]*u[4]+Hinv[2])/b;
        ya = (Hinv[3]*u[3]+Hinv[4]*u[4]+Hinv[5])/b;
        xdiff = u[0]-xa;
        ydiff = u[1]-ya;
        d2= xdiff*xdiff+ydiff*ydiff;
        p[i] = MAX(d1,d2);
    }
}


/* orientation */

int all_Hori_valid (double * us, int *idx)
{
    double p[3], q[3];
    double *a, *b, *c, *d;

    a = us + 6*idx[0];
    b = us + 6*idx[1];
    c = us + 6*idx[2];
    d = us + 6*idx[3];

    crossprod(p,a,b);
    crossprod(q,a+3,b+3);

    if ((p[0]*c[0]+p[1]*c[1]+p[2]*c[2])*(q[0]*c[3]+q[1]*c[4]+q[2]*c[5])<0)
        return 0;
    if ((p[0]*d[0]+p[1]*d[1]+p[2]*d[2])*(q[0]*d[3]+q[1]*d[4]+q[2]*d[5])<0)
        return 0;

    crossprod(p,c,d);
    crossprod(q,c+3,d+3);

    if ((p[0]*a[0]+p[1]*a[1]+p[2]*a[2])*(q[0]*a[3]+q[1]*a[4]+q[2]*a[5])<0)
        return 0;
    if ((p[0]*b[0]+p[1]*b[1]+p[2]*b[2])*(q[0]*b[3]+q[1]*b[4]+q[2]*b[5])<0)
        return 0;

    return 1;
}

//https://stackoverflow.com/a/3813723/1983544
inline int collinear(double x1, double y1, double x2, double y2, double x3, double y3) {
    const double area2 = fabs(x1 * (y2 - y3) + x2*(y3 - y1) - x3 *(y1 - y2));
    printf("%f\n", area2);
    return  area2<= 10.0;
}

int all_not_collinear(double * us, int *idx)
{
    double *a, *b, *c, *d;

    a = us + 6*idx[0];
    b = us + 6*idx[1];
    c = us + 6*idx[2];
    d = us + 6*idx[3];
    // 1st image
    if (collinear(a[0], a[1],
                  b[0], b[1],
                  c[0], c[1]))
        return 0;
    if (collinear(a[0], a[1],
                  b[0], b[1],
                  d[0], d[1]))
        return 0;
    if (collinear(b[0], b[1],
                  c[0], c[1],
                  d[0], d[1]))
        return 0;
    if (collinear(a[3], a[4],
                  c[3], c[4],
                  d[3], d[4]))
        return 0;

    // 2nd image
    if (collinear(a[3], a[4],
                  b[3], b[4],
                  c[3], c[4]))
        return 0;
    if (collinear(a[3], a[4],
                  b[3], b[4],
                  d[3], d[4]))
        return 0;
    if (collinear(b[3], b[4],
                  c[3], c[4],
                  d[3], d[4]))
        return 0;
    if (collinear(a[3], a[4],
                  c[3], c[4],
                  d[3], d[4]))
        return 0;
    return 1;
}
