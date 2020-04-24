#ifndef __FTOOLS_H__
#define __FTOOLS_H__

#define u1 (*(u))
#define u2 (*(u+1))
#define u3 (*(u+2))
#define u4 (*(u+3))
#define u5 (*(u+4))
#define u6 (*(u+5))

#define a11 (*A)
#define a12 (*(A+1))
#define a13 (*(A+2))
#define a21 (*(A+3))
#define a22 (*(A+4))
#define a23 (*(A+5))
#define a31 (*(A+6))
#define a32 (*(A+7))
#define a33 (*(A+8))

#define b11 (*B)
#define b12 (*(B+1))
#define b13 (*(B+2))
#define b21 (*(B+3))
#define b22 (*(B+4))
#define b23 (*(B+5))
#define b31 (*(B+6))
#define b32 (*(B+7))
#define b33 (*(B+8))

#define rr_a (*po)
#define rr_d (*(po + 3))

#include "Fcustomdef.h"
void lin_fm(const double *u, double *p, const int* inl, const int len);

void slcm(double *A, double *B, double *p);

#ifdef __cplusplus
extern "C"
#endif
void FDs (const double *u, const double *F, double *p, int len);


#ifdef __cplusplus
extern "C"
#endif
void FDsSym (const double *u, const double *F, double *p, int len);

#ifdef __cplusplus
extern "C"
#endif
void FDsidx (const double *mu, const double *F, double *p, int len,  int *idx, int siz);


#ifdef __cplusplus
extern "C"
#endif
void FDsSymidx (const double *mu, const double *F, double *p, int len,  int *idx, int siz);



#ifdef __cplusplus
extern "C"
#endif
void FDsfull (const double *u, const double *F, double *p, int len);


#ifdef __cplusplus
extern "C"
#endif
void exFDs (const double *u, const double *F, double *p, double *w, int len);


#ifdef __cplusplus
extern "C"
#endif
void exFDsSym (const double *u, const double *F, double *p, double *w, int len);

int rroots3 (double *po, double *r);

void lin_fmN(const double *u, double *p, const int *inl, int len,
             double *A1, double *A2);

void singulF(double *F);

void u2f(const double *u, const int *inl, int len,
         double *F, double *buffer);

void u2fw(const double *u, const int *inl, const double * w,
          int len, double *F, double *buffer);

void epipole(double *ec, const double *F);

int all_ori_valid(double *F, double *us, int *idx, int N);

int exFDso (const double *u, const double *F, double *p, double *w, int len,
            double th, int * inl_buff, int **inls);

void FDso (const double *u, const double *F, double *p, double *sgn, int len);

int nullspace_qr7x9(const double *A, double *N);

#endif /* __FTOOLS_H__ */
