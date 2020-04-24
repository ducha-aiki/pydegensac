#ifndef _RANH2EL_H_
#define _RANH2EL_H_

/* X'Y'A'B'C'XYABC, first is u', second u! (Chum's notation :-) */
#define _u1 (u[0])
#define _u2 (u[1])
#define _u4 (u[5])
#define _u5 (u[6])
/* transposed! */
#define _a1 (A[0])
#define _a2 (A[3])
#define _a3 (A[6])
#define _a4 (A[1])
#define _a5 (A[4])
#define _a6 (A[7])
#define _a7 (A[2])
#define _a8 (A[5])
#define _a9 (A[8])
#define _b1 (B[0])
#define _b2 (B[3])
#define _b3 (B[6])
#define _b4 (B[1])
#define _b5 (B[4])
#define _b6 (B[7])
#define _b7 (B[2])
#define _b8 (B[5])
#define _b9 (B[8])

#define SQRT2 (1.4142135623730950)
/* (18/7)^2 -> ratio of medians RMSE on GT from minimal samples (2el/4pt), multiplier for threshold of points to go to LO */
#define TAU (18.0*18.0/7.0/7.0)

#include "rtools.h"

Score ransacH2el (double *u10, int len, double th, double conf, int max_sam,
                  double *H, unsigned char * inl, int *data_out, int do_lo, int inlLimit);

void getTransf (double *u10, double *N, double *D);

/*% function H = A2toRH (N1, D1, N2, D2)
%
% computes homography H from two ellipse-to-ellipse correspondences
% N1, N2 are 3x3 matrices representing affine transformations normalizing
%   ellipses in the first image to unit circles
% D1, D2 are 3x3 matrices representing affine transformations de-normalizing
%   unit circles to ellipses in the second image
%
% For more details see  Chum, Matas ICPR 2012:.
% Homography Estimation from Correspondences of Local Elliptical Features*/
int A2toRH(double *N1, double *D1, double *N2, double *D2, double *u, int *samidx, double *h);
/* More elliptical correspondences */
int AntoRH(double *u, int *inls, int len, double *h);

/* Fill Z by equation coefficients given by u */
void Zu(double * Z, double * u, int len);

/* Fill Z by equation coefficients given by u, using normalization transformations T1&T2 */
void ZuN(double * Z, double * u, double * T1, double * T2, int len);

/* Fill Z by equation coefficients given by N&D */
void Znd(double * Z, double * A, double * B, int len);

/* Compute normalization transformation from 2 points */
void norm2pt(double x1, double y1, double x2, double y2, double *T, double *iT);

/* Compute normalization transformation from more points */
void norm10(double * u10, int * inls, int len, double * T, double * iT);

/* Perform inner RANSAC followed by iterative least squares.

	For small size of inner sample use model computation by ellipses, for bigger use points. */
Score inHraniEl (double * u10, double *u6, int len, int *inliers, int ninl, double th, double *Z,
                 double **errs, double *buffer, double *H, unsigned inlLimit);



#endif //_RANH2EL_H_
