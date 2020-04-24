#ifndef __DEG_UTILS_H_
#define __DEG_UTILS_H_
/*function [deg, H, inl] = checksample(F, u7, th)
%checksample tests for the degeneracy of 7pt sample
% [deg, H, inl] = checksample(F, u7, th)
% F is the fundamental matrix
% u7 6-by-7 matrix of 2 view correspondences, th - threshold
% when the seven corrs contain at least five point correspondences (inl)
% linked by homography (H) then deg is set to sum(inl)*/

int checksample(double * F, double * u7, double th, double * H);

/*function H = Hdetect (F, u3)

% Hdetect calcucates homography from fund. matrix F and 3 point corrs.
% H = Hdetect (F, u3)
% F :  3-by-3 rank 2 fundamental matrix
% u3 : 6-by-3 correspondences of image pts in homog. coordinates
% see Hartley & Zisserman: Scene planes and homographies (p.318)*/

void Hdetect(double * F, double * u7, unsigned char * IDXS, double * H);

void dHDs(double * H, double * u, unsigned len, double * Ds, int * bufferP, double * bufferZ);

void sortDs(double * Ds, double * sDs, unsigned char * idx);

void skew_sym(double * a, double * ax);

void fillu3(double * u7, unsigned char * IDXS, double * u3a, double * u3b);

void crossp(double *h,double *u,double *v);

/*function [F, inls] = rFtH(u, hinl, th, H)*/
unsigned rFtH(double * u, unsigned char * hinl, double th, double * H, unsigned len,
              double * F, int * bufferP, double * bufferZ);

unsigned dmin(unsigned a, unsigned b);

/*%calculates number of samples needed to be done
function SampleCnt = nsamples(ni, ptNum, pf, conf)*/

/*unsigned nsamples(unsigned ni, unsigned np, unsigned ss, double conf);*/

/*function [F,inl, tinl] = innerFH (uH, uO, u, th, num, sam_sizH, sam_sizO)*/
void innerFH(double * uH, unsigned lenH, double * uO, unsigned lenO,
             double * u, unsigned len, double th, unsigned repCount, unsigned sam_sizH, unsigned sam_sizO,
             double * F, unsigned char * inl);

/*Draw combined sample from two arrays of TCs*/
void dual_sample(double * uA, unsigned lenA, unsigned sA, double * uB, unsigned lenB, unsigned sB, double * usam);

/*Iterative LSQs
function [F,inl] = u2Fit(u, F, th, ths)*/
unsigned u2Fit(double * u, unsigned len, double * F, unsigned char * inl, double th, double ths, unsigned iters);

/*Copied from mex file*/
unsigned innerH(double * H, double * u, unsigned len, double th, unsigned iters, unsigned char * inl, int * pool, double * buffer);

/*Transforms inliers from list to bit-array*/
void transformInliers(int * inl, int * inl2, unsigned inlCount, unsigned len);

#endif //__DEG_UTILS__
