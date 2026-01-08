#ifndef _EXP_RANSAC_RAN_H
#define _EXP_RANSAC_RAN_H

#include "exp_ranF.h"
#include "Htools.h"

#define D3_H_RATIO 1
#define D3_H_MIN 0
/*empirically 7*mss...*/
#define D3_H_MAX 28

int HcloseToSingular(const double *h);


Score exp_iterH(double *u, int len, int *inliers, double th, double ths,
                int steps, double *H, double *Z, double **errs, double *buffer,
                int iterID, unsigned inlLimit, double *resids);

Score exp_inHrani (double *u, int len, int *inliers, int ninl,
                   double th, double *Z, double **errs,
                   double *buffer, double *H, int rep,
                   int * iterID, unsigned inlLimit, double *resids);

#ifdef __cplusplus
extern "C"
#endif
Score exp_ransacHcustomLAF (double *u, double *u_1, double *u_2, int len, double th,
                            double laf_coef,
                            double conf, int max_sam,
                            double *H, unsigned char * inl,
                            int iter_type, int * data_out,
                            int oriented_constraint, unsigned inlLimit, double **resids,
                            HDsPtr HDS1, HDsiPtr HDSi1, HDsidxPtr HDSidx1, double SymCheck_th);
/* Model-Change Error of homography (squared),
        deefined as mean square error of sample points when added the point as hard constraint */
void hMCEs(double *Z, double *u, double *d, int *samidx, int len, double * errs, double thr);

#endif

