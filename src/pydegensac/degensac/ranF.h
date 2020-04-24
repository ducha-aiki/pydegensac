#ifndef __RAN_F_H__
#define __RAN_F_H__

#include "rtools.h"


Score iterF (double *u, int len, int *inliers, double th, double ths, double *F,
             double **errs, double *buffer, unsigned inlLimit);

Score inFrani (double *u, int len, int *inliers, int ninl, double th,
               double **errs, double *buffer, double *F, unsigned inlLimit);

Score ransacF (double *u, int len, double th, double conf, int max_sam,
               double *F, unsigned char * inl, int * data_out, int do_lo, int inlLimit);

/* LO-RANSAC as simply as possible */
void ransacFsimple (double *u, int len, double th, double *F);

int prosacF(double *u, int len, double th, double conf, int * gf,
            double *F, unsigned char * inl, int * data_out, double * outn);

#endif /* __RAN_F_H__ */

