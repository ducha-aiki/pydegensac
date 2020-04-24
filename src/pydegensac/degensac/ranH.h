#ifndef __RAN_H_H__
#define __RAN_H_H__

#include "rtools.h"


Score iterH (double *u, int len, int *inliers, double th, double ths, double *H, double *Z,
             double **errs, double *buffer, unsigned inlLimit);

Score inHrani (double *u, int len, int *inliers, int ninl, double th, double *Z,
               double **errs, double *buffer, double *H, unsigned inlLimit);

/* LO-RANSAC */
#ifdef __cplusplus
extern "C"
#endif
Score ransacH (double *u, int len, double th, double conf, int max_sam,
               double *H, unsigned char * inl, int *data_out, int do_lo, int inlLimit);

/* LO-RANSAC as simply as possible */
void ransacHsimple (double *u, int len, double th, double *H);

#endif /* __RAN_H_H__ */

