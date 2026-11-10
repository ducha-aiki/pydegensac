#ifndef __RAN_H_H__
#define __RAN_H_H__

#include "rtools.h"


Score iterH (double *u, int len, int *inliers, double th, double ths, double *H, double *Z,
             double **errs, double *buffer, unsigned inlLimit);

Score inHrani (double *u, int len, int *inliers, int ninl, double th, double *Z,
               double **errs, double *buffer, double *H, unsigned inlLimit);

/* intbuff contains len integers and is private to this estimator call. */
Score inHrani_workspace(double *u, int len, int *inliers, int ninl, double th, double *Z,
                        double **errs, double *buffer, double *H, unsigned inlLimit,
                        int *intbuff);

#endif /* __RAN_H_H__ */

