#ifndef _EXP_RANSAC_RAN_F
#define _EXP_RANSAC_RAN_F

#include "rtools.h"
#include "Fcustomdef.h"

/*
3D - speeding up LO
D1: Delayed   - replace first 50 rule by smarter solution, e.g. run LO after 5 (or smth) steps after sftb and only if no other sftb appeared - not implemented yet
(IDEA: use two separate heuristics, one to turn LO on after reaching nsamples(N/2) and the second after reaching nsamples(k*I), where k is some empirically/statistically estimated ratio of number of inliers before and after the LO. Use the one which happens the first. That would try enough of promising samples on difficult scenes but would end fast on easy data)
D2: Different - dont repeat LO if inliers set is very similar - implemented as if is the same -> HASHING
D3: Detached  - use only subset of inliers (randomly chosen)
*/

#define D3_F_RATIO 1
#define D3_F_MIN 0
/*empirically 7*mss...*/
#define D3_F_MAX 49

/* Turns on DegenSAC */
#define __DEGEN__
/* Turns on inlier limit */
#define __D3__
/* Turns off oriented constraints */
/*#define __OC_OFF__*/
/* Turns on hashing of already processed inlier sets */
/*#define __HASHING__*/
/* Turns on final least squares on the result of RANSAC */
/*#define __FINAL_LSQ__*/
/* Turns on least squares at the very beginning of LO */
#define __LSQ_BEFORE_LO__
/* LSq Before LO on Model-Change Error inliers */
/*#define __LSBL_MCE__*/
/* Ibase on Model-Change Error  */
/*#define __IB_MCE__*/
/* MCE by soft constraint */
/* (Least Squares with added hi-weighted equations instead of hard constraint) */
#define __MCE_SOFT__
#define MCE_SOFT_WEIGHT 100

#ifdef __linux__
/*microseconds*/
unsigned getticks(void);
#endif /*__linux__*/

Score exp_iterF(double *u, int len, int *inliers, int * inl2, double th, double ths, int iters,
                double *F, double **errs, double *buffer, int * samidx,
                int iterID, unsigned inlLimit, double *resids);


Score exp_inFrani (double *u, int len, int *inliers, int ninl,
                   double th, double **errs, double *buffer,
                   double *F, int * samidx, int * iterID, unsigned inlLimit, double *resids);
#ifdef __cplusplus
extern "C"
#endif
int exp_ransacF (double *u, int len, double th, double conf, int max_sam,
                 double *F, unsigned char * inl, int * data_out, int do_lo, unsigned inlLimit, double **resids, double* H_best, int *Ih);

Score exp_iterFcustom(double *u, int len, int *inliers, int * inl2, double th, double ths, int iters,
          double *F, double **errs, double *buffer, int * samidx, int iterID, unsigned inlLimit, double *resids, exFDsPtr EXFDS1,FDsPtr FDS1);


Score exp_inFranicustom (double *u, int len, int *inliers, int ninl,
             double th, double **errs, double *buffer,
             double *F, int * samidx, int * iterID, unsigned inlLimit, double *resids,exFDsPtr EXFDS1,FDsPtr FDS1);


#ifdef __cplusplus
extern "C"
#endif
int exp_ransacFcustomLAF(double *u, double *u_1, double *u_2,int len, double th,  double laf_coef, double conf, int max_sam,
            double *F, unsigned char * inl,
            int * data_out, int do_lo, unsigned inlLimit, double **resids, double* H_best, int* Ih,exFDsPtr EXFDS1, FDsPtr FDS1, FDsidxPtr FDS1idx,  double SymCheck_th, int enable_degen_check);

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
void exFDs (const double *u, const double *F, double *p, double *w, int len);


#ifdef __cplusplus
extern "C"
#endif
void exFDsSym (const double *u, const double *F, double *p, double *w, int len);;


#endif // _EXP_RANSAC_RAN_F

