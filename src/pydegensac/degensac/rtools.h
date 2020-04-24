#ifndef __RTOOLS_H__
#define __RTOOLS_H__

#define DEGENSAC_EPS 2.2204e-16
#define MAX_SAMPLES 1000000
#define CONFIDENCE 0.95
#define ITER_SAM 50
#define RAN_REP 10
#define ILSQ_ITERS 4
#define TC 4

#define INL_LIMIT_H 28
#define INL_LIMIT_F 49

#define RESIDS_M (2 + RAN_REP*(1+ILSQ_ITERS+1))

/* RANSAC Scoring */
typedef struct
{
    /* number of inliers, rectangular gain function */
    unsigned I;
    /* MSAC scoring, truncated quadratic gain function */
    double J;
    /* number of inliers in both ways */
    unsigned Is;
    /* number of LAF-consistent inliers in both ways */
    unsigned Ilafs;

} Score;

#define SC_R 1
#define SC_M 2
#define SC_H 3
#define __SCORE__ SC_M

/* MSAC Width Multiplier (squared) */
#if __SCORE__ == SC_M || __SCORE__ == SC_H
#define MWM (9/4)
#else
#define MWM (1)
#endif

int sample (int *pool, int max_sz, int i);

int *randsubset (int * pool, int max_sz, int siz);

void rsample (double *data, int dat_siz,
              int *pool, int size, int max_sz, double *dst);

void addcorrT (double *src, int dat_siz, int max_sz, double *dst);

void rsampleT (double *data, int dat_siz,
               int *pool, int size, int max_sz, double *dst);

void rsampleTn (double *data, int dat_siz, int *pool,
                int size, int n, int max_sz, double *dst);

void multirsample (double *data, int dat_siz, int dps,
                   int *pool, int size, int max_sz, double *dst);

void multirsampleT (double *data, int dat_siz, int dps,
                    int *pool, int size, int max_sz, double *dst);

/* Indexes of inliers with error lower than given threshold. Returns RANSAC score. */
Score inlidxs (const double * err, int len, double th, int * inl);

/*Indexes of inliers with error lower than given threshold. Returns inliers count.*/
int inlidxso (const double * err, const double * sgn, int len, double th,
              int * inl_buff, int ** inls);

/*Number of samples to ensure given confidence*/
int nsamples(int ninl, int ptNum, int samsiz, double conf);

double truncQuad(double epsilon, double thr);

/* Score comparator */
int scoreLess(const Score s1, const Score s2);

/*Extract sample given by samidx (list of indexes)*/
void loadSample(double * u, int * samidx, unsigned sample_size, unsigned data_size, double * u_out);

#endif /* __RTOOLS_H__ */

