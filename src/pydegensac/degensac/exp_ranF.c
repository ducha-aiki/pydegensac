#undef __STRICT_ANSI__
////#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "exp_ranF.h"
#include "DegUtils.h"
#include "hash.h"
#include "bsd_random.h"

static DEGENSAC_THREAD_LOCAL HashTable HASH_TABLE_F;

#include "Ftools.h"
#include "rtools.h"
#include "utools.h"
#include "bsd_random.h"
#include "../matutls/matutl.h"
#include <time.h>
//#include <mex.h>

#define CHECK_ORIG(inl,orig) 1
//#define CHECK_ORIG(inl,orig) (inl[orig[0]] && inl[orig[1]] && inl[orig[2]] && inl[orig[3]] && inl[orig[4]] && inl[orig[5]] && inl[orig[6]]) //to be used only with DegenSAC

#define CHECK_COEF 16.0
#define SYMM_COEF 0.6
#define __HASHING__

#ifdef __linux__
#include<sys/time.h>
/*microseconds*/
/*static __inline__*/ unsigned getticks(void) {
    unsigned s, us;
    struct timeval tv;
    gettimeofday(&tv, 0);
    s = tv.tv_sec;
    us = tv.tv_usec;
    return (s%10)*1000*1000 + us;
}
#endif /*__linux__*/

#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

static void reseed_rng(unsigned seed_value) {
    degensac_srand(seed_value);
    degensac_srandom(seed_value);
}

static int init_lo_workspace(int len, double **workspace, double **weights,
                             int **intbuff, int **intbuff2, int **intbuff_best)
{
    int *ints;

    *workspace = (double *) malloc((size_t) len * sizeof(double) +
                                   3 * (size_t) len * sizeof(int));
    if (*workspace == NULL)
        return 0;
    *weights = *workspace;
    ints = (int *) (*workspace + len);
    *intbuff = ints;
    *intbuff2 = ints + len;
    *intbuff_best = ints + 2 * len;
    return 1;
}

int no_mto(double *A)
{
    double x,y;
    int i,j;
    int pts_ok = 1;
    for (i=0;i<7;i++)
        for (j=i+1;j<7;j++)
            pts_ok = 0;//()
}


/******* Custom *********/
static Score exp_iterFcustom_workspace(double *u, int len, int *inliers, int * inl2, double th, double ths, int iters,
                      double *F, double **errs, double *buffer, int * samidx, int iterID, unsigned inlLimit, double *resids, exFDsPtr EXFDS1,FDsPtr FDS1,
                      const double *u_soa, FDsPtr FDS1soa, double *w) {
    double *d = errs[1];
    double f[9], dth;
    unsigned it;
    Score S = {0,0}, Ss, maxS;
#ifdef __D3__
    int * detachedInl;
    unsigned detachedCount;
#endif //__D3__
#ifdef __HASHING__
    int iterIDret;
    uint32_t hash;
#endif //__HASHING__
    dth = (ths - th) / ILSQ_ITERS;

    /* F from the sample inliers by th */
    maxS = inlidxs(errs[4], len, th, inliers);
    if (maxS.I < 8) {
        return S;
    }
    S = inlidxs(errs[4], len, th*MWM, inliers);
#ifdef __D3__
    //D3 - calculate number of inliers detached for LSQ
    detachedCount = (int)(S.I * D3_F_RATIO);
    if (detachedCount < D3_F_MIN) {
        detachedCount = D3_F_MIN;
    }
    if (detachedCount > inlLimit) {
        detachedCount = inlLimit;
    }
    if (detachedCount < 8) {
        detachedCount = 8;
    }
    if (detachedCount >= S.I) { // if we want to use more than (or all) we have, just use what we have without shuffling
        u2f(u, inliers, S.I, f, buffer);
    } else {
        detachedInl = randsubset (inliers, S.I, detachedCount);
        u2f(u, detachedInl, detachedCount, f, buffer);
    }
#else //__D3__
    u2f(u, inliers, S.I, f, buffer);
#endif //__D3__

    /*iterate */
    for (it = 0; it < iters; it ++) {
        EXFDS1 (u, f, d, w, len);
        if (resids != NULL) {
            memcpy(resids + it*len, d, len*sizeof(double));
        }

        S = inlidxs(d, len, th, inliers);

#ifdef __HASHING__
        hash = SuperFastHash((const char *)inliers, S.I * sizeof(*inliers));
        iterIDret = htContains(&HASH_TABLE_F, hash, S.I, iterID);
        if (iterIDret != -1 && iterIDret != iterID) {
            S.I = 0;
            S.J = 0;
            return S;
        }
        if (iterIDret == -1) {
            htInsert(&HASH_TABLE_F, hash, S.I, iterID);
        }
#endif //__HASHING__

        transformInliers(inliers, inl2, S.I, len);

        if (scoreLess(maxS, S) && CHECK_ORIG(inl2,samidx)) {
            maxS = S;
            errs[1] = errs[0];
            errs[0] = d;
            d = errs[1];
            memcpy(F, f, 9*sizeof(double));
        }

        Ss = inlidxs (d, len, ths*MWM, inliers);
        if (Ss.I < 8) {
            return maxS;
        }

#ifdef __D3__
        //D3 - calculate number of inliers detached for LSQ
        detachedCount = (int)(Ss.I * D3_F_RATIO);
        if (detachedCount < D3_F_MIN) {
            detachedCount = D3_F_MIN;
        }
        if (detachedCount > inlLimit) {
            detachedCount = inlLimit;
        }
        if (detachedCount < 8) {
            detachedCount = 8;
        }
        if (detachedCount >= Ss.I) { // if we want to use more than (or all) we have, just use what we have without shuffling
            u2fw(u, inliers, w, Ss.I, f, buffer);
        } else {
            detachedInl = randsubset (inliers, Ss.I, detachedCount);
            u2fw(u, detachedInl, w, detachedCount, f, buffer);
        }
#else //__D3__
        u2fw(u, inliers, w, Ss.I, f, buffer);
#endif //__D3__

        ths -= dth;
    }

    FDS1soa (u_soa, f, d, len);
    if (resids != NULL) {
        memcpy(resids + 4*len, d, len*sizeof(double));
    }
    S = inlidxs (d, len, th, inliers);
    transformInliers(inliers, inl2, S.I, len);
    if (scoreLess(maxS, S) && CHECK_ORIG(inl2,samidx))
    {
        maxS = S;
        errs[1] = errs[0];
        errs[0] = d;
        memcpy(F, f, 9*sizeof(double));
    }

    return maxS;
}

Score exp_iterFcustom(double *u, int len, int *inliers, int *inl2, double th, double ths, int iters,
                      double *F, double **errs, double *buffer, int *samidx, int iterID, unsigned inlLimit, double *resids, exFDsPtr EXFDS1, FDsPtr FDS1,
                      const double *u_soa, FDsPtr FDS1soa) {
    Score S;
    double *w = (double *) malloc((size_t) len * sizeof(double));

    if (w == NULL)
        return (Score) {0, 0};
    S = exp_iterFcustom_workspace(u, len, inliers, inl2, th, ths, iters, F, errs, buffer,
                                  samidx, iterID, inlLimit, resids, EXFDS1, FDS1, u_soa, FDS1soa, w);
    free(w);
    return S;
}

static Score exp_inFranicustom_workspace (double *u, int len, int *inliers, int ninl,
                         double th, double **errs, double *buffer,
                         double *F, int * samidx, int * iterID, unsigned inlLimit, double *resids,exFDsPtr EXFDS1,FDsPtr FDS1,
                         const double *u_soa, FDsPtr FDS1soa, int *intbuff, int *intbuff2, int *intbuff_best, double *w) {
    unsigned ssiz, i;
    Score S = {0, 0}, maxS = {0, 0};
    int jj;
    double *d, f[9];
    int *sample;
    if (ninl < 16) {
        /*//printf("Prematurely escaped LO, not enough inliers (<16)!\n");*/
        if (resids != NULL) {
            memset(resids, 0, (RESIDS_M-2)*len*sizeof(double));
        }
        return maxS; /*Zeros*/
    }
    ssiz = ninl / 2;
    if (ssiz > 14) {
        ssiz = 14;
    }

    d = errs[2];
    errs[2] = errs[0];
    errs[0] = d;

    for (i = 0; i < RAN_REP; i++) {
        sample = randsubset(inliers, ninl, ssiz);
        u2f(u, sample, ssiz, f, buffer);
        FDS1soa (u_soa, f, errs[0], len);
        if (resids != NULL) {
            memcpy(resids + i*6*len, errs[0], len*sizeof(double)); // pointer to resids already moved to the 3rd field of current part
        }
        errs[4] = errs[0];

        S = exp_iterFcustom_workspace(u, len, intbuff, intbuff2, th, TC*th, ILSQ_ITERS, f, errs, buffer, samidx, ++*iterID, inlLimit, resids != NULL ? resids + i*6*len + len : NULL,EXFDS1,FDS1,u_soa,FDS1soa,w);
        if (scoreLess(maxS, S)) {
            maxS = S;
            d = errs[2];
            errs[2] = errs[0];
            errs[0] = d;
            memcpy(F,f,9*sizeof(double));
            for ( jj =0; jj <maxS.I ; jj++) {
                intbuff_best[jj] = intbuff[jj];
            }
        }
    }

    d = errs[2];
    errs[2] = errs[0];
    errs[0] = d;
    for ( jj =0; jj <maxS.I ; jj++) {
        inliers[jj] = intbuff_best[jj];
    }
    return maxS;
}

Score exp_inFranicustom (double *u, int len, int *inliers, int ninl,
                         double th, double **errs, double *buffer,
                         double *F, int *samidx, int *iterID, unsigned inlLimit, double *resids, exFDsPtr EXFDS1, FDsPtr FDS1,
                         const double *u_soa, FDsPtr FDS1soa) {
    Score S;
    double *w;
    int *work;

    if (ninl < 16)
        return exp_inFranicustom_workspace(u, len, inliers, ninl, th, errs, buffer, F, samidx,
                                           iterID, inlLimit, resids, EXFDS1, FDS1, u_soa, FDS1soa,
                                           NULL, NULL, NULL, NULL);
    w = (double *) malloc((size_t) len * sizeof(double));
    work = (int *) malloc(3 * (size_t) len * sizeof(int));
    if (w == NULL || work == NULL) {
        free(w);
        free(work);
        return (Score) {0, 0};
    }
    S = exp_inFranicustom_workspace(u, len, inliers, ninl, th, errs, buffer, F, samidx,
                                    iterID, inlLimit, resids, EXFDS1, FDS1, u_soa, FDS1soa,
                                    work, work + len, work + 2 * len, w);
    free(work);
    free(w);
    return S;
}

static inline void sampleF7(const double *u, int *pool, int len, double *A) {
    int i, k, l;
    for (i = 0; i < 7; i++) {
        int q = sample(pool, len, i);
        const double *s = u + 6 * q;
        double *row = A + 9 * i;
        int pos = 0;
        for (k = 0; k < 3; k++)
            for (l = 0; l < 3; l++)
                row[pos++] = s[k+3] * s[l];
    }
}

int exp_ransacFcustomLAF(double *u, double *u_1, double *u_2, int len, double th, double laf_coef,
                         double conf, int max_sam,
                         double *F, unsigned char * inl,
                         int * data_out, int do_lo, unsigned inlLimit, double **resids, exFDsPtr EXFDS1, FDsPtr FDS1, FDsidxPtr FDS1idx, double SymCheck_th,
                         int enable_degen_check, int seed, FDsPtr FDS1soa) {

    int *pool, no_sam = 0, new_sam;  double *buffer, *u_soa, u7[6*7], H[3*3], FBest[3*3];
    int *bufferP;
    double *f1, *f2;
    int do_update;
    double poly[4], roots[3], f[9], *err, *d, *d_check;
    double *errs[5];
    int nsol, i, j, *inliers, new_max, do_iterate;
    unsigned I;
    Score maxS = {0,0}, maxSs = {0,0}, S = {0,0}, Scheck={0,0};
    int *samidx, samidxBest[7];
    double *errorsBest; /*To store best non-deg inls for ALO LO*/ //TODO replace with errs[3], there are the same data! BUT! do we want BEST data or BEST NON-DEGENERATED data? :(
    /* to eliminate */
    int degen_cnt = 0, iter_cnt = 0, LmaxI, iterID = 0;
    unsigned non_degen_samples_count = 0; //only those with so-far-the-best model
    double jj;
    double *HDs;
    int bad_model = 0;
    const int doSymCheck = SymCheck_th > 0;
    const int DO_LAF_CHECK =  laf_coef > 0;
    const double th_laf_check = laf_coef * th;
    int p1_inliers = 0;
    double *err_laf;
    double *workspace, *double_cursor, *lo_workspace, *lo_weights;
    int *int_cursor, *lo_intbuff, *lo_intbuff2, *lo_intbuff_best;
    size_t workspace_doubles, workspace_ints;
    int a, model_found = 0;

    lo_workspace = NULL;
    if (resids != NULL)
        *resids = NULL;

    if (seed >= 0) {
        reseed_rng((unsigned)seed);
    } else {
        reseed_rng((unsigned)time(NULL)); // Mishkin - randomization
    }

#ifdef USE_QR
    double A[7*9], sol[2*9];
#else
    double A[9*9], sol[9*9];
    int nullspace_buff[2*9];
    int nullsize;
    for (i=7*9; i<9*9; i++) {
        A[i] = 0.0;
    }
#endif

#ifdef __HASHING__
    htInit(&HASH_TABLE_F);
#endif // __HASHING__

    ////printf("__PROFILE: BEFORE ransac: %d\n", getticks()/1000);

    /* Most O(len) scratch lives in one per-call allocation. The LO block is
       allocated only if an LO with enough inliers actually runs. */
    workspace_doubles = 4 + 18 + 4;
    workspace_ints = 2; /* pool and inliers */
    if (do_lo) {
        workspace_doubles++; /* errorsBest */
    }
    if (doSymCheck)
        workspace_doubles++;
    if (DO_LAF_CHECK)
        workspace_doubles++;
    if (enable_degen_check) {
        workspace_doubles++;
        workspace_ints++;
    }
    workspace = (double *) malloc(workspace_doubles * len * sizeof(double) +
                                  workspace_ints * len * sizeof(int));
    if (workspace == NULL)
        goto allocation_failure;
    double_cursor = workspace;
    int_cursor = (int *) (workspace + workspace_doubles * len);

    u_soa = double_cursor;
    double_cursor += 4 * len;
    buffer = double_cursor;
    double_cursor += 18 * len;
    if (do_lo) {
        errorsBest = double_cursor;
        double_cursor += len;
    } else {
        errorsBest = NULL;
    }
    err = double_cursor;
    double_cursor += 4 * len;
    if (DO_LAF_CHECK) {
        err_laf = double_cursor;
        double_cursor += len;
    } else {
        err_laf = NULL;
    }
    if (doSymCheck) {
        d_check = double_cursor;
        double_cursor += len;
    } else {
        d_check = NULL;
    }
    if (enable_degen_check) {
        HDs = double_cursor;
        double_cursor += len;
    } else {
        HDs = NULL;
    }
    lo_weights = NULL;

    pool = int_cursor;
    int_cursor += len;
    if (enable_degen_check) {
        bufferP = int_cursor;
        int_cursor += len;
    } else {
        bufferP = NULL;
    }
    inliers = int_cursor;
    int_cursor += len;
    lo_intbuff = lo_intbuff2 = lo_intbuff_best = NULL;

    for (i = 0; i < len; i ++) {
        pool[i] = i;
    }
    samidx = pool + len - 7;

    /* Coordinates repacked once per call into x1|y1|x2|y2, so the whole-array
       error evaluations read contiguously instead of gathering four doubles
       out of every six. Same arithmetic, same order, same results. */
    soa_pack(u, u_soa, len);
    for (i=0; i<4; i++)
        errs[i] = err + i * len;
    errs[4] = errs[3];


    maxS.I  = 8;
    maxSs.I = 8;
    // max_sam = MAX_SAMPLES;
    no_sam = 0;

    f1 = sol;
    f2 = sol+9;

    /*  srand(RAND_SEED++); */
    while(no_sam < max_sam) {
        no_sam ++;

        sampleF7(u, pool, len, A);
        loadSample(u, samidx, 7, 6, u7);


#if USE_QR
        /* QR. Measured against the LU path below on M1 (2026-08-12): 16%
           slower end to end (288 vs 342 estimator calls per 10 s) and
           statistically indistinguishable in accuracy -- min p = 0.585 over
           40 comparisons at 300 seeds, every KS <= 0.063. LU stays the
           default; this is kept because it is the reference formulation, not
           because it is worth enabling.

           A failed factorisation means a degenerate sample; skip it, as the
           LU path does when the null space is not 2-dimensional. The return
           value used to be discarded, which fed an uninitialised `sol` into
           the model loop. */
        if (nullspace_qr7x9(A, sol) != 0) {
            continue;
        }
#else
        /* use LU */
        for (i = 7*9; i < 9*9; ++i) {
            A[i] = 0.0;
        }

        nullsize = nullspace(A, f1, 9, nullspace_buff);
        if (nullsize != 2) {
            continue;
        }
#endif
        slcm (f1, f2, poly);
        nsol = rroots3(poly, roots);

        new_max = 0; do_iterate = 0;
        LmaxI = 0;
        for (i = 0; i < nsol; i++) {
            for (j = 0; j < 9; j++) {
                f[j] = f1[j] * roots[i] + f2[j] * (1 -roots[i]);
            }

            /* orient. constr. */
#ifndef __OC_OFF__
            if (!all_ori_valid(f, u, samidx, 7))  continue;
#endif

            d = errs[i];
            FDS1soa(u_soa, f, d, len);
            S = inlidxs(d, len, th, inliers);

            if (S.I > LmaxI) LmaxI = S.I;

            if(scoreLess(maxS, S)) {
                ///
                if (doSymCheck) { //Mishkin. Check by symmetrical distance
                    FDsSymidx(u, f, d_check, len,inliers, S.I);
                    S.Is = 0;
                    for (j = 0; j < S.I; j++)
                        if (d_check[inliers[j]] <= SymCheck_th) S.Is++;


                    if (S.Is < maxS.Is)
                        continue;
                }

                if (DO_LAF_CHECK) {
                    //Mishkin. Check LAF (2 points on ellipse) symmetrical distance
                    FDS1idx(u_1, f, err_laf, len,inliers, S.I);
                    p1_inliers = 0;
                    S.Ilafs = 0;
                    for (j = 0; j < S.I; j++)
                        if (err_laf[inliers[j]] <= th_laf_check) p1_inliers++;

                    // end_of_check

                    //Mishkin. Now second point symmetrical distance
                    FDS1idx(u_2, f, err_laf, len,inliers, S.I);
                    for (j = 0; j < S.I; j++)
                        if (err_laf[inliers[j]] <= th_laf_check)  S.Ilafs++;

                    S.Ilafs = min(S.Ilafs, p1_inliers);
                    if (S.Ilafs < maxS.Ilafs)
                        continue; //skip if first is not good
                }

                errs[i] = errs[3];
                errs[3] = d;
                maxS = S;
                ////printf("I risen in main loop to %u.\n", maxS.I);
                memcpy(F,f,9*sizeof(double)); /*!!!*/
                model_found = 1;
                new_max = 1;

            }



            if(scoreLess(maxSs, S)) {
                maxSs = S;
                ////printf("__PROFILE: BEFORE checksample: %d\n", getticks()/1000);
                if (enable_degen_check && checksample(f, u7, 3*th, H)) {
                    ////printf("__PROFILE: AFTER  checksample: %d\n", getticks()/1000);
                    dHDs(H, u, len, HDs, bufferP, buffer);
                    I = 0;
                    for (j = 0; j < len; ++j) {
                        if (HDs[j] < th*3) {
                            ++I;
                        }
                    }
                    if (I < 8) {
                        break;
                    }
                    ////printf("__PROFILE: BEFORE innerH: %d\n", getticks()/1000);
                    ////printf("I before innrH %u.\n", I);

                    I = innerH(H, u, len, 16*th, 10, inl, bufferP, buffer); /*originally was 30 reps, lowered because of bad performance*/

                    ////printf("I after innrH %u.\n", I);

                    ////printf("__PROFILE: AFTER  innerH: %d\n", getticks()/1000);
                    if (I > 6) {
                        /*[aF, v] = rFtH(u, ahi, th, aH);
                                                no_i = sum(v);
                                                f//printf(1,'P+P %d %d\n',sum(ahi), no_i);*/
                        ////printf("__PROFILE: BEFORE rFtH: %d\n", getticks()/1000);
                        I = rFtH(u, inl, th, H, len, f, bufferP, buffer);
                        ////printf("I after rFtH %u.\n", I);

                        ////printf("__PROFILE: AFTER  rFtH: %d\n", getticks()/1000);
                        if(I > maxS.I) { //TODO hybrid scoring down to rFtH?
                            FDS1soa(u_soa, f, errs[3], len);
                            maxS.I = I; /*maxS.J is set later*/
                            ////printf("I risen in degen to %u.\n", maxS.I);
                            memcpy(F,f,3*3*sizeof(double));
                            model_found = 1;
                            new_max = 1;
                            d = errs[3]; /*For IJ calculation*/
                        } else {
                            FDS1soa(u_soa, f, errs[i], len);
                            d = errs[i]; /*For IJ calculation*/
                        }
                        I = 0;
                        jj = 0;
                        for (j = 0; j < len; j++) {
                            if (d[j] <= th) {
                                I++;
                            }
                            jj += truncQuad(d[j],th);
                        }
                        if (new_max) {
                            maxS.J = jj;
                        }
                        ++degen_cnt;
                    }
                } else {
                    ////printf("__PROFILE: AFTER  checksample: %d\n", getticks()/1000);
                    do_iterate = (do_lo>0 && (no_sam > ITER_SAM));
                    errs[4] = d;
                    non_degen_samples_count++;
                    memcpy(samidxBest, samidx, 7 * sizeof(int));
                    if (do_lo)
                        memcpy(errorsBest, d, len * sizeof(double));
                    memcpy(FBest, f, 3*3*sizeof(double));
                }
            }
        }

        if (data_out != NULL)
            data_out[LmaxI+2] ++;

        if (do_lo>0 && (no_sam == ITER_SAM) && non_degen_samples_count) {
            do_iterate = 1;
        }

        if (do_iterate) {
            iter_cnt ++;
            if (resids != NULL) {
                *resids = (double *) realloc(*resids, iter_cnt * RESIDS_M * len * sizeof(double));
                memcpy(*resids + RESIDS_M*(iter_cnt - 1)*len, errs[4], len*sizeof(double));
            }

#ifdef __LSQ_BEFORE_LO__
            d = errs[0];
            S = inlidxs(errs[4], len, TC*th*MWM, inliers);
            u2f(u, inliers, S.I, f, buffer);
            FDS1soa(u_soa, f, d, len);
            S = inlidxs(d, len, th, inliers);
            if (resids != NULL) {
                memcpy(*resids + RESIDS_M*(iter_cnt - 1)*len + len, d, len*sizeof(double));
            }
#else
            S = inlidxs(errs[4], len, th, inliers);
#endif /* __LSQ_BEFORE_LO__ */
            /*******/
            if (S.I >= 16 && lo_workspace == NULL &&
                !init_lo_workspace(len, &lo_workspace, &lo_weights, &lo_intbuff,
                                   &lo_intbuff2, &lo_intbuff_best))
                goto allocation_failure;
            S = exp_inFranicustom_workspace(u, len, inliers, S.I, th, errs, buffer, f, samidx, &iterID, inlLimit,
                                  resids != NULL ? *resids + 2*len + (iter_cnt-1)*RESIDS_M*len : NULL,EXFDS1,FDS1,u_soa,FDS1soa,
                                  lo_intbuff, lo_intbuff2, lo_intbuff_best, lo_weights);
            /*******/
            // minimalistic LO' (just one iterations)
            /*			S = exp_iterF(u, len, inliers, bufferP, th, 16*TC*th, 10, f, errs, buffer, samidx, ++iterID, inlLimit, *resids + 2*len + (iter_cnt-1)*RESIDS_M*len);*/
            /*******/

            if(scoreLess(maxS, S)) {

                char do_update = 1;
                if (doSymCheck) { //Mishkin. Check by symmetrical distance
                    FDsSymidx(u, f, d_check, len,inliers, S.I);
                    S.Is = 0;
                    for (j = 0; j < S.I; j++) {
                        if (d_check[inliers[j]] <= SymCheck_th) S.Is++;
                    }

                    if (S.Is < maxS.Is)
                        do_update = 0;
                }
                if (DO_LAF_CHECK && do_update) {

                    //Mishkin. Check LAF (2 points on ellipse)  distance
                    FDS1idx(u_1, f, err_laf, len,inliers, S.I);
                    p1_inliers = 0;
                    S.Ilafs = 0;
                    for (j = 0; j < S.I; j++)
                        if (err_laf[inliers[j]] <= th_laf_check) p1_inliers++;

                    // end_of_check

                    //Mishkin. Now second point  distance
                    FDS1idx(u_2, f, err_laf, len,inliers, S.I);
                    for (j = 0; j < S.I; j++)
                        if (err_laf[inliers[j]] <= th_laf_check)  S.Ilafs++;

                    S.Ilafs = min(S.Ilafs, p1_inliers);
                    if (S.Ilafs < maxS.Ilafs)
                        do_update = 0; //skip if first is not good
                }



                if (do_update) {
                    ////printf("LO %d\n", I);
                    d = errs[0];
                    errs[0] = errs[3];
                    errs[3] = d;
                    maxS = S;
                    ////printf("I risen in LO to %u.\n", maxI);
                    memcpy(F,f,9*sizeof(double)); /*!!!*/
                    model_found = 1;
                    new_max = 1;
                }
            }

            if (new_max) {
                new_sam = nsamples(maxS.I+1, len, 7, conf);
                if (new_sam < max_sam) {
                    max_sam = new_sam;
                }
            }
        }
    }
    /*If there were no LOs, do at least one NOW!*/
    if (do_lo && (!iter_cnt && !degen_cnt) && non_degen_samples_count) { //TODO maybe not a good idea to supress LO after degen...? full vs. full+
        ////printf("Running ALO LO\n");
        loadSample(u, samidxBest, 7, 6, u7);
        if (enable_degen_check && checksample(FBest, u7, 3*th, H)) { //TODO is this necessary? NO, if running full+ version (degen_cnt > 0 if found deg. sample with big consensus)
            ////printf("__PROFILE: AFTER  checksample: %d\n", getticks()/1000);
            dHDs(H, u, len, HDs, bufferP, buffer);
            I = 0;
            for (j = 0; j < len; ++j) {
                if (HDs[j] < th*3) {
                    ++I;
                }
            }
            if (I >= 8) {
                ////printf("__PROFILE: BEFORE innerH: %d\n", getticks()/1000);
                I = innerH(H, u, len, 16*th, 10, inl, bufferP, buffer); /*originally was 30 reps, lowered because of bad performance*/
                ////printf("__PROFILE: AFTER  innerH: %d\n", getticks()/1000);
            }
            if (I > 6) {
                /*[aF, v] = rFtH(u, ahi, th, aH);
                                no_i = sum(v);
                                f//printf(1,'P+P %d %d\n',sum(ahi), no_i);*/
                ////printf("__PROFILE: BEFORE rFtH: %d\n", getticks()/1000);
                I = rFtH(u, inl, th, H, len, f, bufferP, buffer);
                ////printf("__PROFILE: AFTER  rFtH: %d\n", getticks()/1000);
                if(I > maxS.I) { //TODO hybrid scoring down to rFtH?
                    FDS1soa(u_soa, f, errs[3], len);
                    maxS.I = I; /*maxS.J is set later*/
                    ////printf("I risen in degen to %u.\n", maxS.I);
                    memcpy(F,f,3*3*sizeof(double));
                    model_found = 1;
                    new_max = 1;
                    d = errs[3]; /*For IJ calculation*/
                } else {
                    FDS1soa(u_soa, f, errs[i], len);
                    d = errs[i]; /*For IJ calculation*/
                }
                I = 0;
                jj = 0;
                for (j = 0; j < len; j++) {
                    if (d[j] <= th) {
                        I++;
                    }
                    jj += truncQuad(d[j],th);
                }
                if (new_max) {
                    maxS.J = jj;
                }
                ++degen_cnt;
            }
        } else {
            iter_cnt ++;

            if (resids != NULL) {
                *resids = (double *) realloc(*resids, iter_cnt * RESIDS_M * len * sizeof(double));
                memcpy(*resids + RESIDS_M*(iter_cnt - 1)*len, errs[4], len*sizeof(double));
            }
#ifdef __LSQ_BEFORE_LO__
            d = errs[0];
            S = inlidxs(errorsBest, len, TC*th*MWM, inliers);
            u2f(u, inliers, S.I, f, buffer);
            FDS1soa(u_soa, f, d, len);
            S = inlidxs(d, len, th, inliers);
            if (resids != NULL) {
                memcpy(*resids + RESIDS_M*(iter_cnt - 1)*len + len, d, len*sizeof(double));
            }
#else
            S = inlidxs(errorsBest, len, th, inliers);
#endif /* __LSQ_BEFORE_LO__ */
            /*******/
            if (S.I >= 16 && lo_workspace == NULL &&
                !init_lo_workspace(len, &lo_workspace, &lo_weights, &lo_intbuff,
                                   &lo_intbuff2, &lo_intbuff_best))
                goto allocation_failure;
            S = exp_inFranicustom_workspace (u, len, inliers, S.I, th, errs, buffer, f, samidxBest, &iterID, inlLimit,
                                   resids != NULL ? *resids + 2*len + (iter_cnt-1)*RESIDS_M*len : NULL,EXFDS1,FDS1,u_soa,FDS1soa,
                                   lo_intbuff, lo_intbuff2, lo_intbuff_best, lo_weights);
            /*******/
            // minimalistic LO' (just one iterations)
            /*			S = exp_iterF(u, len, inliers, bufferP, th, 16*TC*th, 10, f, errs, buffer, samidxBest, ++iterID, inlLimit, *resids + 2*len + (iter_cnt-1)*RESIDS_M*len);*/
            /*******/

            if(scoreLess(maxS, S)) {
                char do_update = 1;
                if (doSymCheck) { //Mishkin. Check by symmetrical distance
                    FDsSymidx(u, f, d_check, len,inliers, S.I);
                    S.Is = 0;
                    for (j = 0; j < S.I; j++)
                        if (d_check[inliers[j]] <= SymCheck_th) S.Is++;


                    if (S.Is < maxS.Is)
                        do_update = 0;
                }
                if (DO_LAF_CHECK && do_update) {

                    //Mishkin. Check LAF (2 points on ellipse) symmetrical distance
                    FDS1idx(u_1, f, err_laf, len,inliers, S.I);
                    p1_inliers = 0;
                    S.Ilafs = 0;
                    for (j = 0; j < S.I; j++)
                        if (err_laf[inliers[j]] <= th_laf_check) p1_inliers++;

                    // end_of_check

                    //Mishkin. Now second point symmetrical distance
                    FDS1idx(u_2, f, err_laf, len,inliers, S.I);
                    for (j = 0; j < S.I; j++)
                        if (err_laf[inliers[j]] <= th_laf_check)  S.Ilafs++;

                    S.Ilafs = min(S.Ilafs, p1_inliers);
                    if (S.Ilafs < maxS.Ilafs)
                        do_update = 0; //skip if first is not good
                }

                if (do_update) {
                    ////printf("LO %d\n", I);
                    d = errs[0];
                    errs[0] = errs[3];
                    errs[3] = d;
                    maxS = S;
                    ////printf("I risen in LO to %u.\n", maxI);
                    memcpy(F,f,9*sizeof(double)); /*!!!*/
                    model_found = 1;
                    new_max = 1;
                }
            }
        }
    }

    if (!model_found) {
        memset(F, 0, 9 * sizeof(double));
        memset(inl, 0, len * sizeof(unsigned char));
        maxS.I = 0;
        goto cleanup;
    }

    d = errs[3];

#ifdef __FINAL_LSQ__
    I = inlidxs(d, len, th, inliers); //LSQ in the end
    u2f(u, inliers, I, F, buffer);
    FDS1soa(u_soa, F, d, len);
#endif

    for (j = 0; j < len; j++) {
        if (d[j] <= th) {
            inl[j] = 1;
        } else {
            inl[j] = 0;
        }

    }
    if (doSymCheck) { //Mishkin. Check by symmetrical distance

        S =  inlidxs(d, len, th, inliers);
        FDsSymidx(u, F, d_check, len, inliers, S.I);
        for (j = 0; j < S.I; j++)
            if (d_check[inliers[j]]> SymCheck_th)
                inl[j] = 0;

    }
    if (DO_LAF_CHECK && do_update) {
        S =  inlidxs(d, len, th, inliers);

        //Mishkin. Check LAF (2 points on ellipse) symmetrical distance
        FDS1idx(u_1, F, err_laf, len,inliers, S.I);
        for (j = 0; j < S.I; j++)
            if (err_laf[inliers[j]] > th_laf_check)
                inl[j] = 0;

        // end_of_check

        //Mishkin. Now second point symmetrical distance
        FDS1idx(u_2, F, err_laf, len,inliers, S.I);
        for (j = 0; j < S.I; j++)
            if (err_laf[inliers[j]] > th_laf_check)
                inl[j] = 0;
    }


cleanup:
    /* deallocations */

#ifdef __HASHING__
    htClear(&HASH_TABLE_F);
#endif // __HASHING__
    free(lo_workspace);
    free(workspace);
    if (data_out != NULL) {
        *data_out = no_sam;
        data_out[1] = iter_cnt;
    }

    ////printf("__PROFILE: AFTER ransac: %d\n", getticks()/1000);
    return maxS.I;

allocation_failure:
    memset(F, 0, 9 * sizeof(double));
    memset(inl, 0, len * sizeof(unsigned char));
    if (resids != NULL) {
        free(*resids);
        *resids = NULL;
    }
    no_sam = 0;
    iter_cnt = 0;
    maxS.I = 0;
    goto cleanup;

}
