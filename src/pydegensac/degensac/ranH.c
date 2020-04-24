//#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <limits.h>
#include <math.h>

#include "../matutls/matutl.h"
#include "utools.h"
#include "Htools.h"

#include "ranH.h"

#define SYM_COEF 2
#define DO_SYMMETRY_CHECK
#define FULL_SYMM1


Score iterH (double *u, int len, int *inliers, double th, double ths,
			double *H, double *Z, double **errs, double *buffer, unsigned inlLimit) {
	double *d = errs[1];
	double h[9], dth;
	int it, *inlSubset;
	Score S = {0,0}, Ss, maxS;


	dth = (ths - th) / ILSQ_ITERS; 

	/* H from the sample inliers by th */
	maxS = inlidxs(errs[4], len, th, inliers);
	if (maxS.I < 4) {
		return S;
	}
	if (maxS.I <= inlLimit) { /* if we are under the limit, just use what we have without shuffling */
		u2h(u, inliers, maxS.I, h, buffer);
	} else {
		inlSubset = randsubset (inliers, maxS.I, inlLimit);
		u2h(u, inlSubset, inlLimit, h, buffer);
	}

	/*iterate */
	for (it = 0; it < ILSQ_ITERS; ++it) {
#ifdef FULL_SYMM
	    HDsSym(Z, u, h, d, len);
#else
	    HDs (Z, u, h, d, len);
#endif
	    S = inlidxs(d, len, th, inliers);
		Ss = inlidxs(d, len, ths, inliers);

		if (scoreLess(maxS, S)) {
			maxS = S;
			errs[1] = errs[0];
			errs[0] = d;
			d = errs[1];
			memcpy(H, h, 9*sizeof(double));
		}
		if (Ss.I < 4) {

			return maxS;
		}

		if (Ss.I <= inlLimit) { /* if we are under the limit, just use what we have without shuffling */
			u2h(u, inliers, Ss.I, h, buffer);
		} else {
			inlSubset = randsubset (inliers, Ss.I, inlLimit);
			u2h(u, inlSubset, inlLimit, h, buffer);
		}

		ths -= dth;
	}
#ifdef FULL_SYMM
	HDsSym (Z, u, h, d, len);
#else
	HDs (Z, u, h, d, len);
#endif
	S = inlidxs (d, len, th, inliers);
	if (scoreLess(maxS, S)) {
		maxS = S;
		errs[1] = errs[0];
		errs[0] = d;
		memcpy(H, h, 9*sizeof(double));
	}

	return maxS;
}


Score inHrani (double *u, int len, int *inliers, int ninl, double th, double *Z,
			double **errs, double *buffer, double *H, unsigned inlLimit) {
	int ssiz, i;
	Score S, maxS = {0,0};
	double *d, h[9];
	int *sample;
	int *intbuff;

	if (ninl < 8) {
		return maxS;
	}

	intbuff = (int *) malloc (len * sizeof(int));

	ssiz = ninl / 2;
	if (ssiz > 12) {
		ssiz = 12;
	}

	d = errs[2];
	errs[2] = errs[0];
	errs[0] = d;

	for (i = 0; i < RAN_REP; ++i) {
		sample = randsubset(inliers, ninl, ssiz);
		u2h(u, sample, ssiz, h, buffer);
#ifdef FULL_SYMM
		HDsSym (Z, u, h, errs[0], len);
#else
		HDs (Z, u, h, errs[0], len);
#endif
		errs[4] = errs[0];
		S = iterH(u, len, intbuff, th, TC*th, h, Z, errs, buffer, inlLimit);
		if (scoreLess(maxS, S)) {
			maxS = S;
			d = errs[2];
			errs[2] = errs[0];
			errs[0] = d;
			memcpy(H, h, 9*sizeof(double));
		}
	}
	d = errs[2];
	errs[2] = errs[0];
	errs[0] = d;

	free(intbuff);
	return maxS;
}


Score ransacH (double *u, int len, double th, double conf, int max_sam,
			double *H, unsigned char * inl, int *data_out, int do_lo, int inlLimit) {
	int *pool, no_sam, new_sam, *samidx;
	double *Z, *buffer;
	double *err, *d, *h;
	double *errs[5];
	double *dsym;
	int i, j, *inliers,*inliers_sym, new_max, do_iterate, iter_cnt = 0, rej_cnt = 0;
	Score maxS = {0,0}, maxSs = {0,0}, S,Ssym;
	unsigned seed;

	double tol, v, M[9*9], sol[9*9];
	int nullspace_buff[2*9], nullsize;

	if (inlLimit == 0) { /* in the case of unlimited least squares */
		inlLimit = INT_MAX;
	}
	h = sol;


	/* allocations */
	pool = (int *)malloc(len * sizeof(int));
	for (i = 0; i < len; i++) {
		pool[i] = i;
	}
	samidx = pool + len - 4; /* drawn sample (indexes) is moved to the back of the pool */

	Z = (double *) malloc(len * 18 * sizeof(double));
	lin_hg(u, Z, pool, len);

	buffer = (double *) malloc(len * 18 * sizeof(double));

	err = (double *) malloc(len * 4 * sizeof(double));
	dsym = (double *) malloc(len * 1 * sizeof(double));

	for (i = 0; i < 4; i++) {
		errs[i] = err + i * len;
	}
	errs[4] = errs[3];

	inliers = (int *) malloc(len * sizeof(int));
	inliers_sym = (int *) malloc(len * sizeof(int));

	no_sam = 0;
	seed = rand();

	/* main RANSAC loop */
	while(no_sam < max_sam) {
		no_sam ++;
		new_max = 0; do_iterate = 0;

		srand(seed); /* to keep the same samples regardless any random sampling in the LO */

		/* random minimal sample */
		multirsampleT(Z, 9, 2, pool, 4, len, M); /* nullspace function expects M row-wise, thus 'T' */

		seed = rand();

		/* orientation check */
		if (!all_Hori_valid (u, samidx)) {
			++rej_cnt;
			continue;
		}

		/* use LU */
		for (i = 9*8; i < 9*9; ++i) { /* Fill with zeros to square */
			M[i] = 0.0;
		}
		nullsize = nullspace(M, sol, 9, nullspace_buff);
		if (nullsize != 1) {
			continue;
		}

		v = det3(h); tol = h[8];
		if (tol == 0) {
			for (i = 0; i < 9; ++i) { /* Frobenius norm */
				tol += h[i]*h[i];
			}
			tol = sqrt(tol);
			tol *= 0.001; /* typical ratio H(3,3)/||H||_F */
		}
		tol = tol*tol*tol;
		if (fabs(v/tol) < 10e-2) {
			continue;
		}


		/* consensus */
		d = errs[0];
#ifdef FULL_SYMM
		HDsSym(Z, u, h, d, len);
#else
		HDs(Z, u, h, d, len);
#endif
		S = inlidxs(d, len, th, inliers);

		if (scoreLess(maxS, S)) { /* so-far-the-best */
#ifdef DO_SYMMETRY_CHECK1 // Mishkin. degeneracy test by symmetric geom.error
	//	 printf("Number of inliers before test %d. MaxS = %d\n",S.I, maxS.I);
		 HDsSym(Z, u, h, dsym, len);
		 Ssym = inlidxs(dsym, len, th*SYM_COEF, inliers_sym);
	//	 printf("Number of inliers after test %d. MaxS = %d\n",Ssym.I, maxS.I);
		 if (!scoreLess(maxS, Ssym)) continue;
		 S = Ssym;
#endif
		maxS = S;
		errs[0] = errs[3];
		errs[3] = d;
		memcpy(H,h,9*sizeof(double));
		new_max = 1;
		}
		if (scoreLess(maxSs, S)) { /* so-far-the-best from sample */
			maxSs = S;
			do_iterate = no_sam > ITER_SAM;
			errs[4] = d;
		}

		if (no_sam >= ITER_SAM && iter_cnt == 0 && maxSs.I > 4) { /* after blocking, run LO on sftb sample */
			do_iterate = 1;
		}

		/* Local Optimisation */
		if (do_iterate && do_lo) {
			iter_cnt ++;
			/*******/
			/* minimalistic LO' (just one iterations) */
			/* S = iterH(u, len, inliers, th, TC*th, 4, h, Z, errs, buffer, inlLimit); */
			/*******/
			/* full LO (subsampling and iterations) */
			d = errs[0];
			S = inlidxs(errs[4], len, TC*th, inliers);
			u2h(u, inliers, S.I, h, buffer);
#ifdef FULL_SYMM
			HDsSym(Z, u, h, d, len);
#else
			HDs(Z, u, h, d, len);
#endif

			S = inlidxs(d, len, th, inliers);
			S = inHrani (u, len, inliers, S.I, th, Z, errs, buffer, h, inlLimit);
			/*******/
			tol = h[8];
			if (tol == 0) {
				for (i = 0; i < 9; ++i) { /* Frobenius norm */
					tol += h[i]*h[i];
				}
				tol = sqrt(tol);
				tol *= 0.001; /* typical ratio H(3,3)/||H||_F */
			}
			tol = tol*tol*tol;
			if (scoreLess(maxS, S) && (fabs(det3(h)/tol) > 10e-2)) {
				maxS = S;
				d = errs[0];
				errs[0] = errs[3];
				errs[3] = d;
				memcpy(H, h, 9*sizeof(double));
				new_max = 1;
			}
		}

		if (new_max) { /* updating number of samples needed */
			new_sam = nsamples(maxS.I+1, len, 4, conf);
			if (new_sam < max_sam) {
				max_sam = new_sam;
			}
		}
	}

	/* If there were no LOs, do at least one NOW! */
	if (do_lo && !iter_cnt) {
		++iter_cnt;
		/*******/
		/* minimalistic LO' (just one iterations) */
		/* S = iterH(u, len, inliers, th, TC*th, 4, h, Z, errs, buffer, inlLimit); */
		/*******/
		/* full LO (subsampling and iterations) */
		d = errs[0];
		S = inlidxs(errs[4], len, TC*th, inliers);
		u2h(u, inliers, S.I, h, buffer);
#ifdef FULL_SYMM
		HDsSym(Z, u, h, d, len);
#else
		HDs(Z, u, h, d, len);
#endif
		S = inlidxs(d, len, th, inliers);
		S = inHrani (u, len, inliers, S.I, th, Z, errs, buffer, h, inlLimit);
		/*******/
		tol = h[8];
		if (tol == 0) {
			for (i = 0; i < 9; ++i) { /* Frobenius norm */
				tol += h[i]*h[i];
			}
			tol = sqrt(tol);
			tol *= 0.001; /* typical ratio H(3,3)/||H||_F */
		}
		tol = tol*tol*tol;
		if(scoreLess(maxS, S) && (fabs(det3(h)/tol) > 10e-2)) {
			maxS = S;
			d = errs[0];
			errs[0] = errs[3];
			errs[3] = d;
			memcpy(H, h, 9*sizeof(double));
		}
	}

	if (inl) { /* set output field of inliers (binary this time) */
		d = errs[3];
		for (j = 0; j < len; j++) {
			if (d[j] <= th) {
				inl[j] = 1;
			} else {
				inl[j] = 0;
			}
		}
	}

	if (data_out) {
		data_out[0] = no_sam;
		data_out[1] = iter_cnt;
		data_out[2] = rej_cnt;
	}
	/* deallocations */
	free(pool);
	free(Z);
	free(buffer);
	free(err);
	free(inliers);
	free(inliers_sym);
	free(dsym);
	return maxS;
}

void ransacHsimple (double *u, int len, double th, double *H) {
	/* default settings, LO turned on with default inlier limit, only H returned */
	ransacH (u, len, th, CONFIDENCE, MAX_SAMPLES, H, 0, 0, 1, INL_LIMIT_H);
}

