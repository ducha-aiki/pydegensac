////#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <limits.h>


#include "../matutls/matutl.h"
#include "utools.h"
#include "Ftools.h"

#include "ranF.h"


Score iterF (double *u, int len, int *inliers, double th, double ths,
			double *F, double **errs, double *buffer, unsigned inlLimit) {
	double *d = errs[1], *w;
	double f[9], dth;
	int it, *inlSubset;
	Score S = {0,0}, Ss, maxS;

	w = (double *) malloc(len * sizeof(double));
	dth = (ths - th) / ILSQ_ITERS; 

	/* F from the sample inliers by th */
	maxS = inlidxs(errs[4], len, th, inliers);
	if (maxS.I < 8) {
		free(w);
		return S;
	}
	if (maxS.I <= inlLimit) { /* if we are under the limit, just use what we have without shuffling */
		u2f(u, inliers, maxS.I, f, buffer);
	} else {
		inlSubset = randsubset (inliers, maxS.I, inlLimit);
		u2f(u, inlSubset, inlLimit, f, buffer);
	}

	/*iterate */
	for (it = 0; it < ILSQ_ITERS; it ++) {
		exFDs (u, f, d, w, len);
		S = inlidxs(d, len, th, inliers);
		Ss = inlidxs(d, len, ths, inliers);

		if (scoreLess(maxS, S)) {
			maxS = S;
			errs[1] = errs[0];
			errs[0] = d;
			d = errs[1];
			memcpy(F, f, 9*sizeof(double));
		}
		if (Ss.I < 8) {
			free(w);
			return maxS;
		}

		if (Ss.I <= inlLimit) { /* if we are under the limit, just use what we have without shuffling */
			u2fw(u, inliers, w, Ss.I, f, buffer);
		} else {
			inlSubset = randsubset (inliers, Ss.I, inlLimit);
			u2fw(u, inlSubset, w, inlLimit, f, buffer);
		}

		ths -= dth;
	}

	FDs (u, f, d, len);
	S = inlidxs (d, len, th, inliers);
	if (scoreLess(maxS, S)) {
		maxS = S;
		errs[1] = errs[0];
		errs[0] = d;
		memcpy(F, f, 9*sizeof(double));
	}
	free(w);
	return maxS;
}


Score inFrani (double *u, int len, int *inliers, int ninl, double th,
			double **errs, double *buffer, double *F, unsigned inlLimit) {
	int ssiz, i;
	Score S, maxS = {0,0};
	double *d, f[9];
	int *sample;
	int *intbuff;

	if (ninl < 16) {
		return maxS;
	}

	intbuff = (int *) malloc (len * sizeof(int));

	ssiz = ninl / 2;
	if (ssiz > 14) {
		ssiz = 14;
	}

	d = errs[2];
	errs[2] = errs[0];
	errs[0] = d;

	for (i = 0; i < RAN_REP; ++i) {
		sample = randsubset(inliers, ninl, ssiz);
		u2f(u, sample, ssiz, f, buffer);
		FDs (u, f, errs[0], len);
		errs[4] = errs[0];

		S = iterF(u, len, intbuff, th, TC*th, f, errs, buffer, inlLimit);

		if (scoreLess(maxS, S)) {
			maxS = S;
			d = errs[2];
			errs[2] = errs[0];
			errs[0] = d;
			memcpy(F, f, 9*sizeof(double));
		}
	}
	d = errs[2];
	errs[2] = errs[0];
	errs[0] = d;

	free(intbuff);
	return maxS;
}


Score ransacF (double *u, int len, double th, double conf, int max_sam,
			double *F, unsigned char * inl, int *data_out, int do_lo, int inlLimit) {
	int *pool, no_sam, new_sam, *samidx;
	double *Z, *buffer;
	double *err, *d, f[9];
	double *errs[5];
	int i, j, *inliers, new_max, do_iterate, iter_cnt = 0, rej_cnt = 0, nsol;
	Score maxS = {0,0}, maxSs = {0,0}, S;
	unsigned seed;

	double *f1, *f2, poly[4], roots[3], A[9*9], sol[9*9];
	int nullspace_buff[2*9], nullsize;

	if (inlLimit == 0) { /* in the case of unlimited least squares */
		inlLimit = INT_MAX;
	}
	f1 = sol;
	f2 = sol+9;

	/* allocations */
	pool = (int *)malloc(len * sizeof(int));
	for (i = 0; i < len; i++) {
		pool[i] = i;
	}
	samidx = pool + len - 7; /* drawn sample (indexes) is moved to the back of the pool */

	Z = (double *) malloc(len * 9 * sizeof(double));
	lin_fm(u, Z, pool, len);

	buffer = (double *) malloc(len * 9 * sizeof(double));

	err = (double *) malloc(len * 4 * sizeof(double));
	for (i = 0; i < 4; i++) {
		errs[i] = err + i * len;
	}
	errs[4] = errs[3];

	inliers = (int *) malloc(len * sizeof(int));

	no_sam = 0;
	seed = rand();

	/* main RANSAC loop */
	while(no_sam < max_sam) {
		no_sam ++;
		new_max = 0; do_iterate = 0;

		srand(seed); /* to keep the same samples regardless any random sampling in the LO */

		/* random minimal sample */
		rsampleT(Z, 9, pool, 7, len, A);

		seed = rand();

		/* use LU */
		for (i = 7*9; i < 9*9; ++i) { /* Fill with zeros to square */
			A[i] = 0.0;
		}
		nullsize = nullspace(A, f1, 9, nullspace_buff);
		if (nullsize != 2) {
			continue;
		}

		slcm (f1, f2, poly);  
		nsol = rroots3(poly, roots);

		for (i = 0; i < nsol; i++) { /* 1 or 3 hypotheses per sample */
			for (j = 0; j < 9; j++) {
				f[j] = f1[j] * roots[i] + f2[j] * (1 -roots[i]);
			}

			/* orient. constr. */
			if (!all_ori_valid(f, u, samidx, 7)) {
				++rej_cnt;
				continue;
			}

			/* consensus */
			d = errs[i];
			FDs(u, f, d, len);
			S = inlidxs(d, len, th, inliers);

			if(scoreLess(maxS, S)) { /* so-far-the-best */
				maxS = S;
				errs[i] = errs[3];
				errs[3] = d;
				memcpy(F,f,9*sizeof(double));
				new_max = 1;
			}
			if(scoreLess(maxSs, S)) { /* so-far-the-best from sample */
				maxSs = S;
				do_iterate = no_sam > ITER_SAM;
				errs[4] = d;
			}
		}

		if (no_sam >= ITER_SAM && iter_cnt == 0 && maxSs.I > 7) { /* after blocking, run LO on sftb sample */
			do_iterate = 1;
		}

		/* Local Optimisation */
		if (do_iterate && do_lo) {
			iter_cnt ++;
			/*******/
			/* minimalistic LO' (just one iterations) */
			/* S = iterF(u, len, inliers, th, 16*TC*th, f, errs, buffer, inlLimit);*/
			/*******/
			/* full LO (subsampling and iterations) */
			d = errs[0];
			S = inlidxs(errs[4], len, TC*th, inliers);
			u2f(u, inliers, S.I, f, buffer);
			FDs(u, f, d, len);
			S = inlidxs(d, len, th, inliers);
			S = inFrani(u, len, inliers, S.I, th, errs, buffer, f, inlLimit);
			/*******/

			if(scoreLess(maxS, S)) {
				maxS = S;
				d = errs[0];
				errs[0] = errs[3];
				errs[3] = d;
				memcpy(F, f, 9*sizeof(double));
				new_max = 1;
			}
		}

		if (new_max) { /* updating number of samples needed */
			new_sam = nsamples(maxS.I+1, len, 7, conf);
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
		/* S = iterF(u, len, inliers, intbuff2, th, 16*TC*th, f, errs, buffer, inlLimit); */
		/*******/
		/* full LO (subsampling and iterations) */
		d = errs[0];
		S = inlidxs(errs[4], len, TC*th, inliers);
		u2f(u, inliers, S.I, f, buffer);
		FDs(u, f, d, len);
		S = inlidxs(d, len, th, inliers);
		S = inFrani (u, len, inliers, S.I, th, errs, buffer, f, inlLimit);
		/*******/

		if(scoreLess(maxS, S)) {
			maxS = S;
			d = errs[0];
			errs[0] = errs[3];
			errs[3] = d;
			memcpy(F, f, 9*sizeof(double));
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

	data_out[0] = no_sam;
	data_out[1] = iter_cnt;
	data_out[2] = rej_cnt;

	/* deallocations */
	free(pool);
	free(Z);
	free(buffer);
	free(err);
	free(inliers);

	return maxS;
}

void ransacFsimple (double *u, int len, double th, double *F) {
	int data_out[3];
	/* default settings, LO turned on with default inlier limit, only F returned */
	ransacF (u, len, th, CONFIDENCE, MAX_SAMPLES, F, 0, data_out, 1, INL_LIMIT_F);
}

