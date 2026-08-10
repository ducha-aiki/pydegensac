//#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "../matutls/matutl.h"
#include "Htools.h"

#include "ranH.h"

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

