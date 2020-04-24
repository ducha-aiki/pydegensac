#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <memory.h>
#include "lapwrap.h"
/*#include <mex.h> */

#include "../matutls/matutl.h"
#include "Htools.h"
#include "rtools.h"
#include "utools.h"
#include "ranH.h"


#include "ranH2el.h"


Score ransacH2el (double *u10, int len, double th, double conf, int max_sam,
			double *H, unsigned char * inl, int *data_out, int do_lo, int inlLimit) {
	int *pool, no_sam, new_sam, *samidx;
	double *u6, *Z, *buffer;
	double *err, *d, h[9];
	double *errs[5];
	int i, j, *inliers, new_max, do_iterate, iter_cnt = 0, rej_cnt = 0;
	Score maxS = {0,0}, maxSs = {0,0}, S;
	unsigned seed;
	double N1[9], D1[9], N2[9], D2[9]; /* stored column-wise! */
	
	double tol, v;
	
	
	if (inlLimit == 0) { /* in the case of unlimited least squares */
		inlLimit = INT_MAX;
	}

	/* allocations */
	u6 = (double *) malloc(6 * len * sizeof(double));
	for (i = 0; i < len; ++i) {
		u6[i*6 + 0] = u10[i*10 + 0];
		u6[i*6 + 1] = u10[i*10 + 1];
		u6[i*6 + 2] = 1;
		u6[i*6 + 3] = u10[i*10 + 5];
		u6[i*6 + 4] = u10[i*10 + 6];
		u6[i*6 + 5] = 1;
	}
	
	pool = (int *)malloc(len * sizeof(int));
	for (i = 0; i < len; i++) {
		pool[i] = i;
	}
	samidx = pool + len - 2; /* drawn sample (indexes) is moved to the back of the pool */

	Z = (double *) malloc(len * 18 * sizeof(double));
	lin_hg(u6, Z, pool, len);

	buffer = (double *) malloc(len * 18 * sizeof(double));

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

		srand(seed);

		/* random minimal sample */
		randsubset(pool, len, 2);
		
		seed = rand();

		/* model */
		getTransf(u10 + 10*samidx[0], N1, D1);
		getTransf(u10 + 10*samidx[1], N2, D2);
		if (A2toRH(N1, D1, N2, D2, u10, samidx, h)) {
			continue;
		}
		
		v = det3(h); tol = h[8]; tol = tol*tol*tol; //TODO if tol == 0
		if (fabs(v/tol) < 10e-2) {
			continue;
		}


		/* consensus */
		d = errs[0];
		HDs(Z, u6, h, d, len);
		S = inlidxs(d, len, th, inliers);

		if (scoreLess(maxS, S)) { /* so-far-the-best */
			maxS = S;
			errs[0] = errs[3];
			errs[3] = d;
			memcpy(H,h,9*sizeof(double));
			new_max = 1;
		}
		
		S = inlidxs(d, len, th*TAU, inliers);
		if (scoreLess(maxSs, S)) { /* so-far-the-best from sample */
			maxSs = S;
			do_iterate = no_sam > ITER_SAM;
			if (!new_max) {
				errs[0] = errs[2];
				errs[2] = d;
			}
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
			/* S = iterH(u6, len, inliers, th, TC*th, 4, h, Z, errs, buffer, inlLimit); */
			/*******/
			/* full LO (subsampling and iterations) */
			d = errs[0];
			S = inlidxs(errs[4], len, TC*th*TAU, inliers);
			u2h(u6, inliers, S.I, h, buffer);
			HDs(Z, u6, h, d, len);
			S = inlidxs(d, len, th, inliers);
			S = inHraniEl (u10, u6, len, inliers, S.I, th, Z, errs, buffer, h, inlLimit);
			/*******/
			tol = h[8]; tol = tol*tol*tol;
			if (scoreLess(maxS, S) && (fabs(det3(h)/tol) > 10e-2)) {
				maxS = S;
				d = errs[0];
				errs[0] = errs[3];
				errs[3] = d;
				memcpy(H, h, 9*sizeof(double));
				new_max = 1;
			}
		}

		if (new_max) { /* update number of samples needed */
			new_sam = nsamples(maxS.I+1, len, 2, conf);
			if (new_sam < max_sam) {
				max_sam = new_sam;
			}
		}
	}

	/* If there were no LO's, make at least one NOW! */
	if (do_lo && !iter_cnt) {
		++iter_cnt;
		/*******/
		/* minimalistic LO' (just one iterations) */
		/* S = iterH(u6, len, inliers, th, TC*th, 4, h, Z, errs, buffer, inlLimit); */
		/*******/
		/* full LO (subsampling and iterations) */
		d = errs[0];
		S = inlidxs(errs[4], len, TC*th*TAU, inliers);
		u2h(u6, inliers, S.I, h, buffer);
		HDs(Z, u6, h, d, len);
		S = inlidxs(d, len, th, inliers);
		S = inHraniEl (u10, u6, len, inliers, S.I, th, Z, errs, buffer, h, inlLimit);
		/*******/
		tol = h[8]; tol = tol*tol*tol;
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
	free(u6);

	return maxS;
}


void getTransf (double *u10, double *N, double *D) {
	D[0] = u10[2]; /* a */
	D[1] = u10[3]; /* b */
	D[2] = 0;
	D[3] = 0;
	D[4] = u10[4]; /* c */
	D[5] = 0;
	D[6] = u10[0]; /* x */
	D[7] = u10[1]; /* y */
	D[8] = 1;
	N[0] = 1 / u10[7]; /* 1/a */
	N[1] = - u10[8] / u10[7] / u10[9]; /* -b/(ac) */
	N[2] = 0;
	N[3] = 0;
	N[4] = 1 / u10[9]; /* 1/c */
	N[5] = 0;
	N[6] = - u10[5] / u10[7]; /* -x/a */
	N[7] = (u10[8]*u10[5] - u10[7]*u10[6]) / u10[7] / u10[9]; /* (bx - ay)/(ac) */
	N[8] = 1;
}

int A2toRH(double *N1, double *D1, double *N2, double *D2, double *u, int *samidx, double *h) {
	int do_norm = 0, i;
	double Z[15*15], ZT[15*15]; /* everything here is stored column-wise, unless noted */
	double U[15*15];
	double T1[3*3], iT1[3*3], T2[3*3], iT2[3*3];
	double N1N[3*3], N2N[3*3], D1N[3*3], D2N[3*3], temp[3*3];
	int nullsize, nullspace_buff[2*15];
	
	for (i = 0; i < 2*7*15; ++i) {
		Z[i] = 0.0;
	}
	
	if (do_norm) {
		norm2pt(u[10*samidx[0] + 0], u[10*samidx[0] + 1], u[10*samidx[1] + 0], u[10*samidx[1] + 1], T1, iT1);
		norm2pt(u[10*samidx[0] + 5], u[10*samidx[0] + 6], u[10*samidx[1] + 5], u[10*samidx[1] + 6], T2, iT2);
		ZuN(Z, u + 10*samidx[0], T1, T2, 2);
		ZuN(Z+7, u + 10*samidx[1], T1, T2, 2);
		mmul(N1N, iT2, N1, 3); /* CCMath works with row-wise stored matrices, so we use reverted order: A = B*C -> A^T = C^T*B^T */
		mmul(N2N, iT2, N2, 3);
		mmul(D1N, D1, T1, 3);
		mmul(D2N, D2, T1, 3);
		Znd(Z + 2*7*9, D1N, N1N, 2);
		Znd(Z + 2*7*9 + 2*7*3 + 7, D2N, N2N, 2);
	} else {
		Zu(Z, u + 10*samidx[0], 2);
		Zu(Z+7, u + 10*samidx[1], 2);
		Znd(Z + 2*7*9, D1, N1, 2);
		Znd(Z + 2*7*9 + 2*7*3 + 7, D2, N2, 2);
	}

	mattr(ZT, Z, 15, 2*7);
	for (i = 14*15; i < 15*15; ++i) {
		ZT[i] = 0;
	}
	nullsize = nullspace(ZT, U, 15, nullspace_buff);
	memcpy(h, U, 3*3 * sizeof(double));
	trnm(h, 3); /* Equations are made for H stored row-wise, so transpose now */
	
	if (do_norm) {
		/* Hn = iT1 * H * T2 -> Hn^T = T2^T * H^T * iT1^T */
		mmul(temp, T2, h, 3);
		mmul(h, temp, iT1, 3);
	}
	
	return nullsize != 1;
}

int AntoRH(double *u, int *inls, int len, double *h) {
	if (len < 2) {
		return 1;
	}
	double Ni[3*3], Di[3*3];
	if (len == 2) {
		double N2[3*3], D2[3*3];
		getTransf(u + 10*inls[0], Ni, Di);
		getTransf(u + 10*inls[1], N2, D2);
		return A2toRH(Ni, Di, N2, D2, u, inls, h);
	}
	int do_norm = 1, i;//, k, l;
	double *Z; /* everything here is stored column-wise, unless noted */
	double *U, *VT, *D;
	double T1[3*3], iT1[3*3], T2[3*3], iT2[3*3];
	double Nn[3*3], Dn[3*3], temp[3*3];
	int res;

	int noRows = 7 * len;
	int noCols = 9 + 3 * len;

	Z = (double *) malloc (noRows * noCols * sizeof(double));
	U = (double *) malloc (noRows * noRows * sizeof(double));
	VT = (double *) malloc (noCols * noCols * sizeof(double));
	D = (double *) malloc (noCols * sizeof(double));
	
	for (i = 0; i < noRows * noCols; ++i) {
		Z[i] = 0.0;
	}
	
	if (do_norm) {
		norm10(u, inls, len, T1, iT1);
		norm10(u+5, inls, len, T2, iT2);
	}
	
	for (i = 0; i < len; ++i) {
		getTransf(u + 10*inls[i], Ni, Di);
		/*printf("\n%dth correspondence:\n", inls[i]);
		printf("N_%d:\n", i);
		for (k = 0; k < 3; ++k) {
			for (l = 0; l < 3; ++l) {
				printf("%9.4f", Ni[k + l*3]);
			}
			printf("\n");
		}
		printf("D_%d:\n", i);
		for (k = 0; k < 3; ++k) {
			for (l = 0; l < 3; ++l) {
				printf("%9.4f", Di[k + l*3]);
			}
			printf("\n");
		}*/
		if (do_norm) {
			mmul(Nn, iT2, Ni, 3); /* CCMath works with row-wise stored matrices, so we use reverted order: A = B*C -> A^T = C^T*B^T */
			mmul(Dn, Di, T1, 3);
			ZuN(Z + 7*i, u + 10*inls[i], T1, T2, len);
			Znd(Z + len*7*9 + len*7*3*i + 7*i, Dn, Nn, len);
		} else {
			Zu(Z + 7*i, u + 10*inls[i], len);
			Znd(Z + len*7*9 + len*7*3*i + 7*i, Di, Ni, len);
		}
	}

	res = lap_SVD (D, Z, U, noRows, VT, noCols);

	for (i = 0; i < 9; ++i) {
		h[i] = VT[(i+1)*noCols - 1];
	}
	trnm(h, 3);
	
	if (do_norm) {
		/* Hn = iT1 * H * T2 -> Hn^T = T2^T * H^T * iT1^T */
		mmul(temp, T2, h, 3);
		mmul(h, temp, iT1, 3);
	}
	
	free(Z);
	free(U);
	free(D);
	free(VT);
	return res;
}

void Zu(double * Z, double * u, int len) {
	Z[0 + 0*len*7] = -1;
	Z[0 + 6*len*7] = _u1;
	Z[1 + 1*len*7] = -1;
	Z[1 + 7*len*7] = _u1;
	Z[2 + 2*len*7] = -1;
	Z[2 + 6*len*7] = - _u1 * _u4;
	Z[2 + 7*len*7] = - _u1 * _u5;
	Z[3 + 3*len*7] = -1;
	Z[3 + 6*len*7] = _u2;
	Z[4 + 4*len*7] = -1;
	Z[4 + 7*len*7] = _u2;
	Z[5 + 5*len*7] = -1;
	Z[5 + 6*len*7] = - _u2 * _u4;
	Z[5 + 7*len*7] = - _u2 * _u5;
	Z[6 + 8*len*7] = -1;
	Z[6 + 6*len*7] = - _u4;
	Z[6 + 7*len*7] = - _u5;
}

void ZuN(double * Z, double * u, double * T1, double * T2, int len) {
	Z[0 + 0*len*7] = -1;
	Z[0 + 6*len*7] = _u1*T1[0] + _u2*T1[3] + T1[6];
	Z[1 + 1*len*7] = -1;
	Z[1 + 7*len*7] = _u1*T1[0] + _u2*T1[3] + T1[6];
	Z[2 + 2*len*7] = -1;
	Z[2 + 6*len*7] = - (_u1*T1[0] + _u2*T1[3] + T1[6]) * (_u4*T2[0] + _u5*T2[3] + T2[6]);
	Z[2 + 7*len*7] = - (_u1*T1[0] + _u2*T1[3] + T1[6]) * (_u4*T2[1] + _u5*T2[4] + T2[7]);
	Z[3 + 3*len*7] = -1;
	Z[3 + 6*len*7] = _u1*T1[1] + _u2*T1[4] + T1[7];
	Z[4 + 4*len*7] = -1;
	Z[4 + 7*len*7] = _u1*T1[1] + _u2*T1[4] + T1[7];
	Z[5 + 5*len*7] = -1;
	Z[5 + 6*len*7] = - (_u1*T1[1] + _u2*T1[4] + T1[7]) * (_u4*T2[0] + _u5*T2[3] + T2[6]);
	Z[5 + 7*len*7] = - (_u1*T1[1] + _u2*T1[4] + T1[7]) * (_u4*T2[1] + _u5*T2[4] + T2[7]);
	Z[6 + 8*len*7] = -1;
	Z[6 + 6*len*7] = - (_u4*T2[0] + _u5*T2[3] + T2[6]);
	Z[6 + 7*len*7] = - (_u4*T2[1] + _u5*T2[4] + T2[7]);
}

void Znd(double * Z, double * A, double * B, int len) {
	/* A&B transposed by #defines! */
	Z[2 + 2*len*7] = _a3;
	Z[5 + 2*len*7] = _a6;
	Z[6 + 2*len*7] = 1;
	Z[0 + 0*len*7] = _a2*_b1 - _a1*_b4;
	Z[1 + 0*len*7] = _a2*_b2 - _a1*_b5;
	Z[2 + 0*len*7] = _a2*_b3 - _a1*_b6;
	Z[3 + 0*len*7] = _a5*_b1 - _a4*_b4;
	Z[4 + 0*len*7] = _a5*_b2 - _a4*_b5;
	Z[5 + 0*len*7] = _a5*_b3 - _a4*_b6;
	Z[0 + 1*len*7] = _a1*_b1 + _a2*_b4;
	Z[1 + 1*len*7] = _a1*_b2 + _a2*_b5;
	Z[2 + 1*len*7] = _a1*_b3 + _a2*_b6;
	Z[3 + 1*len*7] = _a4*_b1 + _a5*_b4;
	Z[4 + 1*len*7] = _a4*_b2 + _a5*_b5;
	Z[5 + 1*len*7] = _a4*_b3 + _a5*_b6;
}

void norm2pt(double x1, double y1, double x2, double y2, double *T, double *iT) {
	double xm = (x1 + x2) / 2;
	double ym = (y1 + y2) / 2;
	double dx = (x1 - x2) / 2;
	double dy = (y1 - y2) / 2;
	double sc = sqrt(dx*dx + dy*dy);
	
	if (sc < 1) {
		sc = 1;
	}
	
	iT[0] = sc;
	iT[1] = 0;
	iT[2] = 0;
	iT[3] = 0;
	iT[4] = sc;
	iT[5] = 0;
	iT[6] = xm;
	iT[7] = ym;
	iT[8] = 1;
	
	sc = 1 / sc;
	
	T[0] = sc;
	T[1] = 0;
	T[2] = 0;
	T[3] = 0;
	T[4] = sc;
	T[5] = 0;
	T[6] = - xm * sc;
	T[7] = - ym * sc;
	T[8] = 1;
}

void norm10(double * u10, int * inls, int len, double * T, double * iT) {
	double xm = 0, ym = 0, sc = 0;
	int i;
	
	for (i = 0; i < len; ++i) {
		xm += u10[10 * inls[i]] / len;
		ym += u10[10 * inls[i] + 1] / len;
	}
	
	for (i = 0; i < len; ++i) {
		sc += sqrt((u10[10*inls[i]] - xm) * (u10[10*inls[i]] - xm) + (u10[10*inls[i]+1] - ym) * (u10[10*inls[i]+1] - ym)) / len;
	}
	sc /= SQRT2;
	
	iT[0] = sc;
	iT[1] = 0;
	iT[2] = 0;
	iT[3] = 0;
	iT[4] = sc;
	iT[5] = 0;
	iT[6] = xm;
	iT[7] = ym;
	iT[8] = 1;
	
	sc = 1 / sc;
	
	T[0] = sc;
	T[1] = 0;
	T[2] = 0;
	T[3] = 0;
	T[4] = sc;
	T[5] = 0;
	T[6] = - xm * sc;
	T[7] = - ym * sc;
	T[8] = 1;
}

Score inHraniEl (double * u10, double *u6, int len, int *inliers, int ninl, double th, double *Z,
			double **errs, double *buffer, double *H, unsigned inlLimit) {
	int ssiz, i;
	Score S, maxS = {0,0};
	double *d, h[9];
	int *sample;
	int *intbuff;
	int minPts = 4, loLimit = 8;

	intbuff = (int *) malloc (len * sizeof(int));

	if (ninl < loLimit) {
		return maxS;
	}
	ssiz = ninl / 2;
	if (ssiz > 12) {
		ssiz = 12;
	}

	d = errs[2];
	errs[2] = errs[0];
	errs[0] = d;

	for (i = 0; i < RAN_REP; ++i) {
		sample = randsubset(inliers, ninl, ssiz);
		if (ssiz < minPts) {
			AntoRH(u10, sample, ssiz, h);
		} else {
			u2h(u6, sample, ssiz, h, buffer);
		}
		HDs (Z, u6, h, errs[0], len);
		errs[4] = errs[0];

		S = iterH(u6, len, intbuff, th, TC*th, h, Z, errs, buffer, inlLimit);

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







