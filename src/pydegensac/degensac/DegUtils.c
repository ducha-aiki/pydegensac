#undef __STRICT_ANSI__
#include "DegUtils.h"
#include "rtools.h"
#include "Htools.h"
#include "Ftools.h"
#include "ranH.h"
#include "exp_ranF.h"
#include "matutls/matutl.h"
#include <stdlib.h>
#include <string.h>
//#include <stdio.h>
#include <math.h>
//#include <mex.h>
#include <time.h>

/*function [deg, H, inl] = checksample(F, u7, th)

%checksample tests for the degeneracy of 7pt sample
% [deg, H, inl] = checksample(F, u7, th)
% F is the fundamental matrix
% u7 6-by-7 matrix of 2 view correspondences, th - threshold
% when the seven corrs contain at least five point correspondences (inl)
% linked by homography (H) then deg is set to sum(inl)

IDXS = [1,2,3; 4,5,6; 1,2,7; 4,5,7; 3,6,7];
%IDXS = [IDXS; 2,3,4; 3,4,5; 1,3,5; 2,3,5];

for i = 1:5
  H = Hdetect(F,u7(:,IDXS(i,:)));
  Ds = fHDs(H, u7);
  [sDs, idx] = sort(Ds);
  H = u2H(u7(:,idx(1:5)));
  inl = fHDs(H,u7) < th;
  deg = sum(inl) > 4;
  if deg
    return
  end
end
*/

int checksample(double * F, double * u7, double th, double * H) {
  unsigned char IDXS[5][3] = {{0,1,2}, {3,4,5}, {0,1,6}, {3,4,6}, {2,5,6}};
  int i, j;
  double Ds[7], sDs[7];
  unsigned char idx[7], inlCount = 0;
  int inl[7];
  double buffer[5*18*sizeof(double)], bufferZ[7*18*sizeof(double)];
  int bufferP[7*sizeof(int)];

#ifndef __DEGEN__
  return 0;
#endif

  for (i = 0; i < 5; ++i) {
      //printf("Checking corrs: {%d,%d,%d}...\n",IDXS[i][0],IDXS[i][1],IDXS[i][2]);

      Hdetect(F, u7, IDXS[i], H);
      dHDs(H, u7, 7, Ds, bufferP, bufferZ);
      sortDs(Ds, sDs, idx);

      for (j = 0; j < 5; ++j) {
          inl[j] = idx[j];
        }

      u2h(u7, inl, 5, H, buffer);

      dHDs(H, u7, 7, Ds, bufferP, bufferZ);
      inlCount = 0;
      for (j = 0; j < 7; ++j) {
          /*printf("Error[%d] = %.2f %c %.2f\n", j, Ds[j], (Ds[j]>th?'>':'<'), th);*/
          if (Ds[j] < th) {
              ++inlCount;
            }
        }
      //printf("Found homography with %d inliers.\n", inlCount);
      if (inlCount > 4) {
          return 1;
        }
    }
  return 0;
}


/*function H = Hdetect (F, u3)

% Hdetect calcucates homography from fund. matrix F and 3 point corrs.
% H = Hdetect (F, u3)
% F :  3-by-3 rank 2 fundamental matrix
% u3 : 6-by-3 correspondences of image pts in homog. coordinates
% see Hartley & Zisserman: Scene planes and homographies (p.318)*/

void Hdetect(double * F, double * u7, unsigned char * IDXS, double * H) {
  double D[3], U[3*3], V[3*3], ec[3], Ex[3*3], A[3*3],
      u3a[3*3], u3b[3*3], u3aT[3*3], u3bT[3*3], Au3b[3*3],
      Ft[3*3], F1[3*3], p1[3*3], p1T[3*3], p2[3*3], b[3];
  unsigned char i, j, sing;

  /*[U,D,V] = svd(F');
        ec = V(:,3);
        Ex = skew_sym(ec);
        A  = Ex * F;*/

  /*Transpose F*/
  mattr(Ft,F,3,3);
  /*F would be altered during computation of SVD*/
  memcpy(F1, F, 3*3*sizeof(double));
  /*SVD - F should not be transposed, because it's stored column-wise*/
  svduv(D,F1,U,3,V,3);
  ec[0] = V[2]; ec[1] = V[5]; ec[2] = V[8];
  skew_sym(ec, Ex);
  /*Ft row-wise = F column-wise :-)*/
  mmul(A, Ex, Ft, 3);

  /*p1 = cross(u3(1:3,1:3), A * (u3(4:6,1:3)));
        p2 = -Ex * u3(1:3,1:3);
        b  = sum(p1 .* p2) ./ sum(p2.^2);*/

  fillu3(u7, IDXS, u3a, u3b);
  mmul(Au3b, A, u3b, 3);
  /*transpose for crossprod*/
  mattr(u3aT, u3a, 3, 3);
  mattr(u3bT, Au3b, 3, 3);
  crossp(p1T, u3aT, u3bT);
  crossp(p1T+3, u3aT+3, u3bT+3);
  crossp(p1T+6, u3aT+6, u3bT+6);
  mattr(p1, p1T, 3, 3);
  for (i = 0; i < 9; ++i) {
      Ex[i] *= -1;
    }
  mmul(p2, Ex, u3a, 3);
  b[0] = (p1[0]*p2[0] + p1[3]*p2[3] + p1[6]*p2[6]) / (p2[0]*p2[0] + p2[3]*p2[3] + p2[6]*p2[6]);
  b[1] = (p1[1]*p2[1] + p1[4]*p2[4] + p1[7]*p2[7]) / (p2[1]*p2[1] + p2[4]*p2[4] + p2[7]*p2[7]);
  b[2] = (p1[2]*p2[2] + p1[5]*p2[5] + p1[8]*p2[8]) / (p2[2]*p2[2] + p2[5]*p2[5] + p2[8]*p2[8]);

  /*M  = u3(4:6,1:3)';
        H = A - ec*(inv(M)*b(1:3)')';*/

  mattr(u3bT, u3b, 3, 3);
  /*What to do in singular case?*/
  sing = minv(u3bT, 3);
  /*Dim 3x1*/
  rmmult(u3b, u3bT, b, 3, 3, 1);
  mattr(u3bT, u3b, 3, 1);
  rmmult(u3b, ec, u3bT, 3, 1, 3);
  /*H must be stored column-wise*/
  for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
          H[i+j*3] = A[i*3+j] - u3b[i*3+j];
        }
    }

  /*if isnan(H(1)) || isinf(H(1))
                H = eye(3);
        end*/

  if ( isnan(*H) || isinf(*H) || sing ) {
      H[1] = H[2] = H[3] = H[5] = H[6] = H[7] = 0;
      H[0] = H[4] = H[8] = 1;
    }
}


void sortDs(double * Ds, double * sDs, unsigned char * idx){
  unsigned char i, j, auxI;
  double auxD;
  memcpy(sDs, Ds, 7 * sizeof(double));
  for (i = 0; i < 7; ++i) {
      idx[i] = i;
    }
  for (i = 0; i < 7; ++i) {
      for (j = i+1; j < 7; ++j) {
          if(sDs[j] < sDs[i]) {
              auxD = sDs[j];
              sDs[j] = sDs[i];
              sDs[i] = auxD;
              auxI = idx[j];
              idx[j] = idx[i];
              idx[i] = auxI;
            }
        }
    }
}


void dHDs(double * H, double * u, unsigned len, double * Ds, int * bufferP, double * bufferZ) {
  unsigned i;
  int * p;
  double * Z;

  if (bufferP) {
      p = bufferP;
    } else {
      p = (int *)malloc(len * sizeof(int));
    }
  if (bufferZ) {
      Z = bufferZ;
    } else {
      Z = (double *) malloc(len * 18 * sizeof(double));
    }
  for (i = 0; i < len; i ++) p[i] = i;
  lin_hg(u, Z, p, len);
  HDs(Z, u, H, Ds, len);
  if (!bufferP) {
      free(p);
    }
  if (!bufferZ) {
      free(Z);
    }
}

/*ax = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0]; */
void skew_sym(double * a, double * ax){
  ax[0] = 0;
  ax[1] = -a[2];
  ax[2] = a[1];
  ax[3] = a[2];
  ax[4] = 0;
  ax[5] = -a[0];
  ax[6] = -a[1];
  ax[7] = a[0];
  ax[8] = 0;
}


/*Prepares triplets of points for cross products*/
/**u7 stored column-wise, u3s row-wise*/
void fillu3(double * u7, unsigned char * IDXS, double * u3a, double * u3b){
  unsigned char i,j;
  for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
          u3a[i+j*3] = u7[IDXS[i]*6+j];
          u3b[i+j*3] = u7[IDXS[i]*6+j+3];
        }
    }
}


/*  crossp.c    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */
void crossp(double *h,double *u,double *v)
{ h[0]=u[1]*v[2]-u[2]*v[1];
  h[1]=u[2]*v[0]-u[0]*v[2];
  h[2]=u[0]*v[1]-u[1]*v[0];
}


/*function [F, inls] = rFtH(u, hinl, th, H)*/
unsigned rFtH(double * u, unsigned char * hinl, double th, double * H, unsigned len,
              double * F, int * bufferP, double * bufferZ) {

  unsigned char * nhinl, * v, * inl;
  double * Ds, * uN, * us, * uH, *uV;
  unsigned i, nhinlCount = 0, hinlCount = 0, ninl, maxni;
  unsigned * ptr, max_i, m_i, max_sam, s_size, pos, idx, auxI, no_i;
  double ec[3], ecNorm, c1[3], c2[3], aFt[3*3], aFtH[3*3], aF[3*3], Ht[3*3];

  /*MAX_SAM = 10000;
        conf = Ds[j] < th.999;
        sam_sizH = 6;
        sam_sizO = 4;*/
  unsigned MAX_SAM = 10000, no_sam;
  double conf = .999;
  unsigned char sam_sizH = 6;
  unsigned char sam_sizO = 4;

  /*nhinl = fHDs(H,u) > 100 * th;*/
  inl = (unsigned char *) malloc(len*sizeof(unsigned char));
  Ds = (double *) malloc(len*sizeof(double));
  dHDs(H, u, len, Ds, bufferP, bufferZ);
  nhinl = (unsigned char *) malloc(len*sizeof(unsigned char));
  for (i = 0; i < len; ++i) {
      if (Ds[i] > 100*th) {
          nhinl[i] = 1;
          ++nhinlCount;
        } else {
          nhinl[i] = 0;
        }
      if (hinl[i]) {
          ++hinlCount;
        }
    }
  free(Ds);
  Ds = (double *) malloc(nhinlCount*sizeof(double));
  v = (unsigned char *) malloc(nhinlCount*sizeof(unsigned char));

  /*uN = u(:,nhinl);
        us = [uN(1:3,:); H * uN(4:6,:)];
        uH = u(:,hinl);*/
  uN = (double *) malloc(6*nhinlCount*sizeof(double));
  us = (double *) malloc(6*nhinlCount*sizeof(double));
  uV = (double *) malloc(6*nhinlCount*sizeof(double));
  uH = (double *) malloc(6*hinlCount*sizeof(double));
  nhinlCount = 0, hinlCount = 0;
  for (i = 0; i < len; ++i) {
      if(nhinl[i]) {
          memcpy(uN+6*nhinlCount, u+6*i, 6*sizeof(double));
          memcpy(us+6*nhinlCount, u+6*i, 3*sizeof(double));
          /*Product of matrices :-)*/
          us[6*nhinlCount+3] = H[0]*u[6*i+3] + H[3]*u[6*i+4] + H[6]*u[6*i+5];
          us[6*nhinlCount+4] = H[1]*u[6*i+3] + H[4]*u[6*i+4] + H[7]*u[6*i+5];
          us[6*nhinlCount+5] = H[2]*u[6*i+3] + H[5]*u[6*i+4] + H[8]*u[6*i+5];
          ++nhinlCount;
        }
      if(hinl[i]) {
          memcpy(uH+6*hinlCount, u+6*i, 6*sizeof(double));
          ++hinlCount;
        }
    }
  //printf("nhinlCount = %u\n", nhinlCount);

  /*len = size(us,2);
        ptr = [1:len];
        max_i = 3;
        m_i = sam_sizO;
        max_sam = MAX_SAM;
        s_size = 2;*/
  ptr = (unsigned *) malloc(nhinlCount*sizeof(unsigned));
  for (i = 0; i < nhinlCount; ++i) {
      ptr[i] = i;
    }
  max_i = 3;
  m_i = sam_sizO;
  max_sam = MAX_SAM;
  s_size = 2;

  /*no_sam = 0;
        no_mod = 0;
        while no_sam < 2*max_sam
          for pos = 1:s_size
                  idx = pos + ceil(rand * (len-pos));
                  ptr([pos, idx]) = ptr([idx, pos]);
          end;

          no_sam = no_sam +1;*/
  //printf("Loop...\n");
  if (nhinlCount < 4 || hinlCount < 6) {
      max_i = 0;
    } else {
      for(no_sam = 1; no_sam < 2*max_sam; ++no_sam) {
          for (pos = 0; pos < s_size; ++pos) {
              idx = pos + 1 + rand()%(nhinlCount-pos-1);
              auxI = ptr[pos];
              ptr[pos] = ptr[idx];
              ptr[idx] = auxI;
            }

          /*ec = cross(cross(us(1:3,ptr(1)),(us(4:6,ptr(1)))), ...
                                                 cross(us(1:3,ptr(2)),(us(4:6,ptr(2)))));
                        aFt = skew_sym(ec / norm(ec));*/
          crossp(c1, us+6*ptr[0], us+6*ptr[0]+3);
          crossp(c2, us+6*ptr[1], us+6*ptr[1]+3);
          crossp(ec, c1, c2);
          ecNorm = sqrt(ec[0]*ec[0] + ec[1]*ec[1] + ec[2]*ec[2]);
          ec[0] = ec[0]/ecNorm;
          ec[1] = ec[1]/ecNorm;
          ec[2] = ec[2]/ecNorm;
          skew_sym(ec, aFt);

          /*Ds = FDs(aFt*H, uN);
                        v  = Ds < th*2;
                        no_i  = sum(v);*/

          mattr(Ht, H, 3, 3);
          mmul(aFtH, aFt, Ht, 3);
          mattr(aFt, aFtH, 3, 3);

          FDs(uN, aFt, Ds, nhinlCount);
          no_i = 0;
          for (i = 0; i < nhinlCount; ++i) {
              if (Ds[i] < th*2) {
                  ++no_i;
                  v[i] = 1;
                } else {
                  v[i] = 0;
                }
            }

          /*if no_i > m_i
                                m_i = no_i;
                                [aF,inl] = innerFH (uH, uN(:,v), u, th, 15, sam_sizH, sam_sizO);*/
          if (no_i > m_i) {
              //printf("no_i = %d\n", no_i);
              no_i = 0;
              for (i = 0; i < nhinlCount; ++i) {
                  if (v[i]) {
                      memcpy(uV+6*no_i, uN+6*i, 6*sizeof(double));
                      ++no_i;
                    }
                }
              m_i = no_i;
              innerFH(uH, hinlCount, uV, no_i, u, len, th, 15, sam_sizH, sam_sizO, aF, inl);

              /*if sum(inl) > max_i
                                          max_i = sum(inl);
                                          inls = inl;
                                          F = aF;
                                          maxni = sum(inl & nhinl);
                                          max_sam = min([max_sam,nsamples(maxni, len, 2, conf)]);
                                end
                        end
                end*/
              ninl = 0;
              for (i = 0; i < len; ++i) {
                  if (inl[i]) {
                      ++ninl;
                    }
                }
              if (ninl > max_i) {
                  max_i = ninl;
                  memcpy(F, aF, 3*3*sizeof(double));
                  maxni = 0;
                  for (i = 0; i < len; ++i) {
                      if(inl[i] && nhinl[i]) {
                          ++maxni;
                        }
                    }
                  max_sam = dmin(max_sam, nsamples(maxni, nhinlCount, 2, conf));
                }
            }
        }
    }

  /*fprintf(1,'in P+P %d samples\n',no_sam);
        fprintf(1,'   %d, %d\n', sum(hinl), max_i);*/
  //printf("in P+P %d samples\n",no_sam);
  //printf("   %d, %d\n", hinlCount, max_i);
  //printf("DeAllocations...\n");
  free(inl);
  free(Ds);
  free(nhinl);
  free(uN);
  free(us);
  free(uV);
  free(uH);
  free(ptr);
  free(v);
  return max_i;
}

unsigned dmin(unsigned a, unsigned b) {
  return a>b ? b : a;
}

/*%SampleCnt calculates number of samples needed to be done

function SampleCnt = nsamples(ni, ptNum, pf, conf)

if conf > 1
  error('Conf must be less then or equal to 1');
end

q  = prod ([(ni-pf+1) : ni] ./ [(ptNum-pf+1) : ptNum]);

if q < eps
   SampleCnt = Inf;
else  
   SampleCnt  = log(1 - conf) / log(1 - q);
end

if SampleCnt < 1
   SampleCnt = 1;
end */

/*unsigned nsamples(unsigned ni, unsigned np, unsigned ss, double conf) {
	unsigned i;
	double q = 1;
	if (conf >= 1) {
		fprintf(stderr, "Conf must be less than 1!\n");
		return (unsigned)-1;
	}
	for (i = 0; i < ss; ++i) {
		q *= (((double)(ni-i))/(np-i));
	}
	if (q < 1e-15) {
		return (unsigned)-1;
	}
	return ceil(log(1-conf) / log(1-q));
}*/


/*function [F,inl, tinl] = innerFH (uH, uO, u, th, num, sam_sizH, sam_sizO)*/
void innerFH(double * uH, unsigned lenH, double * uO, unsigned lenO,
	     double * u, unsigned len, double th, unsigned repCount, unsigned sam_sizH, unsigned sam_sizO,
	     double * F, unsigned char * inl) {

  unsigned i, rep, max_i, max_s, no_i;
  double aF[3*3];

  unsigned char * v = (unsigned char *) malloc(len*sizeof(unsigned char));
  double * usam = (double *) malloc(6*(sam_sizH+sam_sizO)*sizeof(double));
  double * Ds = (double *) malloc(len*sizeof(double));
  int * allInl = (int *) malloc((sam_sizH+sam_sizO)*sizeof(int));
  double * buffer = (double *) malloc(9*(sam_sizH+sam_sizO)*sizeof(double));

  for (i = 0; i < sam_sizH+sam_sizO; ++i) {
      allInl[i] = i;
    }

  /*lenH = size(uH,2);
        lenO = size(uO,2);

        ptrH = int32([1:lenH]-1);
        ptrO = int32([1:lenO]-1);

        sH = int32(sam_sizH);
        sO = int32(sam_sizO);

        F = ones(3);
        inl = zeros(1,size(u,2));
        tinl = inl;
        max_i = 0;
        max_s = 0;*/
  for (i = 0; i < 3*3; ++i) {
      F[i] = 1;
    }
  for (i = 0; i < len; ++i) {
      inl[i] = 0;
    }
  max_i = 0;
  max_s = 0;

  /*for rep = 1:num*/
  for (rep = 0; rep < repCount; ++rep) {

      /*usam = [rsample(ptrH, uH, sH),rsample(ptrO, uO, sO)];*/
      dual_sample(uH, lenH, sam_sizH, uO, lenO, sam_sizO, usam);

      /*aF = fu2F(usam);
                Ds = fFDs(aF,u);*/
      u2f(usam, allInl, sam_sizH+sam_sizO, aF, buffer);

      FDs(u, aF, Ds, len);

      /*v  = Ds < th;
                no_i  = sum(v);*/
      no_i = 0;
      for (i = 0; i < len; ++i) {
          if (Ds[i] < th) {
              v[i] = 1;
              ++no_i;
            } else {
              v[i] = 0;
            }
        }

      /*if max_i < no_i
                  inl = v;
                  F = aF;
                  max_i = no_i;
                end*/
      if (max_i < no_i) {
          memcpy(inl, v, len*sizeof(unsigned char));
          memcpy(F, aF, 3*3*sizeof(double));
          max_i = no_i;
        }

      /*if no_i > max_s
                        max_s = no_i;
                        [aF, v] = u2Fit(u,aF,th, th*3);
                        no_i  = sum(v);*/
      if (no_i > max_s) {
          max_s = no_i;
          no_i = u2Fit(u, len, aF, v, th, th*3, 4);

          /*if max_i < no_i
                                inl = v;
                                F = aF;
                                max_i = no_i;
                        end
                end*/
          if (max_i < no_i) {
              memcpy(inl, v, len*sizeof(unsigned char));
              memcpy(F, aF, 3*3*sizeof(double));
              max_i = no_i;
            }
        }

      /*tinl = tinl + inl;
        end*/
    }

  free(usam);
  free(Ds);
  free(allInl);
  free(buffer);
  free(v);
}

/*Draw combined sample from two arrays of TCs*/
void dual_sample(double * uA, unsigned lenA, unsigned sA, double * uB, unsigned lenB, unsigned sB, double * usam) {
  unsigned idx, pos, i;
  unsigned * ptrA = (unsigned *) malloc(lenA*sizeof(unsigned));
  unsigned * ptrB = (unsigned *) malloc(lenB*sizeof(unsigned));

  for (i = 0; i < lenA; ++i) {
      ptrA[i] = i;
    }
  for (i = 0; i < lenB; ++i) {
      ptrB[i] = i;
    }

  /*Shuffle pointers - at least the first 's̈́'*/
  for (pos = 0; pos < sA; ++pos) {
      idx = rand() % lenA;
      i = ptrA[pos];
      ptrA[pos] = ptrA[idx];
      ptrA[idx] = i;
    }
  for (pos = 0; pos < sB; ++pos) {
      idx = rand() % lenB;
      i = ptrB[pos];
      ptrB[pos] = ptrB[idx];
      ptrB[idx] = i;
    }

  /*Fill samples*/
  for (i = 0; i < sA; ++i) {
      memcpy(usam+6*i, uA+6*ptrA[i], 6*sizeof(double));
    }
  for (i = 0; i < sB; ++i) {
      memcpy(usam+6*(i+sA), uB+6*ptrB[i], 6*sizeof(double));
    }

  free(ptrA);
  free(ptrB);
}

/*function [F,inl] = u2Fit(u, F, th, ths)*/
unsigned u2Fit(double * u, unsigned len, double * F, unsigned char * inl, double th, double ths, unsigned iters) {
  /*Error threshold step*/
  double dth = (ths - th) / (iters - 1);
  unsigned iter, i, no_i;
  int * inlI = (int *) malloc(len*sizeof(int));
  double * Ds = (double *) malloc(len*sizeof(double));
  double * buffer = (double *) malloc(9*len*sizeof(double));

  for (iter = 0; iter < iters; ++iter) {
      /*New inliers from old F*/
      FDs(u, F, Ds, len);
      no_i = 0;
      for (i = 0; i < len; ++i) {
          if (Ds[i] < ths) {
              inl[i] = 1;
              ++no_i;
            } else {
              inl[i] = 0;
            }
        }
      /*Not enough inliers for F estimation*/
      if (no_i < 8) {
          free(inlI);
          free(Ds);
          free(buffer);
          return no_i;
        }
      /*New F from new inliers*/
      no_i = 0;
      for (i = 0; i < len; ++i) {
          if (inl[i]) {
              inlI[no_i++] = i;
            }
        }
      u2f(u, inlI, no_i, F, buffer);

      ths -= dth;
    }

  /*New inliers from new F*/
  FDs(u, F, Ds, len);
  no_i = 0;
  for (i = 0; i < len; ++i) {
      if (Ds[i] < th) {
          inl[i] = 1;
          ++no_i;
        } else {
          inl[i] = 0;
        }
    }

  free(inlI);
  free(Ds);
  free(buffer);
  return no_i;
}


unsigned innerH(double * H, double * u, unsigned len, double th, unsigned iters, unsigned char * inl, int * pool, double * buffer)
{
  double *err, *d, *Z;
  double *errs[5];
  int i, j, I, *inliers;
  Score S;

  err = (double *) malloc(len * 4 * sizeof(double));
  for (i=0; i<4; i++) {
      errs[i] = err + i * len;
    }
  inliers = (int *) malloc(sizeof(int) * len);
  for (i=0; i<len; i++) pool[i] = i;
  Z = (double *) malloc(len * 18 * sizeof(double));
  lin_hg(u, Z, pool, len);

  d = errs[0];

  HDs(Z, u, H, d, len);
  S = inlidxs(d, len, th, inliers);
  S = inHrani (u, len, inliers, S.I, th, Z, errs, buffer, H, iters);

  d = errs[0];
  I = 0;
  for (j = 0; j < len; j++) {
      if (d[j] <= th) {
          ++I;
          inl[j] = 1;
        } else {
          inl[j] = 0;
        }
    }

  free(err);
  free(Z);
  free(inliers);
  
  return I;
}

void transformInliers(int * inl, int * inl2, unsigned inlCount, unsigned len) {
  unsigned i;
  for (i = 0; i < len; ++i) {
      inl2[i] = 0;
    }
  for (i = 0; i < inlCount; ++i) {
      inl2[inl[i]] = 1;
    }
}



