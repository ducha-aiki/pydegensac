#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rtools.h"

#ifdef WIN32
#define random rand
#else
#include "bsd_random.h"
#define random degensac_random
#endif

/* inline int sample (int *pool, int max_sz, int i) */
int sample (int *pool, int max_sz, int i)
{
  int j,q,s;

  s = random() % (max_sz - i);
  j = max_sz - i - 1;
  q = pool[s];
  pool[s] = pool[j];
  pool[j] = q;

  return q;
}

int * randsubset (int * pool, int max_sz, int siz)
{
  int i,j,q,s;

  for (i = 0; i < siz; i++)
    {
      s = random() % (max_sz - i);
      j = max_sz - i - 1;
      q = pool[s];
      pool[s] = pool[j];
      pool[j] = q;
    }

  return pool + max_sz - siz;
}

void addcorrT (double *src, int dat_siz, int max_sz, double *dst)
{
  int j;
  
  for (j = 0; j < dat_siz; j++)
    {
      *dst = *src;
      dst ++;
      src += max_sz;
    }
}

void rsample (double *data, int dat_siz, 
              int *pool, int size, int max_sz, double *dst)
{
  double *src, *p;
  int i, j, q;
  
  for (i = 0; i < size; i++)
    {
      q = sample (pool, max_sz, i);
      src = data + q;
      p = dst + i;
      for (j = 0; j < dat_siz; j++)
        {
          *p = *src;
          p += size;
          src += max_sz;
        }
    }
}  


void rsampleT (double *data, int dat_siz, 
               int *pool, int size, int max_sz, double *dst)
{
  double *src, *p;
  int i, j, q;
  p = dst;
  
  for (i = 0; i < size; i++)
    {
      q = sample (pool, max_sz, i);
      src = data + q;
      for (j = 0; j < dat_siz; j++)
        {
          *p = *src;
          p ++;
          src += max_sz;
        }
    }
}

void rsampleTn (double *data, int dat_siz, int *pool, 
                int size, int n, int max_sz, double *dst)
{
  double *src, *p;
  int i, j, q;
  p = dst;
  
  for (i = 0; i < size; i++)
    {
      q = sample (pool, n, i);
      src = data + q;
      for (j = 0; j < dat_siz; j++)
        {
          *p = *src;
          p ++;
          src += max_sz;
        }
    }
}

void multirsample (double *data, int dat_siz, int dps, 
                   int *pool, int size, int max_sz, double *dst)
{
  double *src, *p;
  int i, j, k, q;
  int  ssh = dps * max_sz, dsh = dps * size;
  
  for (i = 0; i < size; i++)
    {
      q = sample (pool, max_sz, i);
      src = data + dps*q;
      p = dst + i * dps;
      for (j = 0; j < dat_siz; j++)
        {
          for (k = 0; k < dps; k++)
            p[k] = src[k];
          p += dsh;
          src += ssh;
        }
    }
}

void multirsampleT (double *data, int dat_siz, int dps, 
                    int *pool, int size, int max_sz, double *dst)
{
  double *src, *p;
  int i, j, k, q;
  int  ssh = dps * max_sz, bsz = dat_siz*dps;
  
  for (i = 0; i < size; i++)
    {
      q = sample (pool, max_sz, i);
      src = data + dps*q;
      p = dst + i *bsz;
      for (j = 0; j < dat_siz; j++)
        {
          for (k = 0; k < dps; k++)
            p[k*dat_siz] = src[k];
          p += 1;
          src += ssh;
        }
    }
}


/*Indexes of inliers with error lower than given threshold. Returns RANSAC score.*/
Score inlidxs (const double * err, int len, double th, int * inl) {
  int i;
  Score s = {0,0,0,0};
  /* The MSAC gain needs 1/(th*9/4) per correspondence; computing it once and
     multiplying is a numerics change (see truncQuadInv). th == 0 has no usable
     reciprocal and scores 0 for every point, so it keeps the scalar path --
     the flag is loop-invariant, so the branch hoists out. */
  const int score_j = (th != 0);
  const double inv = score_j ? 1/(th*9/4) : 0;
  /* Four partial sums rather than one: the gain accumulator was a serial FP
     dependency chain over every correspondence, running at add latency when
     the work is trivially parallel. Reassociation changes rounding.
     Branchless compress store: the write always happens, the index only
     advances for inliers. Every caller allocates `inl` at `len` entries, so
     the speculative write at s.I is in bounds. */
  double j0 = 0, j1 = 0, j2 = 0, j3 = 0;
  int n4 = len & ~3;
  if (score_j) {
      for (i = 0; i < n4; i += 4) {
          j0 += truncQuadInv(err[i],   inv);
          j1 += truncQuadInv(err[i+1], inv);
          j2 += truncQuadInv(err[i+2], inv);
          j3 += truncQuadInv(err[i+3], inv);
          inl[s.I] = i;     s.I += (err[i]   <= th);
          inl[s.I] = i+1;   s.I += (err[i+1] <= th);
          inl[s.I] = i+2;   s.I += (err[i+2] <= th);
          inl[s.I] = i+3;   s.I += (err[i+3] <= th);
        }
      for (i = n4; i < len; ++i) {
          j0 += truncQuadInv(err[i], inv);
          inl[s.I] = i;
          s.I += (err[i] <= th);
        }
      s.J = (j0 + j1) + (j2 + j3);
    } else {
      for (i = 0; i < len; ++i) {
          inl[s.I] = i;
          s.I += (err[i] <= th);
        }
    }
  return s;
}


int inlidxso (const double * err, const double * sgn, int len, double th,
              int * inl_buff, int ** inls)
{
  int i, po = 0, ne = 0;
  
  for(i = 0; i < len; i ++)
    if (err[i] <= th)
      {
        if (sgn[i] > 0)
          {
            inl_buff[po] = i;
            po ++;
          } else
          {
            ne ++;
            inl_buff[len-ne] = i;
          }
      }

  if (po >= ne)
    {
      *inls = inl_buff;
      return po;
    }
  *inls = inl_buff + len - ne;
  return ne;
}

int nsamples(int ninl, int ptNum, int samsiz, double conf)
{
  double a = 1, b = 1;
  int i;

  for (i = 0; i < samsiz; i++)
    {
      a *= ninl-i;
      b *= ptNum -i;
    }
  a = a/b;
  if (a < DEGENSAC_EPS)
    return MAX_SAMPLES;
  a = 1-a;
  if (a < DEGENSAC_EPS)
    return 1;
  else
    {
      b = log(1-conf) / log(a);
      if (b > MAX_SAMPLES)
        return MAX_SAMPLES; else
        return (int) ceil(b);
    }
}


int scoreLess(const Score s1, const Score s2) {
#if __SCORE__ == SC_M
  return s1.J < s2.J;
#endif

#if __SCORE__ == SC_H
  if (s1.I == s2.I) { /*Key feature of hybrid scoring - compare MLE score in case of equal inliers count*/
      return s1.J < s2.J;
    }
#endif
  return s1.I < s2.I;
}


void loadSample(double * u, int * samidx, unsigned sample_size, unsigned data_size, double * u_out) {
  unsigned i;
  for (i = 0; i < sample_size; ++i) {
      memcpy(u_out + data_size*i, u + data_size*samidx[i], data_size * sizeof(double));
    }
}
