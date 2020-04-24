typedef void (*HDsPtr) (const double*, const double *,const double *, double *, int);
typedef void (*HDsiPtr) (const double*, const double *,const double *, double *, int, int*, int);
typedef void (*HDsidxPtr) (const double*, const double *,const double *, double *, int, int*, int);

#ifdef __cplusplus
extern "C"
#endif
void lin_hg(const double *u, double *dst, const int* inl, int len);

void lin_hgN(const double *u, double *p, const int* inl, int len,
             double *A1, double *A2);

void u2h(const double *u, const int *inl, int len, double *H, double * buffer);

void pinvJ (double a, double b, double c, double d, double e, double *pJ);



#ifdef __cplusplus
extern "C"
#endif
void HDs(const double *lin, const double * u,
         const double *H, double *p, int len);

#ifdef __cplusplus
extern "C"
#endif
void HDsSymSumSq(const double *lin, const double * u,
            const double *H, double *p, int len);

#ifdef __cplusplus
extern "C"
#endif
void HDsSymSum(const double *lin, const double * u,
            const double *H, double *p, int len);

#ifdef __cplusplus
extern "C"
#endif
void HDsSymMax(const double *lin, const double * u,
            const double *H, double *p, int len);


#ifdef __cplusplus
extern "C"
#endif
void HDsSymMaxSq(const double *lin, const double * u,
            const double *H, double *p, int len);


/* Sampson error for homography and point correspondences, computed only on a subset pts */
#ifdef __cplusplus
extern "C"
#endif
void HDsi(const double *lin, const double * u6,
          const double *H, double *p, int len, int *pts, int ni);
#ifdef __cplusplus
extern "C"
#endif
void HDsidx(const double *lin, const double * mu,
            const double *H, double *p, int len, int *idx, int siz);

//!!!! The diffence between HDsi and HDsidx is that HDsi stores result into p[0..size] (so, just beginning of the array)
//! while the HDsidx stores result into p[idx[0..size]] - stores in "correct" indexes

#ifdef __cplusplus
extern "C"
#endif
void HDsSymSumidx(const double *lin, const double * mu, const double *H,
               double *p, int len, int *idx, int siz);
#ifdef __cplusplus
extern "C"
#endif
void HDsSymSumSqidx(const double *lin, const double * mu, const double *H,
               double *p, int len, int *idx, int siz);

#ifdef __cplusplus
extern "C"
#endif
void HDsSymMaxidx(const double *lin, const double * mu, const double *H,
               double *p, int len, int *idx, int siz);
#ifdef __cplusplus
extern "C"
#endif
void HDsSymMaxSqidx(const double *lin, const double * mu, const double *H,
               double *p, int len, int *idx, int siz);

/**/
#ifdef __cplusplus
extern "C"
#endif
void HDsiSymSum(const double *lin, const double * u6,
          const double *H, double *p, int len, int *pts, int ni);

#ifdef __cplusplus
extern "C"
#endif
void HDsiSymSumSq(const double *lin, const double * u6,
          const double *H, double *p, int len, int *pts, int ni);


#ifdef __cplusplus
extern "C"
#endif
void HDsiSymMax(const double *lin, const double * u6,
          const double *H, double *p, int len, int *pts, int ni);

#ifdef __cplusplus
extern "C"
#endif
void HDsiSymMaxSq(const double *lin, const double * u6,
          const double *H, double *p, int len, int *pts, int ni);


int all_Hori_valid (double * us, int *idx);
void invertH(const double *H, double *Hinv);
int all_not_collinear(double * us, int *idx);
//int HcloseToSingular(const double *h);

