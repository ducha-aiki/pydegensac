/*  ccmath.h    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */
/*
                               CCM

                Numerical Analysis Toolkit Header File
                      ELF Shared Library Version
*/
/* Required for Shared Library */

#define XMATH 1

/* Define File Pointers and Standard Library */

#include <stdio.h>
#include <stdlib.h>

/* Definitions of Types */

#ifndef NULL
#define NULL ((void *)0
#endif

/* Complex Types */

#ifndef CPX
struct complex
{
    double re,im;
};
typedef struct complex Cpx;
#define CPX  1
#endif

/* Orthogonal Polynomial Type */

#ifndef OPOL
struct orpol
{
    double cf,hs,df;
};
typedef struct orpol Opol;
#define OPOL 1
#endif

/* Tree Types */

#ifdef BAL
struct tnode
{
    char *key,*rec;
    int bal;
    struct tnode *pr,*pl;
};
#else
struct tnode
{
    char *key,*rec;
    struct tnode *pr,*pl;
};
#endif

/* Time Series Types */

struct mcof
{
    double cf;
    int lag;
};
struct fmod
{
    int fac;
    double val;
};

/* List Definition */

struct llst
{
    char *pls;
    struct llst *pt;
};

/* Hash Table Definition */

struct tabl
{
    char *key,*val;
    struct tabl *pt;
};

/* Extended Precision Types */

/* XMATH must be defined to use extended precision functions */
#ifdef XMATH
#ifndef XPR
#define XDIM 7
struct xpr
{
    unsigned short nmm[XDIM+1];
};
extern unsigned short m_sgn,m_exp;
extern short bias;
extern int itt_div,k_tanh;
extern int ms_exp,ms_trg,ms_hyp;
extern short max_p,k_lin;
extern short d_bias,d_max,d_lex;
extern struct xpr zero,one,two,ten;
extern struct xpr x_huge;

/* Variables used in extended precision arithmetic */

unsigned short m_sgn=0x8000,m_exp=0x7fff;
short bias=16383;
int itt_div=2,k_tanh=5;
int ms_exp=21,ms_hyp=25,ms_trg=31;
short max_p=16*XDIM,k_lin= -8*XDIM;
short d_bias=15360,d_max=2047,d_lex=12;
struct xpr zero= {{0x0,0x0}};
struct xpr one= {{0x3fff,0x8000}};
struct xpr two= {{0x4000,0x8000}};
struct xpr ten= {{0x4002,0xa000}};
struct xpr x_huge= {{0x7fff,0x0}};

/* Variables used in the extended precision math functions */

struct xpr pi4= {{0x3FFE,0xC90F,0xDAA2,0x2168,0xC234,0xC4C6,0x628B,0x80DC}};
struct xpr pi2= {{0x3FFF,0xC90F,0xDAA2,0x2168,0xC234,0xC4C6,0x628B,0x80DC}};
struct xpr pi= {{0x4000,0xC90F,0xDAA2,0x2168,0xC234,0xC4C6,0x628B,0x80DC}};
struct xpr ee= {{0x4000,0xADF8,0x5458,0xA2BB,0x4A9A,0xAFDC,0x5620,0x273D}};
struct xpr ln2= {{0x3FFE,0xB172,0x17F7,0xD1CF,0x79AB,0xC9E3,0xB398,0x3F3}};
struct xpr srt2= {{0x3FFF,0xB504,0xF333,0xF9DE,0x6484,0x597D,0x89B3,0x754B}};
#define XPR 1
#endif
#endif


/*     FUNCTION DECLARATIONS   */


/*   Linear Algebra     */


/* Real Linear Systems */


int minv(double *a,int n) ;

int psinv(double *v,int n) ;

int ruinv(double *a,int n) ;

int solv(double *a,double *b,int n) ;

int solvps(double *s,double *x,int n) ;

int solvru(double *a,double *b,int n) ;

void solvtd(double *a,double *b,double *c,double *x,int m) ;

void eigen(double *a,double *eval,int n) ;

void eigval(double *a,double *eval,int n) ;

double evmax(double *a,double *u,int n) ;

int svdval(double *d,double *a,int m,int n) ;

int sv2val(double *d,double *a,int m,int n) ;

int svduv(double *d,double *a,double *u,int m,double *v,int n) ;

int sv2uv(double *d,double *a,double *u,int m,double *v,int n) ;

int svdu1v(double *d,double *a,int m,double *v,int n) ;

int sv2u1v(double *d,double *a,int m,double *v,int n) ;

void mmul(double *mat,double *a,double *b,int n) ;

void rmmult(double *mat,double *a,double *b,int m,int k,int n) ;

void vmul(double *vp,double *mat,double *v,int n) ;

double vnrm(double *u,double *v,int n) ;

void matprt(double *a,int n,int m,char *fmt) ;

void fmatprt(FILE *fp,double *a,int n,int m,char *fmt) ;

void trnm(double *a,int n) ;

void mattr(double *a,double *b,int m,int n) ;

void otrma(double *at,double *u,double *a,int n) ;

void otrsm(double *st,double *u,double *s0,int n) ;

void mcopy(double *a,double *b,int m) ;

void ortho(double *evc,int n) ;

void smgen(double *a,double *eval,double *evec,int n) ;

/* utility routines for real symmertic eigensystems */

void house(double *a,double *d,double *ud,int n) ;

void housev(double *a,double *d,double *ud,int n) ;

int qreval(double *eval,double *ud,int n) ;

int qrevec(double *eval,double *evec,double *dp,int n) ;

/* utility routines for singular value decomposition */

int qrbdi(double *d, double *e,int n) ;

int qrbdv(double *d, double *e,double *u,int m,double *v,int n) ;

int qrbdu1(double *d, double *e,double *u,int m,double *v,int n) ;

void ldumat(double *a,double *u,int m,int n) ;

void ldvmat(double *a,double *v,int n) ;

void atou1(double *a,int m,int n) ;

void atovm(double *v,int n) ;


/* Complex Matrix Algebra */


int cminv(Cpx *a,int n) ;

int csolv(Cpx *a,Cpx *b,int n) ;

void heigvec(Cpx *a,double *eval,int n) ;

void heigval(Cpx *a,double *eval,int n) ;

double hevmax(Cpx *a,Cpx *u,int n) ;

void cmmul(Cpx *c,Cpx *a,Cpx *b,int n) ;

void cmmult(Cpx *c,Cpx *a,Cpx *b,int m,int k,int n) ;

void cvmul(Cpx *vp,Cpx *mat,Cpx *v,int n) ;

Cpx cvnrm(Cpx *u,Cpx *v,int n) ;

void cmprt(Cpx *a,int n,int m,char *fmt) ;

void trncm(Cpx *a,int n) ;

void hconj(Cpx *u,int n) ;

void cmattr(Cpx *a,Cpx *b,int m,int n) ;

void utrncm(Cpx *at,Cpx *u,Cpx *a,int n) ;

void utrnhm(Cpx *ht,Cpx *u,Cpx *h0,int n) ;

void cmcpy(Cpx *a,Cpx *b,int n) ;

void unitary(Cpx *u,int n) ;

void hmgen(Cpx *h,double *eval,Cpx *u,int n) ;


/* utility routines for hermitian eigen problems */

void chouse(Cpx *a,double *d,double *ud,int n) ;

void chousv(Cpx *a,double *d,double *ud,int n) ;

void qrecvc(double *eval,Cpx *evec,double *ud,int n) ;



/*   Geometry    */



void crossp(double *h,double *u,double *v) ;

double dotp(double *u,double *v,int m) ;

double metpr(double *u,double *a,double *v,int n) ;

void scalv(double *r,double s,int n) ;

void trvec(double *c,double *a,double *b,int n) ;

double leng(double *a,double *b,int n)  ;

void rotax(double *v,double az,double pa,double ang,int k) ;

void euler(double *pv,int m,double a,double b,double c) ;

/*    plane trigonometry   */

void trgsas(double a,double g,double b,double *ans);

int trgasa(double a,double ss,double b,double *asn);

double trgarea(double a,double b,double c);

int trgsss(double a,double b,double c,double *ang);

int trgssa(double a,double b,double ba,double *an);

/*    spherical trigonometry  */

void stgsas(double a,double g,double b,double *ang);

int stgasa(double a,double c,double b,double *ang);

int stgsss(double a,double b,double c,double *ang);

int stgaaa(double a,double b,double c,double *ang);

double stgarea(double a,double b,double c);

/*    hyperbolic trigonometry  */

void htgsas(double a,double g,double b,double *an);

int htgasa(double a,double cc,double b,double *ans);

int htgsss(double a,double b,double c,double *ang);

int htgaaa(double a,double b,double c,double *as);

double htgarea(double a,double b,double c);



/*   Numerical Integration    */



double fintg(double a,double b,int n,double te,double (*func)()) ;

/* functional form: double (*func)(double) */

double chintg(double *a,int m,double (*func)()) ;

/* functional form: double (*func)(double) */

double fchb(double x,double *a,int m) ;

int deqsy(double *y,int n,double a,double b,int nd,double te,
          int (*fsys)()) ;

/* functional form: int (*fsys)(double x,double *y,double *dp) */



/*   Optimization and Roots   */



int optmiz(double *x,int n,double (*func)(),double de,
           double test,int max) ;

/* functional form: double (*func)(double *x) */

double optsch(double (*func)(),double a,double b,double test) ;

/* functional form: double (*func)(double) */

int plrt(double *cof,int n,struct complex *root,double ra,double rb) ;

struct complex polyc(struct complex z,double *cof,int n) ;

double secrt(double (*func)(),double x,double dx,double test) ;

/* functional form: double (*func)(double) */

int solnl(double *x,double *f,double (*fvec[])(),int n,double test) ;

/* functional form: double (*fvec[])(double *x) */

int solnx(double *x,double *f,double (*fvec[])(),double *jm,
          int n,double test) ;

/* functional form: double (*fvec[])(double *x) */




/*   Curve Fitting and Least Squares   */



void chcof(double *c,int m,double (*func)()) ;

/* functional form: double (*func)(double) */

void chpade(double *c,double *a,int m,double *b,int n) ;

double ftch(double x,double *a,int m,double *b,int n) ;

void cspl(double *x,double *y,double *z,int m,double tn) ;

void csplp(double *x,double *y,double *z,int m,double tn) ;

double csfit(double w,double *x,double *y,double *z,int m) ;

double tnsfit(double w,double *x,double *y,double *z,
              int m,double tn) ;

double dcspl(double x,double *u,double *v,double *z,int m) ;


/* polynominal least squares functions use the Opol structure. */

void plsq(double *x,double *y,int n,Opol *c,double *ssq,int m) ;

double pplsq(double *x,double *y,int n,double *b,int m) ;

double evpsq(double x,Opol *c,int m) ;

double evpsqv(double x, Opol *c,int m,double *sig,double sqv) ;

void psqcf(double *pc,Opol *c,int m) ;

void psqvar(double *var,double s,Opol *c,int m) ;


/* QR transformation for linear least squares. */

double qrlsq(double *a,double *b,int m,int n,int *f) ;

double qrvar(double *v,int m,int n,double ssq) ;


/* singular value decomposition least squares. */

double lsqsv(double *x,int *pr,double *var,double *d,double *b,
             double *v,int m,int n,double th) ;

int svdlsq(double *d,double *a,double *b,int m,double *v,int n) ;

int sv2lsq(double *d,double *a,double *b,int m,double *v,int n) ;


/* utility called by svdlsq and sv2lsq. */

int qrbdbv(double *d,double *e,double *b,double *v,int n) ;


/* nonlinear least squares */

double seqlsq(double *x,double *y,int n,double *par,double *var,
              int m,double de,double (*func)(),int kf) ;

/* functional form: double (*func)(double x,double *par) */

double gnlsq(double *x,double *y,int n,double *par,
             double *var,int m,double de,double (*func)()) ;

/* functional form: double (*func)(double x,double *par) */

double fitval(double x,double *s,double *par,double (*fun)(),
              double *v,int n) ;

/* functional form: double (*func)(double x,double *par) */

void setfval(int i,int n) ;



/*    Fourier Analysis    */



void fft2(struct complex *ft,int m,int inv) ;

void fft2_d(struct complex *a,int m,int n,int f) ;

void fftgc(struct complex **pc,struct complex *ft,int n,
           int *kk,int inv) ;

void fftgr(double *x,struct complex *ft,int n,int *kk,int inv) ;

void ftuns(struct complex **pt,int n) ;

int pfac(int n,int *kk,int fe) ;

void pshuf(Cpx **pa,Cpx **pb,int *kk,int n) ;

int pwspec(double *x,int n,int m) ;

void smoo(double *x,int n,int m) ;



/*   Simulation Support    */


double *autcor(double *x,int n,int lag) ;

int *hist(double *x,int n,double xmin,double xmax,
          int kbin,double *bin) ;

unsigned int lran1() ;

void setlran1(unsigned int seed) ;

unsigned int lrand() ;

void setlrand(unsigned int seed) ;

int bran(int n) ;

void setbran(unsigned int seed) ;

int bran2(int n) ;

void setbran2(unsigned int seed) ;

double unfl() ;

void setunfl(unsigned int seed) ;

double unfl2() ;

void setunfl2(unsigned int seed) ;

double nrml() ;

void setnrml(unsigned int seed) ;

void norm(double *err) ;

void setnorm(unsigned int seed) ;

void norm2(double *err) ;

void setnorm2(unsigned int seed) ;

void sampl(void **s,int n,void **d,int m) ;

void shuffl(void **s,int n) ;

/*      utility routines used for 2^31 - 1 modular arithmetic   */

unsigned int lrana(unsigned int s) ;

unsigned int lranb(unsigned int s) ;



/*   Sorts and Searches      */



int batdel(char *kin,struct tnode *hd) ;

struct tnode *batins(char *kin,struct tnode *hd) ;

struct tnode *btsearch(char *kin,struct tnode *hd) ;

void btsort(struct tnode *hd,struct tnode **ar) ;

void prbtree(struct tnode *hd,int m) ;

int btdel(char *kin,struct tnode *hd) ;

struct tnode *btins(char *kin,struct tnode *hd) ;

struct tnode *tsearch(char *kin,struct tnode *hd) ;

void tsort(struct tnode *hd,struct tnode **ar) ;

void prtree(struct tnode *hd,int m) ;

int hashdel(char *kin,struct tabl *harr[],int mh) ;

struct tabl *hashins(char *kin,struct tabl *harr[],int mh) ;

struct tabl *hfind(char *kin,struct tabl *harr[],int mh) ;

int hval(char *key,int mh) ;

struct llst *msort(struct llst *st,int dim,int (*comp)()) ;

void qsrt(void *v,int i,int j,int (*comp)()) ;

void hsort(void *v,int n,int (*comp)()) ;

void ssort(void *v,int n,int (*comp)()) ;

/* comparison functions for sort routines. */

/* define the functional form of int (*comp)() */

int dubcmp(double *x,double *y) ;

int intcmp(int *x,int *y) ;

int unicmp(unsigned *x,unsigned *y) ;

/* the standard library function strcmp will also work
   with these sorts */



/*   Statistical Distributions    */



double qnorm(double x) ;

double pctn(double pc) ;

double qgama(double x,double a) ;

double pctg(double pc,double a) ;

double qbeta(double x,double a,double b) ;

double pctb(double pc,double a,double b) ;

double qgnc(double x,double a,double d) ;

double pctgn(double pc,double a,double d) ;

double qbnc(double x,double a,double b,double d) ;

double pctbn(double pc,double a,double b,double d) ;



/*    Special Functions    */


/* elliptic integrals and functions */

double nome(double k,double *pk,double *pkp) ;

double amelp(double u,double k) ;

double theta(double u,int n) ;

void stheta(double k) ;

double felp(double an,double k,double *pk,double *pz,double *ph) ;

double gelp(double an,double k,double as,double bs,
            double ds,double *pg,double *pf,double *pk) ;

double g2elp(double an,double bn,double k,double as,
             double bs,double ds) ;


/* bessel functions */

double jbes(double v,double x) ;

double ibes(double v,double x) ;

double kbes(double v,double x) ;

double nbes(double v,double x) ;

double drbes(double x,double v,int f,double *p) ;

double rcbes() ;

void setrcb(double u,double y,int fl,int dr,double *pf,
            double *ph) ;


/* spherical bessel functions */

double jspbes(int n,double x) ;

double kspbes(int n,double x) ;

double yspbes(int n,double x) ;

double drspbes(double x,int n,int f,double *p) ;

double rcspbs() ;

void setrcsb(int n,double y,int fl,int dr,double *pf,double *ph) ;

/* airy functions */

double airy(double x,int df) ;

double biry(double x,int df) ;

/* gamma and related functions */

double gaml(double x) ;

double psi(int m) ;

double psih(double v) ;


/* support routines for evaluation of elliptic integrals */

double gsng(double *pa,double *pb,double *pc,double b,double an) ;

double gsng2(double *pa,double *pb,double *pc,double b,
             double an,double bn) ;



/*    Complex Arithmetic    */



struct complex cmul(struct complex s,struct complex t) ;

struct complex cdiv(struct complex s,struct complex t) ;

struct complex cadd(struct complex s,struct complex t) ;

struct complex csub(struct complex s,struct complex t) ;

struct complex crmu(double a,struct complex z) ;

struct complex cimu(double b,struct complex z) ;

struct complex ccng(struct complex z) ;

struct complex cdef(double r,double i) ;

double cabs(struct complex c) ;

double cnrm(struct complex z) ;

struct complex cexp(struct complex z) ;

struct complex clog(struct complex z) ;

struct complex csinh(struct complex z) ;

struct complex ccosh(struct complex z) ;

struct complex ctanh(struct complex z) ;

struct complex casinh(struct complex z) ;

struct complex cacosh(struct complex z) ;

struct complex catanh(struct complex z) ;

struct complex casin(struct complex z) ;

struct complex cacos(struct complex z) ;

struct complex catan(struct complex z) ;

struct complex csqrt(struct complex z) ;

struct complex csin(struct complex z) ;

struct complex ccos(struct complex z) ;

struct complex ctan(struct complex z) ;



/*    Time Series          */



double sarma(double er) ;

void setsim(int k) ;

double parma(double *x,double *e) ;

double evfmod(struct fmod y) ;

void setevf(int k) ;

double drfmod(struct fmod y,double *dr) ;

void setdrf(int k) ;

double seqtsf(struct fmod *x,int n,double *var,int kf) ;

double fixtsf(struct fmod *x,int n,double *var,double *cr) ;

double evmod(double y) ;

void setev(int k) ;

double drmod(double y,double *dr) ;

void setdr(int k) ;

double seqts(double *x,int n,double *var,int kf) ;

double fixts(double *x,int n,double *var,double *cr) ;

int resid(double *x,int n,int lag,double **pau,int nbin,
          double xa,double xb,int **phs,int *cks) ;

int sany(double *x,int n,double *pm,double *cd,double *ci,
         int nd,int ms,int lag) ;

double sdiff(double y,int nd,int k) ;

double sintg(double y,int nd,int k) ;

double xmean(double *x,int n) ;



/*    Extended Precision Arithmetic   */

/* XMATH must be defined to use these functions */

#ifdef XMATH

struct xpr xadd(struct xpr s,struct xpr t,int f) ;

struct xpr xmul(struct xpr s,struct xpr t) ;

struct xpr xdiv(struct xpr s,struct xpr t) ;

double xtodub(struct xpr s) ;

struct xpr dubtox(double y) ;

struct xpr inttox(int n) ;

int xprcmp(struct xpr *pa,struct xpr *pb) ;

struct xpr xneg(struct xpr s) ;

struct xpr xabs(struct xpr s) ;

int xex(struct xpr *ps) ;

int neg(struct xpr *ps) ;

struct xpr xfrex(struct xpr s,int *p) ;

struct xpr xfmod(struct xpr s,struct xpr t,int *p) ;

struct xpr xsqrt(struct xpr z) ;

struct xpr xexp(struct xpr z) ;

struct xpr xlog(struct xpr z) ;

struct xpr xpwr(struct xpr s,int n) ;

struct xpr xpr2(struct xpr s,int m) ;

struct xpr xtan(struct xpr z) ;

struct xpr xcos(struct xpr z) ;

struct xpr xsin(struct xpr z) ;

struct xpr xatan(struct xpr z) ;

struct xpr xasin(struct xpr z) ;

struct xpr xacos(struct xpr z) ;

struct xpr xtanh(struct xpr z) ;

struct xpr xsinh(struct xpr z) ;

struct xpr xcosh(struct xpr z) ;

struct xpr  atox(char *s) ;

void prxpr(struct xpr u,int lim) ;

void xprint(struct xpr x) ;

/* special applications */

void xchcof(struct xpr *cf,int m,struct xpr (*xfunc)()) ;

/* functional form: xpr (*xfunc)(xpr *cf) */

struct xpr xevtch(struct xpr z,struct xpr *a,int m) ;


/* utility operations on extended precision numbers */

struct xpr sfmod(struct xpr s,int *p) ;

void lshift(int n,unsigned short *pm,int m) ;

void rshift(int n,unsigned short *pm,int m) ;

#endif


/*   Utility Operations (on Bits)  */


unsigned short bset(unsigned short x,unsigned short n) ;

int bget(unsigned short x,unsigned short n) ;

int bcnt(unsigned short x) ;

unsigned int lbset(unsigned int x,int n) ;

int lbget(unsigned int x,int n) ;

int lbcnt(unsigned int x) ;

void bitpc(unsigned char x) ;

void bitps(unsigned short x) ;

void bitpl(unsigned int x) ;

void bitpf(float x);

void bitpd(double x) ;

#ifdef XMATH
void bpatx(struct xpr x) ;
#endif

double pwr(double y,int n) ;


/*
     special declarations required for shared library
*/

int np,nma,nar,nfc,ndif;
struct mcof *par,*pma,*pfc;
