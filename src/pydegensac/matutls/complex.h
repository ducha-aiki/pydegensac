/*  complex.h    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */
#ifndef CPX
struct complex
{
    double re,im;
};
typedef struct complex Cpx;
#define CPX  1
#endif

#include <math.h>
struct complex cadd(Cpx a,Cpx b);
struct complex csub(Cpx a,Cpx b);
struct complex cmul(Cpx a,Cpx b);
struct complex cdiv(Cpx a,Cpx b);
struct complex crmu(double x,Cpx a);
struct complex cimu(double y,Cpx a);
struct complex ccng(Cpx c);
struct complex cdef(double r,double i);
#ifndef WIN32
double cabs(Cpx a);
#endif
struct complex cnrm(Cpx a);
struct complex csqrt(Cpx a);
struct complex cexp(Cpx a);
struct complex clog(Cpx a);
struct complex csin(Cpx a);
struct complex ccos(Cpx a);
struct complex ctan(Cpx a);
struct complex casin(Cpx f);
struct complex cacos(Cpx f);
struct complex catan(Cpx f);
struct complex csinh(Cpx h);
struct complex ccosh(Cpx h);
struct complex ctanh(Cpx h);
struct complex casinh(Cpx g);
struct complex cacosh(Cpx g);
struct complex catanh(Cpx g);
