/*  svduv.c    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */
#include <stdlib.h>
#include <math.h>
#include "matutl.h"
static int svduv_impl(double *d,double *a,double *u,int m,double *v,int n,
                      int null_only)
{ double *p,*p1,*q,*pp,*w,*e;
  double s,h,r,t,sv;
  int i,j,k,mm,nm,ms;
  if(m<n) return -1;
  w=(double *)calloc(m+n,sizeof(double)); e=w+m;
  for(i=0,mm=m,nm=n-1,p=a; i<n ;++i,--mm,--nm,p+=n+1){
    if(mm>1){ sv=h=0.;
      for(j=0,q=p,s=0.; j<mm ;++j,q+=n){
	w[j]= *q; s+= *q* *q;
       }
      if(s>0.){
	h=sqrt(s); if(*p<0.) h= -h;
	s+= *p*h; s=1./s; t=1./(w[0]+=h);
        sv=1.+fabs(*p/h);
	for(k=1,ms=n-i; k<ms ;++k){
	  for(j=0,q=p+k,r=0.; j<mm ;q+=n) r+=w[j++]* *q;
	  r*=s;
	  for(j=0,q=p+k; j<mm ;q+=n) *q-=r*w[j++];
	 }
        for(j=1,q=p; j<mm ;) *(q+=n)=t*w[j++];
       }
      *p=sv; d[i]= -h;
     }
    if(mm==1) d[i]= *p;
    p1=p+1; sv=h=0.;
    if(nm>1){
      for(j=0,q=p1,s=0.; j<nm ;++j,++q) s+= *q* *q;
      if(s>0.){
	h=sqrt(s); if(*p1<0.) h= -h;
        sv=1.+fabs(*p1/h);
	s+= *p1*h; s=1./s; t=1./(*p1+=h);
	for(k=n,ms=n*(m-i); k<ms ;k+=n){
	  for(j=0,q=p1,pp=p1+k,r=0.; j<nm ;++j) r+= *q++ * *pp++;
	  r*=s;
	  for(j=0,q=p1,pp=p1+k; j<nm ;++j) *pp++ -=r* *q++;
	 }
        for(j=1,q=p1+1; j<nm ;++j) *q++ *=t;
       }
      *p1=sv; e[i]= -h;
     }
    if(nm==1) e[i]= *p1;
   }
  /* For a 9x8 input, column 8 of u is already the null vector after
     ldumat. qrbdv rotates only columns 0..7, and the sign repair below
     changes v, not u. Keep the reduction arithmetic shared with svduv
     so the eight-point F solver gets exactly the same null vector. */
  if(null_only){
    ldumat(a,u,m,n);
    free(w);
    return 0;
   }
  ldvmat(a,v,n); ldumat(a,u,m,n);
  qrbdv(d,e,u,m,v,n);
  for(i=0; i<n ;++i){
    if(d[i]<0.){ d[i]= -d[i];
      for(j=0,p=v+i; j<n ;++j,p+=n) *p= - *p;
     }
   }
  free(w);
  return 0;
}

int svduv(double *d,double *a,double *u,int m,double *v,int n)
{
  return svduv_impl(d,a,u,m,v,n,0);
}

int svdnull8(double *a,double *u)
{
  double d[8];
  return svduv_impl(d,a,u,9,NULL,8,1);
}
