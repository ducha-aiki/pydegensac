/*  qrlsq.c    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */
#include <stdlib.h>
#include <math.h>
int solvru(double *a,double *b,int n);
double qrlsq(double *a,double *b,int m,int n,int *f)
{ double *p,*q,*w;
  double s,h,r;
  int i,j,k,mm,ms;
  if(m<n) return -1;
  w=(double *)calloc(m,sizeof(double));
  for(i=0,mm=m,p=a; i<n ;++i,--mm,p+=n+1){
    if(mm>1){
      for(j=0,q=p,s=0.; j<mm ;++j,q+=n){
	w[j]= *q; s+= *q* *q;
       }
      if(s>0.){
	h=sqrt(s); if(*p<0.) h= -h;
	s+= *p*h; s=1./s; w[0]+=h;
	for(k=1,ms=n-i; k<ms ;++k){
	  for(j=0,q=p+k,r=0.; j<mm ;q+=n) r+=w[j++]* *q;
	  r=r*s;
	  for(j=0,q=p+k; j<mm ;q+=n) *q-=r*w[j++];
	 }
        *p= -h;
        for(j=0,q=b+i,r=0.; j<mm ;) r+=w[j++]* *q++;
        for(j=0,q=b+i,r*=s; j<mm ;) *q++ -=r*w[j++];
       }
     }
   }
  *f=solvru(a,b,n);
  for(j=n,q=b+j,s=0.; j<m ;++j,++q) s+= *q* *q;
  free(w);
  return s;
} 




