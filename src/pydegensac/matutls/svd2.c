#include "matutl.h"
#include <stdlib.h>
#include <stdio.h>
#include "svduv.c"
#include "ldvmat.c"
#include "ldumat.c"
#include "atou1.c"
#include "qrbdu1.c"
#include "qrbdv.c"
#include "svdu1v.c"

void pm(double *p, int c, int r)
{
  int i,j;

  for (i=0; i<c; i++)
  {
    for (j=0; j<r; j++, p++)
      printf("%4.2f ", *p);
    printf("\n");
  }
 printf("\n");
}

int main(int argc, char **argv)
{
 double ma[9] = {1,2,3,4,5,6,7,8,9};
 double ou[100];
 double d[9];
 double ov[100];
 int m, n;

 m = 3;
 n = 2;

 pm(ma,m,n);

 /* svduv(d,ma,ou,m,ov,n); */

  svdu1v(d,ma,m,ov,n);

  /*m(ou,m,m); */
 pm(d,1,5); 
 pm(ov,n,n); 

 return 0;
}
