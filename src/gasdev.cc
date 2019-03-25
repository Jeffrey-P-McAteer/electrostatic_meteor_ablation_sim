#include <math.h>
#include "eppic-types.h"
#include "eppic.h"
FTYPE gasdev(int &idum)
{

  /*     Generate a true random gaussian  */
  FTYPE ran3(int *);
  static int iset=0;
  static FTYPE gset;
  FTYPE fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*ran3(&idum)-1.0;
      v2=2.0*ran3(&idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;

    return v2*fac;
  } else {
    iset=0;

    return gset;
  }
}

