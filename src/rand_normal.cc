// This function is taken from pg 289 of Numerical Recipes.  
// It generates a random number from a Normal Distribution.
#include <math.h>
#include "Random.h"

float rand_normal(unsigned long long int iran)
{
  static RanFast rnd(iran);
  static int iset = 0;
  static float gset;
  float fac, rsq, v1, v2;
  
  if (iset == 0) {
    do {
      v1 = 2.0*rnd.dbl()-1.0;
      v2 = 2.0*rnd.dbl()-1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}
