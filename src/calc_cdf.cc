
/** Calculates the Cumulative Distribution Function (CDF) for a given
    probability distribution (dist).

**/

#include "eppic-types.h"


void calc_cdf(FTYPEAVec &dist, FTYPEAVec &cdf) {

  int nx=dist.size();
  if (cdf.size()!=nx)
    cdf = FTYPEAVec(nx);
  
  cdf(0)=0;
  for(int ix=1;ix<nx;ix++) 
    cdf(ix)=cdf(ix-1)+dist(ix);
  for(int ix=1;ix<nx;ix++) 
    cdf(ix)/=cdf(nx-1);
  

}
