/* take a ranged array and sum from 0 to nx and same for y,z

 */

#include<eppic.h>

FTYPE sum_noguard(FArrayND_ranged &in) 
{

  FTYPE sum = 0;

  for (int ix=0; ix < nx; ix++) {
#if NDIM >= 2    
    for (int iy=0; iy<ny; iy++) {
#endif      
#if NDIM == 3	
      for (int iz=0; iz<nz; iz++) {
#endif
	sum += in(INDICIES(ix,iy,iz));
#if NDIM == 3
	  }
#endif
#if NDIM >= 2
    }
#endif
  }
  return sum;
  
}
