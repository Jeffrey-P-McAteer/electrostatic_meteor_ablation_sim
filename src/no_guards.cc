/* returns the FArrayND corresponding to a given FArrayND_ranged without it's
   guard cells.
*/


#include<eppic.h>

FArrayND no_guards(FArrayND_ranged &in) 
{
  
  FArrayND sub_array = FArrayND(INDICIES(nx,ny,nz));
  for (int ix=0; ix<nx; ix++) 
    {
#if NDIM >=2
      for (int iy=0; iy<ny; iy++) 
	{
#endif
#if NDIM == 3
	  for (int iz=0; iz<nz; iz++)
	    {
#endif
	      sub_array(INDICIES(ix,iy,iz)) = in(INDICIES(ix,iy,iz));

#if NDIM >=2
	    }
#endif
#if NDIM == 3
	}
#endif
    }
  
  return sub_array;
}
