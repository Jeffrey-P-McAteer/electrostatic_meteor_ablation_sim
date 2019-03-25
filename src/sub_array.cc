/* returns the FArrayND corresponding to a given FArrayND_ranged without it's
   guard cells.
*/


#include<eppic.h>

void sub_array(FArrayND_ranged &in,FArrayND &out,
	       INDICIES(int startx, int starty, int startz),
	       INDICIES(int endbeforex, int endbeforey, int endbeforez)) 
{
  
  for (int ix=startx; ix<endbeforex; ix++) 
    {
#if NDIM >=2
      for (int iy=starty; iy<endbeforey; iy++) 
	{
#endif
#if NDIM == 3
	  for (int iz=startz; iz<endbeforez; iz++)
	    {
#endif
	      out(INDICIES(ix,iy,iz)) = in(INDICIES(ix,iy,iz));
#if NDIM >=2
	    }
#endif
#if NDIM == 3
	}
#endif
    }
  
}

void sub_array_add(FArrayND_ranged &in,FArrayND &out,
		   INDICIES(int startx, int starty, int startz),
		   INDICIES(int endbeforex, int endbeforey, int endbeforez)) 
{
  
  for (int ix=startx; ix<endbeforex; ix++) 
    {
#if NDIM >=2
      for (int iy=starty; iy<endbeforey; iy++) 
	{
#endif
#if NDIM == 3
	  for (int iz=startz; iz<endbeforez; iz++)
	    {
#endif
	      out(INDICIES(ix,iy,iz)) += in(INDICIES(ix,iy,iz));
#if NDIM >=2
	    }
#endif
#if NDIM == 3
	}
#endif
    }
  
}

void sub_array_mult(FArrayND_ranged &in,FArrayND &out,
		    INDICIES(int startx, int starty, int startz),
		    INDICIES(int endbeforex, int endbeforey, int endbeforez)) 
{
  
  for (int ix=startx; ix<endbeforex; ix++) 
    {
#if NDIM >=2
      for (int iy=starty; iy<endbeforey; iy++) 
	{
#endif
#if NDIM == 3
	  for (int iz=startz; iz<endbeforez; iz++)
	    {
#endif
	      out(INDICIES(ix,iy,iz)) *= in(INDICIES(ix,iy,iz));
#if NDIM >=2
	    }
#endif
#if NDIM == 3
	}
#endif
    }
  
}
