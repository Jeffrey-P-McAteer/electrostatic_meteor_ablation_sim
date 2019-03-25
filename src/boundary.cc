/* Handle boundary condition issues */

#include "eppic.h"

void boundary(particle_dist *pic, field &Efield, particle_misc *misc, 
	      int it) 
{
  static int todo=1;
  
  if (todo==0) {
    return;
  }
  
  /* Check to see if the charge_on boundary condition is still necessary 
     every eighth timestep*/
  if (it%8 ==0) {
    int all_done=TRUE;
    for (int id=0; id<ndist; ++id) 
      if (chargeon(id)==1) {
	int finished=TRUE;
	for (int i = 0; i < pic[id].np; ++i) 
	  if (misc[id].charge(i) == 0) finished=FALSE;
	if (finished == TRUE) {
	  chargeon(id)=0;
	  misc[id].charge.~PTYPEAVec(); /* Clear the charge array */
	} else all_done=FALSE;
      }
    if (all_done==TRUE) todo=0;
  }
}

      
	


					 
	    

	
