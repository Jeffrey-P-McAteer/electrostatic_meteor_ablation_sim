/* A routine to remove particles that have been placed in the nabsent array, 
   indicating they should be marked as outside the simulation. 
   This applies for open boundary conditions (not periodic). 
   
   Meers Oppenheim;   July, 2016 */    

#include "eppic.h"

void remove_particles(PTYPEAVec &pos, particle_dist  &adv, int vel_dim, int nabsent_start, 
		      bool left_edge_out, bool right_edge_out)
{
  //int nremoved = 0;
  for (int i=nabsent_start;i<adv.nabsent;i++) {
    int ip = adv.absent(i);
    PTYPE p=pos(ip);
    // mark as outside if conditions met:
    if ( (left_edge_out  && p < (PTYPE) 0.) ||  
	 (right_edge_out && p >= 0) ){  
      
#ifdef DEBUG
      // Test if the particle is already absent...
      if (p<=fabsent2) {
	terminate(-1,"remove_particles: particle being removed is already absent");
      }
#endif
      adv.x(ip)=fabsent2;
      if (NDIM >= 2) adv.y(ip)=fabsent2;
      if (NDIM >= 3) adv.z(ip)=fabsent2;
      // Insure that it does not move (also removes the energy)
      adv.vx(ip)=0.0;
      if (vel_dim >= 2) adv.vy(ip)=0.0;
      if (vel_dim == 3) adv.vz(ip)=0.0;
      adv.np_all--;
    }
  }

  //if (nremoved > 0){
  //printf("proc %d time %d: Removed %d particles\n", mpi_rank, it, nremoved);
  //}
}

