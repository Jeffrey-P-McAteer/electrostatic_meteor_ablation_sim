

#include "eppic-pic.h"
#include "eppic-math.h"

/** This routine initializes a particle distribution so that its
    velocity distribution matches a Maxwellian distribution in 
    3D. 
    
*/

void maxwel_vel_dist_loading(particle_dist &pic, 
			     int vel_dim,
			     FTYPE vthx, FTYPE vthy, FTYPE vthz,
			     int common_seed) 
{

  int randseed = common_seed+10;

  if (vel_dim==1) {
    for (long long int ipart=0; ipart<pic.np; ipart++) {
      pic.vx(ipart) = vthx*gasdev(randseed);
    }
  } 
  if(vel_dim==2) {
    for (long long int ipart=0; ipart<pic.np; ipart++) {
      pic.vx(ipart) = vthx*gasdev(randseed);
      pic.vy(ipart) = vthy*gasdev(randseed);
    }
  } 
  if (vel_dim==3) {
    for (long long int ipart=0; ipart<pic.np; ipart++) {
      pic.vx(ipart) = vthx*gasdev(randseed);
      pic.vy(ipart) = vthy*gasdev(randseed);
      pic.vz(ipart) = vthz*gasdev(randseed);

    }
    
  }


}
