/* Routine to modify the particle velocities for special purpose simulations
   This routine changes Vx as a function of Vy to introduce velocity shears 
   into the simulation */

#include <stdio.h>
#include <math.h>
#include "eppic.h"

void velocity_mod(particle_dist *&pic)
{
  
  // Loop through all the distributions
  for (int id=0; id<ndist; ++id) {
    // User sets param2 to the maximum Vx value to adjust to
    if (param2[id] != 0) {
      FTYPE omega=pic[id].q/pic[id].m*Bx; // Note: carries sign
      FTYPE avg_radius=fabs(sqrt((Sqr(vythd[id])+Sqr(vzthd[id]))/2.)/omega);
      // Set 3*avg_radius as the width of each shear region and 
      // set the remaining part of the simulation to be 
      // 1/2 at the min current and 1/2 at the max current
      /*      FTYPE shear_dy=3*avg_radius;
	      FTYPE flat=(ny*dy-2.*shear_dy);*/

      // Set the ramp-up and ramp-down regions to be 1/6 of the box.
      FTYPE shear_dy=ny*dy/6.;

      FTYPE flat=(ny*dy-2.*shear_dy);
      FTYPE y_start=flat/4.;
      FTYPE y_stop=y_start+flat/2.+2*shear_dy;

      if (flat < shear_dy)
	terminate(1,"Flat region less than shear region in velocity_mod");

      if (mpi_rank == 0) 
	printf("Distribution %d \tShear start: %g Shear stop: %g Shear size: %g\n",
	       id,y_start,y_stop,shear_dy);
      
      // Loop through particles of this distribution
      for (int i = 0; i < pic[id].np; ++i) {

	//Find the guiding center of the magnetized particle
	FTYPE vperp=0;
	if (species_dim[id] != 1) vperp=sqrt(Sqr(pic[id].vy[i]*(dy/dt))+
					     Sqr(pic[id].vz[i]*(dz/dt)));

	// FTYPE phase=tan2(pic[id].vy[i]*dy/dt,pic[id].vz[i]*dz*dt)
	FTYPE sin_phase=0;
	//FTYPE cos_phase=0;
	if (vperp >0 ) {
	  sin_phase=-pic[id].vz[i]*dz/dt/vperp;
	  // cos_phase=+pic[id].vy[i]*dy/dt/vperp;
	} else {
	  sin_phase=0;
	  // cos_phase=0;
	}
	FTYPE radius=vperp/omega; // Carries sign
	FTYPE yc=fmod(pic[id].y[i]*dy - radius*sin_phase + ny*dy, ny*dy);
	// FTYPE zc=fmod(pic[id].z[i]*dz - radius*cos_phase + nz*dz, nz*dz);
	
	//Use the guiding center to adjust Vx 
	// Four regions:
	if (yc <= y_start || yc >= y_stop) {
	  // Outside shearing regions subtract DC component.
	  pic[id].vx[i] -= (param2[id]/2.)*dt/dx;
	} else if (yc >= y_start && yc <= y_start+shear_dy){
	  // In left sheared region
	  FTYPE m=param2[id]/shear_dy;
	  FTYPE b=-param2[id]/2. - m*y_start;
	  pic[id].vx[i] += (m*yc + b)*dt/dx;
	} else if (yc >= y_start+shear_dy && yc <= y_stop-shear_dy){
	  // On the crest (which is flat and unsheared)
	  pic[id].vx[i] += (param2[id]/2.)*dt/dx;
	} else if (yc >= y_stop-shear_dy && yc <= y_stop){
	  // In right sheared region
	  FTYPE m=-param2[id]/shear_dy;
	  FTYPE b=-param2[id]/2. - m*y_stop;
	  pic[id].vx[i] += (m*yc + b)*dt/dx;
	} else {
	  // If we miss any particles, we may catch them here:
	  terminate(1,"Particle out of range in velocity_mod\n");
	}
	
      }
    }
  }
  
}
