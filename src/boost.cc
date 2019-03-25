/* Routine to take a set of particles from a distribution and boost them to a
   higher energy. This mimics photo-ionization and annihilation in a steady
   state process such that for every particle destroyed, one would be created.
   We could also randomize the positions. */ 

#include "eppic.h"
#include "eppic-math.h"

void boost(int id, particle_dist &adv)
{
// Based on the boost_rate (particles/m^(-3)/s, decide how many particles will receive a boost: 
  static const PTYPE fnx = (PTYPE) nx;

  const unsigned long long int iran = 10430+mpi_rank;
  static RanFast rnd(iran);
  PTYPE fnboost=boost_rate[id]*(nx*dx*ny*dy*nz*dz)*dt;
  int nboost=static_cast<int>(fnboost);
  // The fractional particle will contribute:
  if ( (fnboost-nboost) > (rnd.dbl()) ) nboost += 1;

// pick these particles at random from all the particles to boost:
    for (int ib=0; ib<nboost; ib++) {

// Change the velocity of this particle into that of an shell one

      void create_shell_particle(PTYPE vradius, PTYPE vth, 
				 PTYPE &vx, PTYPE &vy, PTYPE &vz);
      PTYPE vx,vy,vz;
      create_shell_particle(v0_boost[id], vth_boost[id], vx, vy, vz);

      //Choose which particle to boost at random
      long i;
      do { i=static_cast<long>(rnd.dbl()*adv.np); } 
      while (adv.x(i)<-fnx);
      adv.vx(i)=vx*dt/dx;
      if (vel_dim[id] >= 2) adv.vy(i)=vy*dt/dy;
      if (vel_dim[id] >= 3) adv.vz(i)=vz*dt/dz;
	
    }
}
