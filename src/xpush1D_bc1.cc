// Push the particle positions to the next timestep in 1D only: 

#include <csignal>
#include <math.h>
//#include <signal.h>
#include "eppic.h"

void xpush1D_bc1(particle_dist &adv, particle_dist &old, particle_dist &cur, 
		 particle_misc &misc, FTYPE alpha, FTYPE vx0, FTYPE vy0)
{
  
  int i;
  int np=old.x.length();
  PTYPE x;

#ifndef __GNUC__
  const PTYPE min_resolvable_value = 0.0;
#else
  // Rounding error problems in some compilers require using this value 
  const int other_bits = 12;
  const PTYPE min_resolvable_value = 1.0/pow(2.,8.*sizeof(PTYPE)-other_bits);
#endif

  PTYPE fnx = (PTYPE) nx*(1-min_resolvable_value);
    
  // For the sake of speed, distinguish between when alpha == 1 or not: 
  if (alpha != 1.) {
    PTYPE falpha = (PTYPE) alpha;
    
    // Update with the adv velocity: 
    for (i = 0; i < np; ++i) {
      if (misc.charge[i] != 0) {
	x = (PTYPE) (old.x(i) + cur.vx(i)* falpha);
      } else {
	x = (PTYPE) (old.x(i) + vx0 * falpha);
      }
      
      //     Enforce periodic boundary conditions: 
      if (x < (PTYPE) 0.) {
	x += fnx;
	misc.charge[i]=1;
      }
      
      if (x >= fnx) {
	x += -fnx;
	misc.charge[i]=1;
      }
      
      adv.x(i) = x;

#ifdef CHECK
      if (adv.x(i) >= nx || adv.x(i) < (PTYPE)0. ) raise(SIGTRAP);
#endif

    } // END for (i = 0; i < np; ++i) 
    
  } else {
    // Update with the adv velocity: 
    for (i = 0; i < np; ++i) {
      if (misc.charge[i] != 0) {
	x = (PTYPE) (old.x(i) + cur.vx(i));
      } else {
	x = (PTYPE) (old.x(i) + vx0);
      }

      //     Enforce periodic boundary conditions: 
      if (x < (PTYPE)  0.) {
	x += fnx;
	misc.charge[i]=1;
      }
      
      if (x >= fnx) {
	x += -fnx;
	misc.charge[i]=1;
      }
      
      adv.x(i) = x;
      
#ifdef CHECK
      if (adv.x(i) >= fnx || adv.x(i) < (PTYPE)0. ) raise(SIGTRAP);
#endif

    } // END for (i = 0; i < np; ++i) 
  }


  
  return;
} // xpush1D 

