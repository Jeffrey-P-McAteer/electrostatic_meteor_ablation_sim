// Push the particle positions to the next timestep with periodic B.C.s
// where the particle charges only turn on only after passing through a boundary
#include <csignal>
#include <math.h>
//#include <signal.h>
#include "eppic.h"

void xpush_bc1(particle_dist &adv, particle_dist &old, particle_dist &cur, 
	       particle_misc &misc, FTYPE alpha, FTYPE vx0, FTYPE vy0)
{
  
  int i;
  int np=old.x.length();
  PTYPE x, y;

#ifndef __GNUC__
  const PTYPE min_resolvable_value = 0.0;
#else
  // Rounding error problems in some compilers require using this value 
  const int other_bits = 12;
  const PTYPE min_resolvable_value = 1.0/pow(2.,8.*sizeof(PTYPE)-other_bits);
#endif

  PTYPE fnx = (PTYPE) nx*(1-min_resolvable_value);
  PTYPE fny = (PTYPE) ny*(1-min_resolvable_value);
    
  // For the sake of speed, distinguish between when alpha == 1 or not: 
  if (alpha != 1.) {
    PTYPE falpha = (PTYPE) alpha;
    
    // Update with the adv velocity: 
    for (i = 0; i < np; ++i) {
      if (misc.charge[i] != 0) {
	x = (PTYPE) (old.x(i) + cur.vx(i)* falpha);
	y = (PTYPE) (old.y(i) + cur.vy(i)* falpha);
      } else {
	x = (PTYPE) (old.x(i) + vx0 * falpha);
	y = (PTYPE) (old.y(i) + vy0 * falpha);
      }
      
      //     Enforce periodic boundary conditions: 
      if (x < (PTYPE) 0.) {
	x += fnx;
	misc.charge[i]=1;
      }
      if (x >= fnx) {
	x -= fnx;
	misc.charge[i]=1;
      }
      
      if (y < (PTYPE) 0.) y += fny;
      if (y >= fny) y -= fny;

      adv.x(i) = x;
      adv.y(i) = y;

#ifdef CHECK
      if (adv.x(i) >= nx || adv.x(i) < (PTYPE)0. ) raise(SIGTRAP);
      if (adv.y(i) >= ny || adv.y(i) < (PTYPE)0. ) raise(SIGTRAP);
#endif

    } // END for (i = 0; i < np; ++i) 
    
  } else {
      // Update with the adv velocity: 
    for (i = 0; i < np; ++i) {
      if (misc.charge[i] != 0) {
	x = (PTYPE) (old.x(i) + cur.vx(i));
	y = (PTYPE) (old.y(i) + cur.vy(i));
      } else {
	x = (PTYPE) (old.x(i) + vx0);
	y = (PTYPE) (old.y(i) + vy0);
      }
      
      //     Enforce periodic boundary conditions: 
      if (x < (PTYPE) 0.) {
	x += fnx;
	misc.charge[i]=1;
      }
      if (x >= fnx) {
	x -= fnx;
	misc.charge[i]=1;
      }

      if (y < (PTYPE) 0.) y += fny;
      if (y >= fny) y -= fny;

      adv.x(i) = x;
      adv.y(i) = y;

      /*
      if (adv.x(i) >= fnx || adv.x(i) < (PTYPE)0. ) raise(SIGTRAP);
      if (adv.y(i) >= fny || adv.y(i) < (PTYPE)0. ) raise(SIGTRAP);
      */
    } // END for (i = 0; i < np; ++i) 
  }


  
  return;
} // xpush 

