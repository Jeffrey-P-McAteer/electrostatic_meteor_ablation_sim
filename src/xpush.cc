// Push the particle positions to the next timestep: 

#include <math.h>
#include <csignal>
//#include <signal.h>
#include "eppic.h"

void xpush(PTYPEAVec &adv, PTYPEAVec &old, PTYPEAVec &rhs, int np,
	   PTYPE fnx, FTYPE alpha)
{
  
  PTYPE x;
  // For the sake of speed, distinguish between when alpha == 1 or not: 
  if (alpha != (FTYPE) 1.) {
    PTYPE falpha = (PTYPE) alpha;
    
    // Update with the adv velocity: 
    for (int i = 0; i < np; ++i) {
      x = old(i) + rhs(i)* falpha;
      //     Enforce periodic boundary conditions: 
      if (x < (PTYPE) 0.) x += fnx;
      if (x >= fnx) x += -fnx;
      adv(i) = x;

      if ((x < -nx || x >= 2*fnx) &&  x > 2*(-fnx) ) { //Test to see if particle is moving too fast but not maked as absent (last condition)
	terminate(-1,"Error: Particle jumps more than entire mesh\n");
      }
    } // END for (i = 0; i < np; ++i) 
    
  } else {
    // Update with the adv velocity: 
    for (int i = 0; i < np; ++i) {
      x = old(i) + rhs(i);
      //     Enforce periodic boundary conditions: 
      if (x < (PTYPE)  0.) x += fnx;
      if (x >= fnx) x += -fnx;
      adv(i) = x;
      
      if ((x < -nx || x >= 2*fnx) &&  x > 2*(-fnx) ) { //Test to see if particle is moving too fast but not maked as absent (last condition)
	terminate(-1,"Error: Particle jumps more than entire mesh\n");
      }

    } // END for (i = 0; i < np; ++i) 
  }

  return;
} // xpush 

