// Push the particle_dist positions to the next timestep: 

#include <math.h>
//#include <signal.h>
#include "eppic.h"
#include "eppic-mpi.h"

void xpush_domain(PTYPEAVec &adv, PTYPEAVec &rhs, int np, PTYPE fnx, FTYPE alpha,
		  intAVec &absent, int &nabsent, int dim)
{
  PTYPE x;
  //#ifdef DEBUG
  //  PTYPE xshift=subdomain.id_number*nx;
  // #endif
  int i;
  bool inject_bound = FALSE;
  //if (dim==0){
  //inject_bound = (boundary_type[0] == INJECT && subdomain.id_number == 0);
  //}
  //else{
    //inject_bound = boundary_type[dim] == INJECT;
  //}
    
  // For the sake of speed, distinguish between when alpha == 1 or not: 
  if (alpha != (FTYPE) 1.) {
    PTYPE falpha = (PTYPE) alpha;
    
    // Update with the adv velocity: 
    for (i = 0; i < np; ++i) {
      x = adv(i) + rhs(i)* falpha;
      
      // place locations of particle_dists that have crossed east
      // or west borders into "absent" arrays, and add their
      // indices to the "absent" array
      // particles where x<-nx are already absent
      adv(i) = x;
      // See if particle is outside of domain and not already absent
      if ((x < (PTYPE) 0. || x>= fnx) && x > 2*(-fnx))  { 
	// In the case of injection particles on the left most boundary are not made
	// absent until they cross x=-1
	if ( inject_bound && (x > -1 && x< (PTYPE) 0.) ) {/*do nothing*/}
	else {
	  if (nabsent >= absent.size()){
            printf("%s:%d: nabsent: %d, absent.size(): %d", __FILE__, __LINE__, nabsent, absent.size());
	    terminate (-1,"Error: Array absent in xpush_domain too small: increase part_pad");
	  }
	  absent(nabsent++) = i;
	}
	if (x < -fnx || x >= 2*fnx ) { //Test to see if particle is moving too fast
          printf("Original position: %f, Final position: %f", x, x-rhs(i)*falpha);
	  terminate(-1,"Error: Particle jumps more than entire mesh");
	}
      }
    } // END for (i = 0; i < np; ++i) 
  } else {
    // Update with the adv velocity: 
    for (int i = 0; i < np; ++i) {
    //for (int i = 0; i < adv.size(); ++i) {
      // #ifdef DEBUG
      /* this fixes a minor bug with precision. The bug is best explained by
	 an example:
	 In 1 domain code:

	 rhs(i) = 17.00000;   // 0.1700000E2
	 adv(i) =  0.000005; // 0.5000000E-6
	 x = rhs(i) + adv(i) // = 0.170000005E2, but truncted to 0.1700000E2

	 In 2 domain code: (nx=16, lets say)
	 rhs(i) =  1.000000;   // 0.1000000E1
	 adv(i) =  0.000005; // 0.5000000E-6
	 x = rhs(i) + adv(i) // = 0.10000005E1, but rounded to 0.1000001E1
	 
	 shifted, this is 0.1700001E2, not the same as 1 domain code
			     */
      //      x  = xshift + adv(i);
      //      x += rhs(i);
      //      x -= xshift;
      // #else
      x  = rhs(i)+adv(i);
      // #endif
      adv(i) = x;
      
      // place locations of particles that have crossed east
      // or west borders into "absent" array
      // particles where x<-nx are already absent
      //if ((x < (PTYPE) 0. || x>= fnx) && x > 2*(-fnx))  {
      if ((x < (PTYPE) 0. || x>= fnx) && x > fabsent2)  {
	// In the case of injection particles on the left most boundary are not made
	// absent until they cross x=-1
	if ( inject_bound && (x > -1 && x< (PTYPE) 0.) ) {/*do nothing*/}
	else {	
	  if (nabsent >= absent.size()){
            printf("%s:%d: nabsent: %d, absent.size(): %d", __FILE__, __LINE__, nabsent, absent.size());
	    terminate (-1,"Error: Array absent in xpush_domain too small: increase part_pad");
	  }
	  absent(nabsent++) = i;
	}
	if (x < -fnx || x >= 2*fnx ) { //Test to see if particle is moving too fast
          printf("Original position: %f, Final position: %f", x-rhs(i), x);
	  terminate(-1,"Error: Particle jumps more than entire mesh");
	}
      }
    } // END for (i = 0; i < np; ++i) 
  }
  
  return;
} // xpush 
