/* This routine lays down a distribution of particles 
   given a function describing the density at all points 
   within the space.  
   It approximates the function with a series of piecewise-linear
   segments so the function must be smooth compared on the scale
   of the inter-particle spacing.
*/
#include <iostream>
#include <cmath>
#include "eppic.h"

void pwluniform(PTYPEAVec &vec, FTYPE x0, FTYPE x1, int &iran,
		     FTYPE (den_fun)(FTYPE x))
{
    FTYPE ran3(int *);
  // Let's divide this space up so there are MIN_PART_IN_SEG 
  // particles per segment.
  const int MIN_PART_IN_SEG=128;

  // We need to integrate the function between start and finish (x0, x1)
  // to determine the desired particle density.
  int nx=(vec.length()/MIN_PART_IN_SEG);
  FTYPE dx=(x1-x0)/nx;
  FTYPE area=(den_fun(x0)+den_fun(x1))/2.;
  FTYPE x=x0;
  for (int i=1; i<nx; i++) {
    x += dx;
    area += den_fun(x);
  }

  // The number of particles per area is:
  FTYPE part_per_area=vec.length()/area;

  // Lay down the particles:
  x=x0;
  int ic = 0; //Particle counter
  FTYPE part_remainder=0.; // fraction of particle left over from previous set
  FTYPE a = den_fun(x);
  FTYPE b = den_fun(x+dx);

  // Loop through the segments
  for (int ix=0; ix<nx; ix++) {
    FTYPE slope=(b-a)/dx;
    FTYPE c=(a+b-slope*(dx))/2.;
    int npart;
    // Determine the number of particles to place plus 
    // the part_remainder 
    FTYPE part_per_length=part_per_area*(a+b)/2;
    
    // If this is not the last segment, use npart times a random 
    // factor (range:0-2) to avoid periodicity.
    if (ix < nx-1) 
      npart = int(part_per_length + part_remainder*ran3(&iran)*2);
    else {
      // If this is the last segment, let's make sure we 
      // have the correct number of particles.
      npart = int(part_per_length + part_remainder);
      if (abs(npart + ic - int(vec.length()) ) <=  1)
	npart = (vec.length()) - ic;
      else
	terminate(1,"Correct number of particles +/-1 not found in pwluniform");
    }      
      part_remainder += part_per_length - npart;
    
    
    // If a==b then we have a uniform density and we should just 
    // lay them down uniformly
    if (a == b || fabs( slope*dx ) < 1e-6) {
      // Place the particles with a little randomness added
      for (int i=0; i<npart ;i++) vec[ic++]=x+dx*(i+ran3(&iran))/npart;
      
    } else { // Place particles in a linearly varying region:
      FTYPE iscale=(1./npart)*(slope*dx/2.+c)*dx;
      for (int i=0; i<npart ;i++) {
	vec[ic++]=(-c+sqrt(Sqr(c)+2*slope*(i*iscale)))/slope+x;
      }
    }
    // Update x for the next iteration
    x+=dx;
    a=b;
    b=den_fun(x+dx);
  }

  if (ic != vec.length()) 
    terminate(1,"Unable to obtain correct number of particles in pwluniform");

}
