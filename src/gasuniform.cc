/* This subroutine returns an evenly spread Gaussian distribution in 2D.
   There will be a strong corelation between vx and vy, so shuffling the 
   arrays may be necessary to decorrelate. 

   istart is the starting point for the bitreversed distribution.

   prime1 and prime2 are two unique primes needed to insure a unique 
   distribution of uniformly distributed particles.
*/

#include <iostream>
#include <cmath>
#include "eppic.h"

void gasuniform(PTYPEAVec &vx, PTYPEAVec &vy, int istart, 
		int prime1, int prime2)
{

  //     Generate a uniform gaussian  
  int N=vx.size();
  double fac, theta, vmag;
    
  // We need a uniformly distributed circle of points with exactly N points 
  
  for (int i=0; i<N; i++) {
    vmag=reverse(i*nsubdomains+1+istart,prime1);
    if (vmag==0) terminate(1," Error in gasuniform: infinite velocity");
    theta = (2.*M_PI) * reverse(i*nsubdomains+1+istart,prime2);
    fac = sqrt ((-2.)*log(vmag));
    vx[i] = fac * cos(theta);
    vy[i] = fac * sin(theta);
  }

}


