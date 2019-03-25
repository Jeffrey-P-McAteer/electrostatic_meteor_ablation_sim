/* Coulomb collision routine:

IN DEVELOPMENT, DEFINITELY NOT WORKING

   This code collides charged particles in the small angle scattering
   framework of the Fokker-Planck equation. The general routine only
   assumes the background population being collided off of is close to
   an isotropic Maxwellian. This routine works for same-species or 
   electron-ion collisions, but those specific routines should be used
   in order to speed up computation time.

   Inputs:
   Distribution of particle velocities: vp
   Random number seed: iran
   Colliding species number: sp
   Field species number: fp

   Outputs:
   Modified velocity: vp

*/

#include <csignal>
#include <cstdio>
#include <math.h>
#include "eppic.h"
#include "eppic-math.h"
#include "eppic-mpi.h"
#include "normal_dev.h"

void coulomb_general_coll(particle_dist &dist, int &iran, int sp, int fp)
{

  FTYPE lambda, nf;    // lambda is the Coulomb Logarithm.
  int iptry, nptry;
  long botch, good;

  // I AM HERE, FIX ME
  // Correct below to use cc_den instead of cc_mass
  if (cc_den[sp] < 0) { nf = dist.n0avg;}
  else { nf = cc_den[sp];}

  ////////////
  // For each particle, collide.
  botch = 0;
  good = 0;
  nptry = dist.np/coulomb_subcycle[sp];
  for(iptry = 0; iptry < nptry; ++iptry){

    // Pick a particle at random to collide
    int ip;
    if(coulomb_subcycle[sp] == 1){
      ip = iptry;
    }
    else{
      do { ip=int(dist.np * ran3_inline(&iran)); }
      while (dist.x(ip) < 0); // This excludes absent particles which have <0
    }

    // Caculate the particle's speed
    FTYPE v_abs;   // The absolute value of a particle's velocity
    v_abs = sqrt(Sqr(dist.vx[ip]*dx/dt)+Sqr(dist.vy[ip]*dy/dt)+Sqr(dist.vz[ip]*dz/dt));

    // Define the Coulomb logarithm
    lambda = log(4*PI*eps*dist.m*Sqr(v_abs)*l_debye*(dist.n0avg/nf)/(dist.q*sign(dist.q)*qd[fp]*sign(qd[fp])));

    // Calculate the diffusion coefficients
    FTYPE d11, d33, fd, q1_coeff, q2_coeff, q3_coeff;

    FTYPE stdev = sqrt(d11*dt*coulomb_subcycle[sp]);

    // Sample Q1 and Q2 from 2D Gaussian defined by d11 and d33 
    q1_coeff = normal_dev(0, stdev, &iran);
    q2_coeff = normal_dev(0, stdev, &iran);
    q3_coeff = normal_dev(0, stdev, &iran);

    // ROTATE INTO FRAME WHERE v = v3
    FTYPE phi, th, sinth, costh, sinphi, cosphi;    // For rotating into frame of reference where v = v3
    phi = atan2(dist.vy[ip],dist.vx[ip]);
    sinphi = sin(phi);
    cosphi = cos(phi);
    th = atan2((cosphi*dist.vx[ip] + sinphi*dist.vy[ip]), dist.vz[ip]);
    sinth = sin(th);
    costh = cos(th);

    // Apply change in velocity due to collision, in reference frame where v = v3
    FTYPE v1r, v2r, v3r;  // velocity components in rotated frame
    FTYPE arg_v3r=Sqr(v_abs) - Sqr(q1_coeff) - Sqr(q2_coeff);    
    
      v1r = q1_coeff;
      v2r = q2_coeff;
      v3r = q3_coeff;

    // ROTATE BACK INTO ORIGINAL FRAME and update velocities
    dist.vx[ip] = (v1r*costh*cosphi - v2r*sinphi + v3r*sinth*cosphi)*dt/dx;
    dist.vy[ip] = (v1r*costh*sinphi + v2r*cosphi + v3r*sinth*sinphi)*dt/dy;
    dist.vz[ip] = (-v1r*sinth + v3r*costh)*dt/dz;


  } // end particle for loop

} // end collision function
   
