/* Coulomb collision routine:

   This code collides charged particles in the small angle scattering
   framework of the Fokker-Planck equation. The ei routine assumes the
   scattering is for electrons off of infinitely massive, stationary
   ions. This assumption greatly simplfies the Rosenbluth potential 
   integral into a rational function of the electron velocity. This
   routine is called at the end of vadvance

   Inputs:
   Distribution of particle velocities: vp
   Random number seed: iran
   Colliding species number: sp, may not need this
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

void coulomb_scatter_ei(particle_dist &dist, int &iran, int sp, int fp)
{

  FTYPE lambda, nf;
  long botch, good, iptry, nptry;

  // Allows a different field particle density to be used to scale
  // the collision frequency
  if (cc_den[sp] < 0) { nf = dist.n0avg;}
  else { nf = cc_den[sp];}
  FTYPE nrat = pow(dist.n0avg/nf, 0.5);

  botch = 0;
  good = 0;


  /////////////////////////////////
  // For each particle, collide. //
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
    FTYPE v_abs;
    v_abs = sqrt(Sqr(dist.vx[ip]*dx/dt)+Sqr(dist.vy[ip]*dy/dt)+Sqr(dist.vz[ip]*dz/dt));
    if(v_abs == 0) v_abs = 1e-30;  // Divide by 0 if quiet start is used

    // Calculate the Coulomb logarithm
    lambda = log(4*PI*eps*dist.m*Sqr(v_abs)*l_debye*nrat/(dist.q*sign(dist.q)*qd[fp]*sign(qd[fp]))); // Correct
    lambda = max(1.0, lambda);
    //    lambda = log(4*PI*eps*dist.m*Sqr(l_debye)*1/(dist.q*sign(dist.q)*qd[fp]*sign(qd[fp]))); // Static CL

    // Calculate the diffusion coefficient, d11, and Langevin coefficients q1, q2
    FTYPE d11, q1_coeff, q2_coeff;
    d11 = (nf*Sqr(dist.q)*Sqr(qd[fp])*lambda/(4*PI*Sqr(eps)*Sqr(dist.m)*v_abs));
    FTYPE stdev = sqrt(d11*dt*coulomb_subcycle[sp]);

    // Sample Q1 and Q2 from 2D Gaussian defined by d11 and d33 
    q1_coeff = normal_dev(0, stdev, &iran);
    q2_coeff = normal_dev(0, stdev, &iran);

    // ROTATE INTO FRAME WHERE v_abs = v3
    FTYPE phi, th, sinth, costh, sinphi, cosphi;    // For rotating into frame of reference where v = v3
    phi = atan2(dist.vy[ip]*dy, dist.vx[ip]*dx);
    sinphi = sin(phi);
    cosphi = cos(phi);
    //old    th = atan2((cosphi*dist.vx[ip]*dx + sinphi*dist.vy[ip]*dy), dist.vz[ip]*dz);
    th = acos(dist.vz[ip]*dz/dt/v_abs);
    sinth = sin(th);
    costh = cos(th);

    // Apply change in velocity due to collision, in reference frame where v = v3
    FTYPE v1r, v2r, v3r;  // velocity components in rotated frame
    FTYPE arg_v3r=Sqr(v_abs) - Sqr(q1_coeff) - Sqr(q2_coeff);    
    
    if(arg_v3r >= 0){
      good += 1;
      v1r = q1_coeff;
      v2r = q2_coeff;
      // Calculate v3 (parallel component) such that energy is exactly conserved 
      v3r = sqrt(arg_v3r);
    }
    else {
      botch += 1;
      FTYPE angle;
      // This else statement catches particles trying to scatter to large angles. The theory
      // does not handle these, so we just deflect the particle 90 degrees in a random 
      // azimuthal direction while conserving momentum.
      angle = ran3_inline(&iran)*2*PI;
      v1r = v_abs*sin(angle);
      v2r = v_abs*cos(angle);
      v3r = 0;

      /*      v1r = 0;  // For passing problem particles without colliding
      v2r = 0;
      v3r = v_abs;*/
      
    }

    // ROTATE BACK INTO ORIGINAL FRAME and update velocities

    dist.vx[ip] = (v1r*costh*cosphi - v2r*sinphi + v3r*sinth*cosphi)*dt/dx;
    dist.vy[ip] = (v1r*costh*sinphi + v2r*cosphi + v3r*sinth*sinphi)*dt/dy;
    dist.vz[ip] = (-v1r*sinth + v3r*costh)*dt/dz;

  } // end particle for loop
  // Debug/diagnostic flags to show how many particles hard scattered 90 degrees
  //  if(mpi_rank == 0) printf("\n %li out of %li collisions botched \n", botch, nptry);
  //  if(mpi_rank == 0) printf("\n %i out of %i collisions good \n", good, nptry);
  //  if(mpi_rank == 0) printf("\n Bad angles for %li out of %li collisions \n", bada, nptry);

}
   
