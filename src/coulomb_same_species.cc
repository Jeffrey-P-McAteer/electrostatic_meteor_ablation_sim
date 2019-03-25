/* Coulomb collision routine:

   This code collides charged particles in the small angle scattering
   framework of the Fokker-Planck equation. This routine is designed to
   handle same species collisions. It assumes the species in question is
   close to an isotropic Maxwellian. The main difference between this
   routine and the general one is how momentum conservation is enforced,
   for same species collisions the species cannot heat up or loss energy
   when colliding with itself.
   The routine will need to be slightly modified if a beam or bump on tail
   distribution is trying to collide off the thermal population. This
   modification depends on how the beam population is entered into EPPIC.

   Inputs:
   Distribution of particle velocities: vp
   Random number seed: iran
   Colliding species number: sp

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

void coulomb_same_species(particle_dist &dist, int &iran, int sp)
{

  FTYPE lambda, nf;    // lambda is the Coulomb Logarithm.
  int iptry, nptry;

  // Allows a different field particle density to be used to obtain the
  // correct collision frequency
  if (cc_den[sp] < 0) { nf = dist.n0avg;}
  else { nf = cc_den[sp];}
  FTYPE nrat = pow(dist.n0avg/nf, 0.5);
  FTYPE cth = pow(2*dist.n0avg*dist.q*dist.q/dist.m/eps, 0.5)*l_debye;  // Calculates Sqrt(2KT/m) from l_debye


  ////////////
  // For each particle, collide.

  nptry = dist.np/coulomb_subcycle[sp];
  FTYPEAVec ipind, pc_vx, pc_vy, pc_vz;
  ipind = FTYPEAVec(nptry) = 0;
  pc_vx = FTYPEAVec(nptry) = 0;
  pc_vy = FTYPEAVec(nptry) = 0;
  pc_vz = FTYPEAVec(nptry) = 0;

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
    ipind[iptry] = ip;

    // Caculate the particle's speed
    FTYPE v_abs, vc;   // The absolute value of a particle's velocity
    v_abs = sqrt(Sqr(dist.vx[ip]*dx/dt)+Sqr(dist.vy[ip]*dy/dt)+Sqr(dist.vz[ip]*dz/dt));
    if(v_abs == 0) v_abs = 1e-30; // Divide by 0 if quiet start is used
    vc = v_abs/cth; // Ratio v/sqrt(2KT/m)

    // Define the Coulomb logarithm
    lambda = log(4*PI*eps*dist.m*Sqr(v_abs)*l_debye*nrat/Sqr(dist.q));
    lambda = max(1.0, lambda);
    lambda = nf*pow(dist.q,4)*lambda/(4*PI*v_abs*pow(eps*dist.m,2));  // Pulls ne^4/(m*eps)^2/(4pi*v) factor into lambda

    FTYPE erfv, expv; // Trying to make as few calculations as possible to get d and F
    erfv = erf(vc);
    expv = exp(-vc*vc);

    // Calculate the diffusion coefficients d11 and d33, Drag coefficient Fd, and Langevin q's 
    // Coordinate system was chosen such that d11 = d22 and there are no off-diagonal terms
    FTYPE Fd, d11, d33, q1_coeff, q2_coeff, q3_coeff;
    Fd = -2*lambda/v_abs*(erfv - 1.1283792*vc*expv);
    d11 = lambda*pow(vc,-2)*(erfv - 1.1283792*vc*expv);
    d33 = lambda*erfv -0.5*d11;

    FTYPE stdev11 = sqrt(d11*dt*coulomb_subcycle[sp]);
    FTYPE stdev33 = sqrt(d33*dt*coulomb_subcycle[sp]);

    // Sample Q1 and Q2 from 2D Gaussian defined by d11 and d33 
    q1_coeff = normal_dev(0, stdev11, &iran);
    q2_coeff = normal_dev(0, stdev11, &iran);
    q3_coeff = normal_dev(0, stdev33, &iran);

    // ROTATE INTO FRAME WHERE v = v3
    FTYPE phi, th, sinth, costh, sinphi, cosphi;    // For rotating into frame of reference where v = v3
    phi = atan2(dist.vy[ip]*dy, dist.vx[ip]*dx);
    sinphi = sin(phi);
    cosphi = cos(phi);
    th = acos(dist.vz[ip]*dz/dt/v_abs);
    sinth = sin(th);
    costh = cos(th);

    // Apply change in velocity due to collision, in reference frame where v = v3
    FTYPE v1r, v2r, v3r;  // velocity components in rotated frame    
    v1r = q1_coeff;
    v2r = q2_coeff;
    v3r = v_abs + Fd*dt*coulomb_subcycle[sp] + q3_coeff;

    // ROTATE BACK INTO ORIGINAL FRAME and update velocities
    pc_vx[iptry] = (v1r*costh*cosphi - v2r*sinphi + v3r*sinth*cosphi)*dt/dx;
    pc_vy[iptry] = (v1r*costh*sinphi + v2r*cosphi + v3r*sinth*sinphi)*dt/dy;
    pc_vz[iptry] = (-v1r*sinth + v3r*costh)*dt/dz;

  } // end particle for loop
  //  if(mpi_rank == 0) printf("ee %li collisions tried \n", nptry);


  // Calculate beta vector for momentum conservation
  FTYPE betax, betay, betaz;
  int d;
  betax = 0; betay = 0; betaz = 0;
  for(iptry = 0; iptry < nptry; iptry++){
    d = ipind[iptry];

    pc_vx[iptry] = dist.vx[d] - pc_vx[iptry];  // delta vx
    pc_vy[iptry] = dist.vy[d] - pc_vy[iptry];
    pc_vz[iptry] = dist.vz[d] - pc_vz[iptry];

    betax += pc_vx[iptry];
    betay += pc_vy[iptry];
    betaz += pc_vz[iptry];
  }
  if(mpi_rank == 5) printf("betax = %g \n", betax);
  betax = betax/nptry;
  betay = betay/nptry;
  betaz = betaz/nptry;


  // Calculate alpha constant for momentum conservation
  FTYPE anum, aden;
  anum = 0; aden = 0;
  for (iptry = 0; iptry < nptry; iptry++){
    d = ipind[iptry];
    
    pc_vx[iptry] = pc_vx[iptry] - betax;
    pc_vy[iptry] = pc_vy[iptry] - betay;
    pc_vz[iptry] = pc_vz[iptry] - betaz;

    anum += dist.vx[d]*pc_vx[iptry] + dist.vy[d]*pc_vy[iptry] + dist.vz[d]*pc_vz[iptry];
    aden += pc_vx[iptry]*pc_vx[iptry] + pc_vy[iptry]*pc_vy[iptry] + pc_vz[iptry]*pc_vz[iptry];
  }
  anum = -2*anum/aden;

  if(mpi_rank == 5) printf("Alpha = %g \n", anum);

  // Update post-collision velocities with momentum conservation
  for (iptry = 0; iptry < nptry; iptry++){
    d = ipind[iptry];
    //    if(isnan(pc_vx[iptry])) printf("pc_vx is NaN for id = %d", d);
    // if(isnan(pc_vy[iptry])) printf("pc_vy is NaN for id = %d", d);
    //if(isnan(pc_vz[iptry])) printf("pc_vz is NaN for id = %d", d);
    dist.vx[d] += anum*pc_vx[iptry];
    dist.vy[d] += anum*pc_vy[iptry];
    dist.vz[d] += anum*pc_vz[iptry];
  }

  aden = aden + 0; // FIX ME breakpoint to remove

} // end collision function
   
