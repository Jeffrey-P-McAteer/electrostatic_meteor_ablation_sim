/* Elastic scattering routine:
   
   This code collides a fraction of particles with a randomly generated 
   neutral particle.  It implicitly assumes the colliding particles are 
   close to equilibrium with the neutrals.


   Input: 
   a distribution of particle velocities, vp
   a mean neutral velocity and temperature, v0_neutral and vth_neutal
   a collision rate times the timestep, coll_frac
   the fastest relative velocity, vrelmax (updated by this routine)
   random number seed, iran

   Output: 
   modified vp 
   updated vrelmax (removed 8/01)
   Calculated change in kinetic energy of vp, dW

*/

#include <math.h>
#include "eppic.h"
      
void elastic_scatter_efast(particle_dist &dist, FTYPE &coll_frac, FTYPE *v0_neutral,
		     FTYPE *vth_neutral, FTYPE m_neutral,
		     FTYPE &vrelmax, int &iran, FTYPE &dW) 
{

  // Decide how many particles to collide this
  // timestep.  
  FTYPE coll_num = (dist.np-dist.nabsent) * coll_frac;

  if (coll_num>(dist.np - dist.nabsent)) {
    char message[128];
    sprintf(message,"%s ( %g > %Ld )\n",
	    "elastic_scatter needs to collide more particles than present",
	    coll_num, dist.np-dist.nabsent);
    terminate(-1,message);
  }

  // Keep track of the system-wide energy change (initialize to zero):
  FTYPE dWp=0.;

  //  For each of these particles...
  FTYPE vx0np=v0_neutral[0];
  FTYPE vy0np=v0_neutral[1];
  FTYPE vz0np=v0_neutral[2];

  int n_collide = 0;
  int n_attempt = 0;
  while (n_collide <= coll_num && n_attempt++<(dist.np-dist.nabsent)) {
    //  Pick a particle at random:
    int ip;
    do { ip=int(dist.np * ran3(&iran)); }
    while (dist.x(ip) < fabsent2/4); // This excludes absent particles which have <0
    
    FTYPE vxcur = dist.vx[ip]*dx/dt;
    FTYPE vycur = dist.vy[ip]*dy/dt;
    FTYPE vzcur = dist.vz[ip]*dz/dt;

    //    Find the relative velocity between the particle and avg. neutral 
    FTYPE vxrel = vxcur - vx0np;
    FTYPE vyrel = vycur - vy0np;
    FTYPE vzrel = vzcur - vz0np;
    FTYPE vrelsqr = vxrel*vxrel + vyrel*vyrel + vzrel*vzrel;

    // Test to make sure vrel < vrelmax 
    //    if (vrel > vrelmax*vrelmax) {
    //      static int warn_vrelmax_change=FALSE;
    //      //      vrelmax = vrel;
    //      if (!warn_vrelmax_change) {
    //	fprintf(stderr,"Warning: vrel (%g) > vrelmax (%g) too low in elastic_scatter/n%s",
    //		sqrt(vrel),vrelmax,"Message will not repeat\n\n");
    //	printf("Warning: vrel (%g) > vrelmax (%g) too low in elastic_scatter\n%s",
    //		sqrt(vrel),vrelmax,"Message will not repeat\n\n");
    //    //	warn_vrelmax_change=TRUE;
    //      }
    //    }

    //     Determine whether to scatter this particle:
    if (Sqr(ran3(&iran)) < (vrelsqr)/(vrelmax*vrelmax)) {
      // Particle counter
      n_collide++;

      // Choose a neutral velocity at random:
      FTYPE gasdev(int &idum);
      FTYPE vxnp=(gasdev(iran)*vth_neutral[0]+v0_neutral[0]);
      FTYPE vynp=(gasdev(iran)*vth_neutral[1]+v0_neutral[1]);
      FTYPE vznp=(gasdev(iran)*vth_neutral[2]+v0_neutral[2]);

      // Recalculate the relative velocity between the particle and neutral 
      vxrel = vxcur - vxnp;
      vyrel = vycur - vynp;
      vzrel = vzcur - vznp;
      FTYPE vrel=sqrt(vxrel*vxrel + vyrel*vyrel + vzrel*vzrel);

      // Go to the center-of-mass frame 
      FTYPE vxcm, vycm, vzcm, vi;
      vxcm = (vxcur*dist.m + vxnp*m_neutral)/(dist.m+m_neutral);
      vycm = (vycur*dist.m + vynp*m_neutral)/(dist.m+m_neutral);
      vzcm = (vzcur*dist.m + vznp*m_neutral)/(dist.m+m_neutral);

      // Define the speed with respect to the c.o.m.
      vi = sqrt(Sqr(vxcur-vxcm) + Sqr(vycur-vycm) + Sqr(vzcur-vzcm) );

      // Find the spherical coord. angles which take the z-axis
      // to the incident direction vi:
      FTYPE costh, sinth, cosphi, sinphi;
      costh = vzrel / vrel;
      sinth = sqrt(1.0-costh*costh);
      cosphi = vxrel / (vrel*sinth);
      sinphi = vyrel / (vrel*sinth);

      // Choose the unit scatter vector relative to the z-axis:
      FTYPE uz, uperp, uphi, ux, uy;
      uz = 2*ran3(&iran) - 1.0;
      uperp = sqrt(1.0-uz*uz);
      uphi = 2*M_PI * ran3(&iran);
      ux = uperp*cos(uphi);
      uy = uperp*sin(uphi);

      // Rotate the unit scatter vector to the incident coord. sys.:
      FTYPE uz1, ux1, uy1, ux2, uy2, uz2, vxf, vyf, vzf;
      uz1 = uz*costh-ux*sinth;
      ux1 = uz*sinth+ux*costh;
      uy1 = uy;
      ux2 = ux1*cosphi-uy1*sinphi;
      uy2 = ux1*sinphi+uy1*cosphi;
      uz2 = uz1;
      vxf = vi * ux2;
      vyf = vi * uy2;
      vzf = vi * uz2;

      // Keep track of the system-wide energy change (subtract initial energy):
      dWp -= Sqr(dist.vx[ip]) + Sqr(dist.vy[ip]) + Sqr(dist.vz[ip]);
        
      // Update the particle velocities:
      dist.vx[ip] = vxcm + vxf;
      dist.vy[ip] = vycm + vyf;
      dist.vz[ip] = vzcm + vzf;

      // Keep track of the system-wide energy change (add final energy):
      dWp += Sqr(dist.vx[ip]) + Sqr(dist.vy[ip]) + Sqr(dist.vz[ip]);

      // Renormalize the particle velocities:
      dist.vx[ip] *= dt/dx;
      dist.vy[ip] *= dt/dy;
      dist.vz[ip] *= dt/dz;

    }
  }

  if (n_attempt>(dist.np-dist.nabsent)) {
    char message[128];
    sprintf(message,"%s ( %d > %Ld )\n",
	    "elastic_scatter tested more particles than present",
	    n_attempt, dist.np-dist.nabsent);
    terminate(-1,message);
  }

  // Pass back the collision fraction actually done:
  coll_frac=(1.0*n_collide)/(dist.np-dist.nabsent);

  // Keep track of the system-wide energy change (put in mks units):
  // The particle energy is m/2*n0*mean(v^2) so I have to normalize by
  // n0*area/np to get an energy.
  dW += dWp*dist.m/2.;

}
