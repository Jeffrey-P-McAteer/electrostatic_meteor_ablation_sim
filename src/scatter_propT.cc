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

#include <csignal> /* for notify_and_trap */
#include <math.h>
#include "eppic.h"
#include "eppic-mpi.h"
      
void scatter_propT(particle_dist &dist, FTYPE coll_frac, FTYPE *v0_neutral,
		   FTYPE *vth_neutral, FTYPE m_neutral,
		   FTYPE &vrelmax, int &iran, FTYPE &dW) 
{

  // Decide how many particles to consider for collisions this
  // timestep.  (Note: could improve this by using a Poisson
  //  distribution of particle numbers.) 
  // Calculate the number of particles to collide
  // The 1.0 is a (non)geometric factor
  FTYPE vth_n=Sqr(vth_neutral[0])+Sqr(vth_neutral[1])+Sqr(vth_neutral[2]);
  vth_n=sqrt(vth_n);
  FTYPE nptry = 1.0*(dist.np-dist.nabsent)*coll_frac * vrelmax / 
    (vth_n*sqrt(m_neutral/dist.m));

  if (nptry>(dist.np-dist.nabsent)) {
    char message[128];
    sprintf(message,"%s ( %g > %Ld )\n",
	    "elastic_scatter attempting to collide more particles than present",
	    nptry, dist.np-dist.nabsent);
    terminate(-1,message);
  }

  // Keep track of the system-wide energy change (initialize to zero):
  FTYPE dWp=0.;
  int ncollide=0;

  //  For each of these particles...
  for (int iptry = 0; iptry < nptry; iptry++) {
    // If there is a residule fractional probablility that a particle will collide (nptry-iptry<1), 
    // use a random number to determine if a particle will collide.
    if ( nptry-iptry > 1. || ran3(&iran) < (nptry-iptry) ) {

      //  Pick a particle at random:
      int ip;
      do { ip=int(dist.np * ran3(&iran)); }
      while (dist.x(ip) < fabsent2/4); // This excludes absent particles which have <0

      // Choose a neutral velocity at random:
      FTYPE gasdev(int &idum);
      FTYPE vxnp=(gasdev(iran)*vth_neutral[0]+v0_neutral[0]);
      FTYPE vynp=(gasdev(iran)*vth_neutral[1]+v0_neutral[1]);
      FTYPE vznp=(gasdev(iran)*vth_neutral[2]+v0_neutral[2]);

      //    Find the relative velocity between the particle and neutral 
      FTYPE vxcur = dist.vx[ip]*dx/dt;
      FTYPE vycur = dist.vy[ip]*dy/dt;
      FTYPE vzcur = dist.vz[ip]*dz/dt;

      FTYPE vxrel = vxcur - vxnp;
      FTYPE vyrel = vycur - vynp;
      FTYPE vzrel = vzcur - vznp;
      FTYPE vrelsqr =  vxrel*vxrel + vyrel*vyrel + vzrel*vzrel;
      FTYPE vrelmaxsqr = vrelmax*vrelmax;

      // Test to make sure vrel < vrelmax 
      //    if (vrelsqr > vrelmaxsqr) {
      //      static int warn_vrelmax_change=FALSE;
      //      //      vrelmax = vrel;
      //      if (!warn_vrelmax_change) {
      //	fprintf(stderr,"Warning: vrel (%g) > vrelmax (%g) too low in elastic_scatter\n%s",
      //		sqrt(vrelsqr),vrelmax,"Message will not repeat\n\n");
      //	printf("Warning: vrel (%g) > vrelmax (%g) too low in elastic_scatter\n%s",
      //		sqrt(vrelsqr),vrelmax,"Message will not repeat\n\n");
      //    //	warn_vrelmax_change=TRUE;
      //      }
      //      }

      //     Determine whether to scatter this particle:
      if (ran3(&iran) < vrelsqr/vrelmaxsqr) {

	// Go to the center-of-mass frame 
	FTYPE vxcm, vycm, vzcm, vi;
	vxcm = (vxcur*dist.m + vxnp*m_neutral)/(dist.m+m_neutral);
	vycm = (vycur*dist.m + vynp*m_neutral)/(dist.m+m_neutral);
	vzcm = (vzcur*dist.m + vznp*m_neutral)/(dist.m+m_neutral);

	// Define the speed with respect to the c.o.m.
	vi = sqrt(Sqr(vxcur-vxcm) +
		  Sqr(vycur-vycm) + Sqr(vzcur-vzcm) );

	// Find the spherical coord. angles which take the z-axis
	// to the incident direction vi:
	FTYPE costh, sinth, cosphi, sinphi;
	FTYPE vrel=sqrt(vrelsqr);

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
	dWp -= Sqr(vxcur) + Sqr(vycur) + Sqr(vzcur);
        
	// Update the particle velocities:
	dist.vx[ip] = (vxcm + vxf)*dt/dx;
	dist.vy[ip] = (vycm + vyf)*dt/dy;
	dist.vz[ip] = (vzcm + vzf)*dt/dz;

#ifdef CHECK
	void notify_and_trap( const char *message);
	if (dist.vx[ip] > (vth_neutral[0]*dt/dx)*sqrt(m_neutral/dist.m)*20) 
	  notify_and_trap("Error: Particle found going over 20*thermal after collision in elastic_scatter.cc");
	if (dist.vy[ip] > (vth_neutral[1]*dt/dy)*sqrt(m_neutral/dist.m)*20) 
	  notify_and_trap("Error: Particle found going over 20*thermal after collision in elastic_scatter.cc");
	if (dist.vz[ip] > (vth_neutral[2]*dt/dz)*sqrt(m_neutral/dist.m)*20) 
	  notify_and_trap("Error: Particle found going over 20*thermal after collision in elastic_scatter.cc");
#endif
	// Keep track of the number which collide
	ncollide++;
	// Keep track of the system-wide energy change (add final energy):
	dWp += Sqr(vxcm+vxf) + Sqr(vycm+vyf) + Sqr(vzcm+vzf);

      }
    }
  }
  // Keep track of the system-wide energy change (put in mks units):
  // The particle energy is m/2*n0*mean(v^2) so I have to normalize by
  // n0*area/np to get an energy.
  dW += dWp*dist.m/2.;
      
}

