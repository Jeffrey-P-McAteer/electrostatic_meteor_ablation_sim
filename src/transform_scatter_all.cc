/* A scattering routine that transforms a particle from one type to another by
   making it dissappear from one distribution and reappear in another one (or two):
   
   This code checks every particle and calculates the collision probability.
   If then determines whether or not to collide a particle with a neutral.
   It implicitly assumes the colliding particles are 
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
      
void transform_scatter_all(particle_dist &dist, FTYPE atmden_xsec_dt, FTYPE *v0_neutral,
			   FTYPE *vth_neutral, FTYPE m_neutral,
			   FTYPE &vrelmax, int &iran, FTYPE &dW, 
			   particle_dist &dist_out)
{

  // Make sure there are particles to collide
  if (dist.np > 0) {
    const int vel_dim = 3;  // This routine only works for 3 velocity dimensions

    FTYPE vth_n=sqrt( Sqr(vth_neutral[0])+Sqr(vth_neutral[1])+Sqr(vth_neutral[2]) );
    FTYPE v0_n=sqrt( Sqr(v0_neutral[0])+Sqr(v0_neutral[1])+Sqr(v0_neutral[2]) );
    FTYPE vxnp = v0_neutral[0];
    FTYPE vynp = v0_neutral[1];
    FTYPE vznp = v0_neutral[2];
    
    // Keep track of the system-wide energy change (initialize to zero):
    FTYPE dWp=0.;
    int ncollide=0;
    
    // We need a list of particles that will be placed in the new particle distribution
    int ndim_place=NDIM+3; // Only works for 3-D velocities
    
    PTYPEAVec place(dist.np*ndim_place);
    int nplace=0;
    
    //  For each of these particles...
    for (int ip = 0; ip < dist.np; ip++) {
      
      // Make sure that the particle is not absent
      if (dist.x(ip) > fabsent2/4) {
	
	// Find the relative velocity between the particle and neutral 
	FTYPE vxcur = dist.vx[ip]*dx/dt;
	FTYPE vycur = dist.vy[ip]*dy/dt;
	FTYPE vzcur = dist.vz[ip]*dz/dt;
	
	FTYPE vxrel = vxcur - vxnp;
	FTYPE vyrel = vycur - vynp;
	FTYPE vzrel = vzcur - vznp;
	FTYPE vrel =  sqrt(vxrel*vxrel + vyrel*vyrel + vzrel*vzrel);
	
	// Determine whether to scatter the particle:
	if (ran3(&iran) < (1 - exp(-vrel*atmden_xsec_dt))) {
	  
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
	  //	  dist.vx[ip] = (vxcm + vxf)*dt/dx;
	  //	  dist.vy[ip] = (vycm + vyf)*dt/dy;
	  //	  dist.vz[ip] = (vzcm + vzf)*dt/dz;
	  
	  // Put the current particle in the list for placement in distribution dist_out
	  place(nplace++)=dist.x[ip];
	  if(ndim >= 2) place(nplace++)=dist.y[ip];
	  if(ndim >= 3) place(nplace++)=dist.z[ip];
	  
	  // this only works for 3-D velocities
	  place(nplace++)=(vxcm + vxf)*dt/dx;
	  place(nplace++)=(vycm + vyf)*dt/dy;
	  place(nplace++)=(vzcm + vzf)*dt/dz;
	  
	  // Put these particles in the absent array
	  if (dist.nabsent >= dist.absent.size()) 
	    terminate (-1,"Error: Array absent in xpush_domain too small: increase part_pad");
	  dist.absent(dist.nabsent++)=ip;
	  dist.x(ip)=fabsent2;
	  if (NDIM >= 2) dist.y(ip)=fabsent2;
	  if (NDIM >= 3) dist.z(ip)=fabsent2;
	  // Ensure that it does not move (also removes the energy)
	  dist.vx(ip)=0.0;
	  if (vel_dim >= 2) dist.vy(ip)=0.0;
	  if (vel_dim == 3) dist.vz(ip)=0.0;
	  dist.np_all--;
      
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
    
    if (nplace > 0) {
      void fill_particle_holes(PTYPEAVec &recv, int &nrecv, particle_dist &adv, int vel_dim,
                               bool add_absent=false);
      fill_particle_holes(place, nplace, dist_out,3);
    }      
  }
}
