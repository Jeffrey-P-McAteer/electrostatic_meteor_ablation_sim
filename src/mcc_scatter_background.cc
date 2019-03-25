/* A scattering routine that that checks every particle and calculates the 
   collision probability.  Particles will collide with a background neutral
   that has constant density and velocity.
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
  
void set_momentum_crosssec(FTYPE &crosssec, FTYPE vrel, int type, FTYPE lsqsim_to_msq, FTYPE mass);

void mcc_scatter_background(int id, particle_dist *adv, FTYPE atmden_dt,
                           int &iran, FTYPE &dW)
{

  // Make sure there are particles to collide
  if (adv[id].np > 0) {
    const int vel_dim = 3;  // This routine only works for 3 velocity dimensions

    // Added by Glenn 5/18/2018 to deal with thermal neutrals
    FTYPE vxnp = vx0d_neutral[id]+gasdev(iran)*vxthd_neutral[id];
    FTYPE vynp = vy0d_neutral[id]+gasdev(iran)*vythd_neutral[id];
    FTYPE vznp = vz0d_neutral[id]+gasdev(iran)*vzthd_neutral[id];
    
    // Keep track of the system-wide energy change (initialize to zero):
    FTYPE dWp=0.;
    int ncollide=0;
    
    // We need a list of particles that will be placed in the new particle distribution
    int ndim_place=NDIM+3; // Only works for 3-D velocities
    
    PTYPEAVec place(adv[id].np*ndim_place);
    int nplace=0;

    FTYPE M = adv[id].m+massd_neutral[id];
    FTYPE cross_sec = coll_cross_section[id];

    //  For each of these particles...
    for (int ip = 0; ip < adv[id].np; ip++) {      
      // Make sure that the particle is not absent
      if (adv[id].x(ip) > fabsent2) {
	// Find the relative velocity between the particle and neutral 
        // Convert to length/time from cell_size/time_step
	FTYPE vxcur = adv[id].vx[ip]*dx/dt;
	FTYPE vycur = adv[id].vy[ip]*dy/dt;
	FTYPE vzcur = adv[id].vz[ip]*dz/dt;
	FTYPE vxrel = vxcur - vxnp;
	FTYPE vyrel = vycur - vynp;
	FTYPE vzrel = vzcur - vznp;
	FTYPE vrel =  sqrt(vxrel*vxrel + vyrel*vyrel + vzrel*vzrel);
	// Calculate the cross sections
        if (crosssec_m_model[id] == 1){
          // convert velocity to km/sec for Jones formula
          set_momentum_crosssec(cross_sec, vrel*vsim_to_kmsec,
                                crosssec_m_model[id], lsqsim_to_msq, adv[id].m);
        }
        else{
          // convert velocity to m/sec
          set_momentum_crosssec(cross_sec, vrel*vsim_to_kmsec/1000.0,
                                crosssec_m_model[id], lsqsim_to_msq, adv[id].m);
        }

	// Determine whether to scatter the particle:
	if (ran3(&iran) < (1 - exp(-vrel*cross_sec*atmden_dt))) {
	  // Go to the center-of-mass frame 
	  FTYPE vxcm, vycm, vzcm, vi;
	  vxcm = (vxcur*adv[id].m + vxnp*massd_neutral[id])/M;
	  vycm = (vycur*adv[id].m + vynp*massd_neutral[id])/M;
	  vzcm = (vzcur*adv[id].m + vznp*massd_neutral[id])/M;
	  
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
	  dWp -= Sqr(vxcur) + Sqr(vycur) + Sqr(vzcur);
	  
          adv[id].vx[ip] = (vxcm + vxf)*dt/dx;
          adv[id].vy[ip] = (vycm + vyf)*dt/dx;
          adv[id].vz[ip] = (vzcm + vzf)*dt/dx;
#ifdef CHECK
	  void notify_and_trap( const char *message);
	  if (adv[id].vx[ip] > (vth_neutral[0]*dt/dx)*sqrt(m_neutral/adv[id].m)*20) 
	    notify_and_trap("Error: Particle found going over 20*thermal after collision in momentum_beta_scatter.cc");
	  if (adv[id].vy[ip] > (vth_neutral[1]*dt/dy)*sqrt(m_neutral/adv[id].m)*20) 
	    notify_and_trap("Error: Particle found going over 20*thermal after collision in momentum_beta_scatter.cc");
	  if (adv[id].vz[ip] > (vth_neutral[2]*dt/dz)*sqrt(m_neutral/adv[id].m)*20) 
	    notify_and_trap("Error: Particle found going over 20*thermal after collision in momentum_beta_scatter.cc");
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
    dW += dWp*adv[id].m/2.;
  }
}


/*
void mcc_scatter_background(particle_dist &dist, FTYPE atmden_dt, FTYPE crosssec,
			    int crosssec_m_model, FTYPE vsim_to_kmsec,
			    FTYPE lsqsim_to_msq,
			    FTYPE *v0_neutral, FTYPE *vth_neutral, FTYPE m_neutral,
			    FTYPE &vrelmax, int &iran, FTYPE &dW)
{
  
  // Make sure there are particles to collide
  if (dist.np > 0) {
    void set_momentum_crosssec(FTYPE &crosssec, FTYPE vrel, int type, FTYPE lsqsim_to_msq, FTYPE mass);
    
    const int vel_dim = 3;  // This routine only works for 3 velocity dimensions
    
    FTYPE vxnp = v0_neutral[0]+gasdev(iran)*vth_neutral[0];
    FTYPE vynp = v0_neutral[1]+gasdev(iran)*vth_neutral[1];
    FTYPE vznp = v0_neutral[2]+gasdev(iran)*vth_neutral[2];
    
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
      if (dist.x(ip) > fabsent2) {
	
	// Find the relative velocity between the particle and neutral 
	FTYPE vxcur = dist.vx[ip]*dx/dt;
	FTYPE vycur = dist.vy[ip]*dy/dt;
	FTYPE vzcur = dist.vz[ip]*dz/dt;
	
	FTYPE vxrel = vxcur - vxnp;
	FTYPE vyrel = vycur - vynp;
	FTYPE vzrel = vzcur - vznp;
	FTYPE vrel =  sqrt(vxrel*vxrel + vyrel*vyrel + vzrel*vzrel);
	
	// Calculate the cross sections
        if (crosssec_m_model == 1){
          set_momentum_crosssec(crosssec, vrel*vsim_to_kmsec, crosssec_m_model, lsqsim_to_msq, dist.m);
        }
        else{
          set_momentum_crosssec(crosssec, vrel*vsim_to_kmsec/1000.0, crosssec_m_model, lsqsim_to_msq, dist.m);
        }
	
	// Determine whether to scatter the particle:
	if (ran3(&iran) < (1 - exp(-vrel*crosssec*atmden_dt))) {
	  
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
	  
	  // Only momentum transfer, don't change species!
	  dist.vx[ip] = (vxcm + vxf)*dt/dx;
	  dist.vy[ip] = (vycm + vyf)*dt/dy;
	  dist.vz[ip] = (vzcm + vzf)*dt/dz;

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
}
*/
    
/*void mcc_scatter_background(particle_dist &dist, FTYPE atmden_xsec_dt, FTYPE *v0_neutral,
			    FTYPE *vth_neutral, FTYPE m_neutral,
			    FTYPE &vrelmax, int &iran, FTYPE &dW) 
{

  // Make sure there are particles to collide
  if (dist.np > 0) {
    FTYPE vth_n=sqrt( Sqr(vth_neutral[0])+Sqr(vth_neutral[1])+Sqr(vth_neutral[2]) );
    FTYPE v0_n=sqrt( Sqr(v0_neutral[0])+Sqr(v0_neutral[1])+Sqr(v0_neutral[2]) );
    FTYPE vxnp = v0_neutral[0];
    FTYPE vynp = v0_neutral[1];
    FTYPE vznp = v0_neutral[2];
    
    // Keep track of the system-wide energy change (initialize to zero):
    FTYPE dWp=0.;
    
    // We need a list of particles that will be placed in the new particle distribution
    int ndim_place=NDIM+3; // Only works for 3-D velocities
    
    //  For each of these particles...
    for (int ip = 0; ip < dist.np; ip++) {
      
      // Make sure that the particle is not absent
      if (dist.x(ip) > fabsent2) {
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
	  dist.vx(ip)=(vxcm + vxf)*dt/dx;
	  dist.vy(ip)=(vycm + vyf)*dt/dx;
	  dist.vz(ip)=(vzcm + vzf)*dt/dx;
	  
#ifdef CHECK
	  void notify_and_trap( const char *message);
	  if (dist.vx[ip] > (vth_neutral[0]*dt/dx)*sqrt(m_neutral/dist.m)*20) 
	    notify_and_trap("Error: Particle found going over 20*thermal after collision in mcc_scatter_background.cc");
	  if (dist.vy[ip] > (vth_neutral[1]*dt/dy)*sqrt(m_neutral/dist.m)*20) 
	    notify_and_trap("Error: Particle found going over 20*thermal after collision in mcc_scatter_background.cc");
	  if (dist.vz[ip] > (vth_neutral[2]*dt/dz)*sqrt(m_neutral/dist.m)*20) 
	    notify_and_trap("Error: Particle found going over 20*thermal after collision in mcc_scatter_background.cc");
#endif

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
}
*/
