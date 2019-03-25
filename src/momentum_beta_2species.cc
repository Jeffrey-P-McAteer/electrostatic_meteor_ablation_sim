/* This will either do just a momentum scatter or ionizing scatter of a particle.
   Two cross sections will need to be calculated
*/

#include <csignal> /* for notify_and_trap */
#include <math.h>
#include "eppic.h"
#include "eppic-mpi.h"

FTYPE get_beta(FTYPE v, int beta_model);
void set_momentum_crosssec(FTYPE &crosssec, FTYPE vrel, int type, FTYPE lsqsim_to_msq, FTYPE mass);
void fill_particle_holes(PTYPEAVec &recv, int &nrecv, particle_dist &adv, int vel_dim,
                         bool add_absent=false);

void momentum_beta_2species(int id, particle_dist *adv, FTYPE atmden_dt, FTYPE vbth, 
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
    FTYPE dWpa=0.;
    FTYPE dWpb=0.;
    int ncollide=0;
    
    // We need a list of particles that will be placed in the new particle distribution
    int ndim_place=NDIM+3; // Only works for 3-D velocities
    
    PTYPEAVec place_a(adv[id].np*ndim_place);
    int nplace_a=0;
    PTYPEAVec place_b(adv[id].np*ndim_place);
    int nplace_b=0;
    
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
	  // Calculate beta (chance of ionization)
	  FTYPE beta = get_beta(vrel*vsim_to_kmsec, beta_model[id]);

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
	  
	  // See if we should ionize the particle or keep it in the same distribution
	  if (ran3(&iran) < beta) {
	    // IONIZING COLLISION
	    // Put the a new particle in the list for placement in distribution dist_out_b
	    place_a(nplace_a++)=adv[id].x[ip];
	    if(ndim >= 2) place_a(nplace_a++)=adv[id].y[ip];
	    if(ndim >= 3) place_a(nplace_a++)=adv[id].z[ip];
	    
	    // this only works for 3-D velocities
	    place_a(nplace_a++)=(vxcm + vxf)*dt/dx;
	    place_a(nplace_a++)=(vycm + vyf)*dt/dy;
	    place_a(nplace_a++)=(vzcm + vzf)*dt/dz;

            // Also put a particle (electron persumably) into dist_out_b
	    place_b(nplace_b++)=adv[id].x[ip];
	    if(ndim >= 2) place_b(nplace_b++)=adv[id].y[ip];
	    if(ndim >= 3) place_b(nplace_b++)=adv[id].z[ip];

            // Get the electron postcollision velocity
            FTYPE vx_b, vy_b, vz_b;
            //if (e_collision_model[id] == 0){
            vx_b = (vxcm + vbth*gasdev(iran)) * (dt/dx);
            vy_b = (vycm + vbth*gasdev(iran)) * (dt/dy);
            vz_b = (vzcm + vbth*gasdev(iran)) * (dt/dz);
            //}
            /*
            else if (e_collision_model[id] == 1){
              // Get a gaussian velocity vector
              FTYPE vthmag = gasdev(iran)*vbth;
              // randomize the direction
              uz = 2*ran3(&iran) - 1.0;
              uperp = sqrt(1.0-uz*uz);
              uphi = 2*M_PI * ran3(&iran);
              ux = uperp*cos(uphi);
              uy = uperp*sin(uphi);
              vx_b = (vxcm + vthmag*ux)*dt/dx;
              vy_b = (vycm + vthmag*uy)*dt/dy;
              vz_b = (vzcm + vthmag*uz)*dt/dz;
            }
            else{
              printf("Unknown e_collision_model (%d) for distribution %d!!!\n",
                     e_collision_model[id], id);
              terminate(-1, "ERROR!!");
            }
            */
	    // this only works for 3-D velocities
	    place_b(nplace_b++) = vx_b;
	    place_b(nplace_b++) = vy_b;
	    place_b(nplace_b++) = vz_b;
            
            // Put these particles in the absent array
            if (adv[id].nabsent >= adv[id].absent.size()){
              terminate (-1,"Error: Array absent in xpush_domain too small: increase part_pad");
            }
            adv[id].absent(adv[id].nabsent++)=ip;
            adv[id].x(ip)=fabsent2;
            if (NDIM >= 2) adv[id].y(ip)=fabsent2;
            if (NDIM >= 3) adv[id].z(ip)=fabsent2;
            // Ensure that it does not move (also removes the energy)
            adv[id].vx(ip)=0.0;
            if (vel_dim >= 2) adv[id].vy(ip)=0.0;
            if (vel_dim == 3) adv[id].vz(ip)=0.0;
            adv[id].np_all--;
            // Keep track of the number which collide
            ncollide++;
            // Keep track of the system-wide energy change (add final energy):
            dWpa += Sqr(vxcm+vxf) + Sqr(vycm+vyf) + Sqr(vzcm+vzf);
            dWpb += Sqr(vx_b) + Sqr(vy_b) + Sqr(vz_b);
          }
	}
      }
    }
  
    // Keep track of the system-wide energy change (put in mks units):
    // The particle energy is m/2*n0*mean(v^2) so I have to normalize by
    // n0*area/np to get an energy.
    dW += (dWp*adv[id].m + dWpa*adv[coll_create_a_id[id]].m + 
           dWpb*adv[coll_create_b_id[id]].m)/2.;

    if (nplace_a > 0 ) {
      fill_particle_holes(place_a, nplace_a, adv[coll_create_a_id[id]], 3);
    }

    if (nplace_b > 0 ) {
      fill_particle_holes(place_b, nplace_b, adv[coll_create_b_id[id]], 3);
    }
  }
}
