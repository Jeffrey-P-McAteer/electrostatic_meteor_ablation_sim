// Select the particle postion update algorithm:
#include <cmath>
#include "eppic.h"
#include "eppic-mpi.h"


extern void xpush_domain(PTYPEAVec &, PTYPEAVec &, int,
                         PTYPE, FTYPE, intAVec &, int &, int);

extern void xpush(PTYPEAVec &, PTYPEAVec &, PTYPEAVec &, int, PTYPE , FTYPE);
extern void xpush_bc1(particle_dist &adv, particle_dist &old, particle_dist &cur, 
		      particle_misc &misc, FTYPE alpha, FTYPE, FTYPE);
extern void xpush1D_bc1(particle_dist &, particle_dist &, particle_dist &, 
			particle_misc &, FTYPE, FTYPE, FTYPE);
void pass_particles(int id, particle_dist  *adv, int nabsent_start, int vel_dim);

void remove_particles(PTYPEAVec &pos,particle_dist  &adv, int vel_dim, int nabsent_start,
			      bool left_edge_out, bool right_edge_out);

void inject_particles_x(particle_dist *adv, int id, FArrayND &rho);
void inject_particles_y(particle_dist *adv, int id, FArrayND &rho);
void inject_particles_z(particle_dist *adv, int id, FArrayND &rho);

void xadvance(int id, particle_dist *adv, particle_dist *old, particle_dist *cur, 
	      particle_misc *misc, fluid *fspecie, FTYPE dtime, FArrayND &rho)
{
  
  // Let's time the amount of time in this routine:
  clock_t start_time=times(&times_buf);

  //  Keep track of the starting index in absent array
  int nabsent_start=adv[id].nabsent;

  if (method[id] >= 0) {

    if (chargeon[id] == 0) {

      // Standard periodic particle advance:
      // Advance x location of particles.  If they exit domain,
      // place their indices in send_east_indices and
      // send_west_indices arrays.
      if (nsubdomains > 1 || boundary_type[0] == INJECT || boundary_type[0] == OPEN) {
        // xpush_domain can add to the absent list but doesn't change xyz,vxvyvz
	xpush_domain(adv[id].x, cur[id].vx, adv[id].np, (PTYPE) nx, dtime,
		     adv[id].absent,adv[id].nabsent, 0);
   
	if (boundary_type[0] == OPEN || boundary_type[0] == INJECT) {
          // remove_particles will change xyz,vxvyvz of all new absent particles
          remove_particles(adv[id].x, adv[id], 
                           vel_dim[id], nabsent_start, 
                           subdomain.id_number == 0, 
                           subdomain.id_number == nsubdomains-1 );
        }
        
        // Pass particles between subdomains
	if(nsubdomains > 1) 
	  pass_particles(id, adv, nabsent_start, vel_dim[id]);

        // Inject particles in the x direction if needed
        if (boundary_type[0] == INJECT){
          inject_particles_x(adv, id, rho);
        }

	nabsent_start=adv[id].nabsent; // Indicate that all particles in nabsent list have been accounted for.
      } else {
	xpush(adv[id].x, old[id].x, cur[id].vx, adv[id].np, (PTYPE) nx, dtime);
      }

      if (species_dim[id] > 1) 
        if (boundary_type[1] == OPEN || boundary_type[1] == INJECT) {
          // added yguard_size by Glenn 5/2/2018
          xpush_domain(adv[id].y, cur[id].vy, adv[id].np, (PTYPE) ny, 
		       dtime,adv[id].absent, adv[id].nabsent, 1);
          
	  remove_particles(adv[id].y, adv[id], vel_dim[id], nabsent_start, TRUE, TRUE);
          // Inject particles if we want
          if (boundary_type[1] == INJECT){
            inject_particles_y(adv, id, rho);
          }
	  nabsent_start=adv[id].nabsent;// Indicate that all particles in nabsent list have been accounted for.
	} else {
	  xpush(adv[id].y, old[id].y, cur[id].vy, adv[id].np, (PTYPE) ny, dtime);
	}
      
      if (species_dim[id] > 2) 
        if (boundary_type[2] == OPEN || boundary_type[2] == INJECT) {
          // added zguard_size by Glenn 5/2/2018
	  xpush_domain(adv[id].z, cur[id].vz, adv[id].np, (PTYPE) nz, 
		       dtime,adv[id].absent, adv[id].nabsent, 2);

	  remove_particles(adv[id].z, adv[id], vel_dim[id], nabsent_start, TRUE, TRUE);

          if (boundary_type[2] == INJECT){
            inject_particles_z(adv, id, rho);
          }
	  nabsent_start=adv[id].nabsent;// Indicate that all particles in nabsent list have been accounted for.

	} else {
	  xpush(adv[id].z, old[id].z, cur[id].vz, adv[id].np, (PTYPE) nz, dtime);
	}

    } else {
      if (species_dim[id] >1) 
	xpush_bc1(adv[id], old[id], cur[id], misc[id], dtime, 
		  vx0d[id]*(dt/dx), vy0d[id]*(dt/dy));
      else
	xpush1D_bc1(adv[id], old[id], cur[id], misc[id], dtime, 
		    vx0d[id]*(dt/dx), vy0d[id]*(dt/dy));
      
    }

    void create_particles(int id, int id2, particle_dist *adv, FTYPE creation_rate);
    if (creation_rate[id]!=0 && n0_shell[id] != 0) create_particles(id, create_pair_dist[id], 
								    adv, creation_rate[id]);

    // To create a source of particles coming from a small ablating object:
    void create_particles_ablating(int id, particle_dist *adv);
    if (creation_rate[id]!=0 && (create_vth[id] != 0 || create_v0[id] != 0) )
      create_particles_ablating(id, adv);

  // } else if (method[id] == -2) {// fluids handled elsewhere
  } else if (method[id] < 0) {// fluids handled elsewhere
  }

  //dls end//////////////////////////////////////////

  // Let's time the amount of time in this routine:
  xadvance_time += times(&times_buf)-start_time;
  

}
