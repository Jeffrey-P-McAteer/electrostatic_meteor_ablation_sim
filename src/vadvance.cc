// Advance the particle velocities:

#include <sys/times.h> //Clock related functions
#include <math.h>
#include "eppic.h"
#include "eppic-mpi.h"

/// vpush1D is not domian decomposed...
extern void vpush1D(particle_dist &, particle_dist &, particle_dist &, field &,
		    FTYPE qmrat, FTYPE alpha);
extern void vpush1D_damp(particle_dist &, particle_dist &, particle_dist &, field &, 
			 particle_misc &, FTYPE, FTYPE, FTYPE);
extern void vpush1D_bc1(particle_dist &, particle_dist &, particle_dist &, field &, 
			particle_misc &, FTYPE, FTYPE, FTYPE);

void vadvance(int id, particle_dist *adv, particle_dist *old, particle_dist *cur, 
	      particle_misc *misc, fluid *fspecie, field &Efield, FTYPE dtime)
{

  static intAVec init_collide(ndist);
  init_collide = TRUE;

  // Let's time the amount of time in this routine:
  clock_t start_time=times(&times_buf);

  if (method[id] >= 0) { 
    //PIC distribution:

    extern void vpush(particle_dist &, particle_dist &, particle_dist &, 
		      field &, FTYPE qmrat, FTYPE alpha);
    extern void vpush_domain(particle_dist &, particle_dist &, 
			     particle_dist &, field &, FTYPE qmrat, 
			     FTYPE alpha);
    extern void vpush_damp(particle_dist &, particle_dist &, particle_dist &, field &, 
			   FTYPE, FTYPE, particle_misc &, FTYPE );
    extern void vpush_bc1(particle_dist &, particle_dist &, particle_dist &, field &, 
			  FTYPE, FTYPE, particle_misc &, FTYPE);
    extern void vpushBx(particle_dist &, particle_dist &, particle_dist &, field &, 
			  FTYPE qmrat, FTYPE alpha, FTYPE);
    extern void vpushBx_domain(particle_dist &, particle_dist &, particle_dist &, field &,  
			       FTYPE qmrat, FTYPE alpha, FTYPE);
    extern void vpushBB_domain(particle_dist &, particle_dist &, particle_dist &, field &,  
			       FTYPE qmrat, FTYPE alpha, FTYPE, FTYPE, FTYPE);
    extern void vpushBz(particle_dist &, particle_dist &, particle_dist &, field &, 
			  FTYPE qmrat, FTYPE alpha, FTYPE);
    extern void vpushBz_domain(particle_dist &, particle_dist &, particle_dist &, field &, 
			  FTYPE qmrat, FTYPE alpha, FTYPE);
    
    if (pdamp_nu[id] <= 0) {
      // Initial condition problems:
      
      if (species_Bx[id] == 0 && species_Bz[id] == 0) {
 	// Unmagnetized:
	if (species_dim[id] > 1)
#ifndef USE_DOMAINS
	  vpush(adv[id], old[id], cur[id], Efield,old[id].q/old[id].m,dtime);
#else
	  vpush_domain(adv[id], old[id], cur[id], Efield,old[id].q/old[id].m,dtime);
#endif
	else // 1-D motion only (a specialized routine is necessary):
	  vpush1D(adv[id], old[id], cur[id], Efield, old[id].q/old[id].m,dtime);
      } else if (species_Bx[id] != 0) {
	// Magnetized along x:
	if (species_dim[id] > 1) {
#ifndef USE_DOMAINS
	    vpushBx(adv[id], old[id], cur[id], Efield, 
		  old[id].q/old[id].m, dtime, species_Bx[id]);
#else
	    if (species_Bz[id] == 0 && species_By[id] == 0) {  // if statement catches 3 of 4 cases when vpushBB is needed
              // Only magnetized along x
	      vpushBx_domain(adv[id], old[id], cur[id], Efield, 
                             old[id].q/old[id].m, dtime, species_Bx[id]);
	    } else {
              // Magnetized along x and (y or z)
	      vpushBB_domain(adv[id], old[id], cur[id], Efield,
                             old[id].q/old[id].m, dtime, species_Bx[id], species_By[id], species_Bz[id]);
	    }
#endif
	}
	else vpush1D(adv[id], old[id], cur[id], Efield, 
		     old[id].q/old[id].m,dtime);
      
      } else if (species_Bz[id] != 0) {
	// Magnetized along z (not x):
	if (species_dim[id] > 1) 
#ifndef USE_DOMAINS
	  vpushBz(adv[id], old[id], cur[id], Efield, 
		  old[id].q/old[id].m, dtime, species_Bz[id]);
#else
        if (species_By[id] == 0){
          // Only magnetized along x)
          vpushBz_domain(adv[id], old[id], cur[id], Efield, 
                         old[id].q/old[id].m, dtime, species_Bz[id]);
        } else {
          // Magnetized along z and y
          vpushBB_domain(adv[id], old[id], cur[id], Efield,
                         old[id].q/old[id].m, dtime, species_Bx[id], species_By[id], species_Bz[id]);
        }        
#endif
	else 
	  vpush1D(adv[id], old[id], cur[id], Efield, old[id].q/old[id].m,dtime);

      } else if (species_By[id] != 0) {
        // Magnetized along y only
        vpushBB_domain(adv[id], old[id], cur[id], Efield,
                       old[id].q/old[id].m, dtime, species_Bx[id], species_By[id], species_Bz[id]);
      }
    } else {  	

      //      terminate(1,"Injection methods not yet implemented in 3-D");

      // Non- Initial condition problems:      

      if (species_dim[id] >1) {
	if (species_Bz[id] != 0) 
	  terminate(1,"Injection not implemented for Bz != 0");

	if (species_Bx[id] == 0) { 	      
	  // UnMagnatized:
	  
	  if (chargeon[id] == 0)  
	    // Particles charged at t=0:
	    vpush_damp(adv[id], old[id], cur[id], Efield, old[id].q/old[id].m,
		      dtime, misc[id], pdamp_nu[id]/vxthd[id]*dx);
	  else            
	    // Some particles not charged at t=0:
	    vpush_bc1(adv[id], old[id], cur[id], Efield, old[id].q/old[id].m, 
		      dtime, misc[id], pdamp_nu[id]/vxthd[id]*dx);
	} else {
	  terminate(1,"This injection method not yet implemented in 3-D");
	  // Magnetized:
	  /*
	  if (chargeon[id] == 0) 
	    vpush2dBx_damp(adv[id], old[id], cur[id], Efield, misc[id], 
			 pdamp_nu[id]/vxthd[id]*dx,
			 old[id].q/old[id].m, dtime, species_Bx[id]);
	  else
	    vpush2dBx_bc1(adv[id], old[id], cur[id], Efield, misc[id], 
			pdamp_nu[id]/vxthd[id]*dx, 
			old[id].q/old[id].m, dtime, species_Bx[id]);
	  */
	}
	

      } else { //if (species_dim[id] !> 1)
	if (chargeon[id] == 0)  
	  // Particles charged at t=0:
	  vpush_damp(adv[id], old[id], cur[id], Efield, old[id].q/old[id].m,
		     dtime, misc[id], pdamp_nu[id]/vxthd[id]*dx);
	else            
	  // Some particles not charged at t=0:
	  vpush_bc1(adv[id], old[id], cur[id], Efield, old[id].q/old[id].m, 
		    dtime, misc[id], pdamp_nu[id]/vxthd[id]*dx);

      } //END if (species_dim[id] !> 1)
      
    }

    
  } else if (method[id] <0 ) {
    // fluid method - handled in fluid_advance}
    
  } else {
    terminate(-1,"Method not know in leapadv_subcycle.cc");
  }

  // Collision algorithms which only modify the velocity go here:
  if (coll_rate[id] > 0. &&  it*dt > coll_start_time[id] ) {
    if (init_collide[id]) {
      init_collide[id] = FALSE;
      FTYPE v0=sqrt( Sqr((vx0d[id]-vx0d_neutral[id])) + 
		     Sqr((vy0d[id]-vy0d_neutral[id])) + 
		     Sqr((vz0d[id]-vz0d_neutral[id])) );
      FTYPE vthd=sqrt( Sqr(vxthd[id]) + Sqr(vythd[id]) + Sqr(vzthd[id]) );
      

      FTYPE vr=v0+4.0*vthd+v0_shell[id];
      if (vrelmax(id) == 0) { vrelmax(id)=vr; }
      else if (vrelmax(id) < vr && mpi_rank==0) 
	cout << "WARNING: vrelmax(" << id << ") = " << vrelmax(id) << " less than "
	     << vr <<", the smallest recommended quantity.\n" ;
			      
    }
    static int iran=-3248-mpi_rank;
    void elastic_scatter(particle_dist &dist, FTYPE coll_frac, FTYPE *v0_neutral,
			 FTYPE *vth_neutral, FTYPE m_neutral,
			 FTYPE &vrelmax, int &iran, FTYPE &dW);

    void elastic_scatter_efast(particle_dist &dist, FTYPE &coll_frac, FTYPE *v0_neutral,
			       FTYPE *vth_neutral, FTYPE m_neutral,
			       FTYPE &vrelmax, int &iran, FTYPE &dW) ;

    void scatter_maxwell(particle_dist &dist, FTYPE coll_frac, FTYPE *v0_neutral,
			 FTYPE *vth_neutral, FTYPE m_neutral, int &iran, FTYPE &dW);

    void scatter_propT(particle_dist &dist, FTYPE coll_frac, FTYPE *v0_neutral,
			 FTYPE *vth_neutral, FTYPE m_neutral,
			 FTYPE &vrelmax, int &iran, FTYPE &dW);

    FTYPE energy_change=0;
    FTYPE coll_frac=coll_rate[id]*dt*dtime;
    FTYPE v0_n2[3]={vx0d_neutral[id],vy0d_neutral[id],vz0d_neutral[id]};
    FTYPE vth_n2[3]={vxthd_neutral[id],vythd_neutral[id],vzthd_neutral[id]};
    
    // We have four elastic scatter routines:
    // 0) For similar mass collisions & coll. rate prop. to vel. (constant cross-section)
    // 1) For similar mass collisions & coll. rate constant for all vel. 
    //                                  (cross-section decreases iwth vel).
    // 2) For collisions between particle and neutral of very different masses & cont. cross-section
    // 3) For collisions where the coll_rate should be proportional to T.
    // 4) Ionizing collisions where a subset of particles are taken and ionized
    // 5) Ionizing collisions where collision probability is calculated for individual particles
    // 6) Ionizing collisions where an ionization event creates two particles
    // 7) Ionizing collisions using Yakov's approximation
    // 8) MCC scatter against a background of neutrals (no ionization)
    // 9) Ionizing collisions using a model for beta (ionization probability), requires beta_model
    // 10) Ionizing collisions using a model for beta (ionization probability), requires beta_model
    //     producing two particles
    // Choose between them based on the mass ratio:
    enum {ELASTIC, MAXWELL, FAST_ELASTIC, NU_PROP_T, IONIZING, IONIZING_ALL, 
	  IONIZING_ALL_2SPECIES, IONIZING_ALL_YAKOVAPPROX, MCC_SCATTER, 
          MOMENTUM_BETA, MOMENTUM_BETA_2SPECIES};

    if (coll_type[id] == ELASTIC && massd_neutral[id]/adv[id].m > 1000) coll_type[id]=FAST_ELASTIC;
    if (coll_type[id] ==  ELASTIC) {
      elastic_scatter(adv[id], coll_frac, v0_n2, vth_n2, 
		      massd_neutral[id], vrelmax[id], iran, energy_change);
    } else if (coll_type[id] ==  FAST_ELASTIC) {
      elastic_scatter_efast(adv[id], coll_frac, v0_n2, vth_n2, 
			    massd_neutral[id], vrelmax[id], iran, energy_change);
    } else if (coll_type[id] ==  MAXWELL) {
      scatter_maxwell(adv[id], coll_frac, v0_n2, vth_n2, 
		      massd_neutral[id], iran, energy_change);
    } else if (coll_type[id] ==  NU_PROP_T) {
      scatter_propT(adv[id], coll_frac, v0_n2, vth_n2, 
		    massd_neutral[id], vrelmax[id], iran, energy_change);
    } else if (coll_type[id] == IONIZING) {
      void transform_scatter(particle_dist &dist, FTYPE coll_frac, FTYPE *v0_neutral,
			     FTYPE *vth_neutral, FTYPE m_neutral,
			     FTYPE &vrelmax, int &iran, FTYPE &dW, particle_dist &dist_out);

      transform_scatter(adv[id], coll_frac, v0_n2, vth_n2, massd_neutral[id], vrelmax[id],
			iran, energy_change, adv[coll_create_id[id]]);
    } else if (coll_type[id] == IONIZING_ALL) {
      void transform_scatter_all(particle_dist &dist, FTYPE atmden_xsec_dt, FTYPE *v0_neutral,
				 FTYPE *vth_neutral, FTYPE m_neutral,
				 FTYPE &vrelmax, int &iran, FTYPE &dW, 
				 particle_dist &dist_out);
      coll_frac=background_neutral_dens[id]*coll_cross_section[id]*dt*dtime;
      transform_scatter_all(adv[id], coll_frac, v0_n2, vth_n2, massd_neutral[id], vrelmax[id],
			    iran, energy_change, adv[coll_create_id[id]]);
    } else if (coll_type[id] == IONIZING_ALL_2SPECIES) {
      void transform_scatter_all_2species(particle_dist &dist, FTYPE atmden_xsec_dt, 
					  FTYPE *v0_neutral, FTYPE *vth_neutral,
					  FTYPE m_neutral, FTYPE &vrelmax, int &iran, 
					  FTYPE &dW, particle_dist &dist_a_out,
					  particle_dist &dist_b_out);
      coll_frac=background_neutral_dens[id]*coll_cross_section[id]*dt*dtime;
      transform_scatter_all_2species(adv[id], coll_frac, v0_n2, vth_n2,
				     massd_neutral[id], vrelmax[id],
				     iran, energy_change, adv[coll_create_a_id[id]],
				     adv[coll_create_b_id[id]]);
    } else if (coll_type[id] == IONIZING_ALL_YAKOVAPPROX) {
      void transform_scatter_all_yakovapprox(particle_dist &dist, FTYPE atmden_xsec_dt, 
					     FTYPE *v0_neutral,
					     FTYPE *vth_neutral, FTYPE m_neutral,
					     FTYPE &vrelmax, int &iran, FTYPE &dW, 
					     particle_dist &dist_out);
      coll_frac=background_neutral_dens[id]*coll_cross_section[id]*dt*dtime;
      transform_scatter_all_yakovapprox(adv[id], coll_frac, v0_n2, vth_n2, 
					massd_neutral[id], vrelmax[id],
					iran, energy_change, adv[coll_create_id[id]]);
    } else if (coll_type[id] == MCC_SCATTER) {
      void mcc_scatter_background(int id, particle_dist *adv, FTYPE atmden_dt,
                                  int &iran, FTYPE &dW);      
      coll_frac=background_neutral_dens[id]*dt*dtime;
      mcc_scatter_background(id, adv, coll_frac, iran, energy_change);
      /*
      void mcc_scatter_background(particle_dist &dist, FTYPE atmden_dt, FTYPE crosssec,
				  int crosssec_m_model, FTYPE vsim_to_kmsec,
				  FTYPE lsqsim_to_msq,
				  FTYPE *v0_neutral, FTYPE *vth_neutral, FTYPE m_neutral,
				  FTYPE &vrelmax, int &iran, FTYPE &dW);
      coll_frac=background_neutral_dens[id]*dt*dtime;
      mcc_scatter_background(adv[id], coll_frac, coll_cross_section[id], crosssec_m_model[id], 
			     vsim_to_kmsec, lsqsim_to_msq, v0_n2, vth_n2, 
			     massd_neutral[id], vrelmax[id],
			     iran, energy_change);
      */
    } else if (coll_type[id] == MOMENTUM_BETA) {
      void momentum_beta_scatter(int id, particle_dist *adv, FTYPE atmden_dt,
                                 int &iran, FTYPE &dW);      
      coll_frac=background_neutral_dens[id]*dt*dtime;
      momentum_beta_scatter(id, adv, coll_frac, iran, energy_change);
    } else if (coll_type[id] == MOMENTUM_BETA_2SPECIES) {
      void momentum_beta_2species(int id, particle_dist *adv, FTYPE atmden_dt, FTYPE vbth,
                                  int &iran, FTYPE &dW);
      coll_frac=background_neutral_dens[id]*dt*dtime;
      momentum_beta_2species(id, adv, coll_frac, coll_create_vthb[id], iran, energy_change);
    }
  }


  // Coulomb collision routine goes here
  static int iran=-3248-cc_rand-mpi_rank;
  void coulomb_scatter_ei(particle_dist &dist, int &iran, int sp, int fp);
  void coulomb_same_species(particle_dist &dist, int &iran, int sp);
  void coulomb_general_coll(particle_dist &dist, int &iran, int sp, int fp);

  int icoll=0;
  for(icoll =0; icoll < ndist; icoll++){
    if(coulomb_type[id][icoll] == '1'){  // For electrons colliding off ions
      coulomb_scatter_ei(adv[id], iran, id, icoll);
    }
    else if(coulomb_type[id][icoll] == '2'){ // For same species collisions
      coulomb_same_species(adv[id], iran, id);
    }
    else if(coulomb_type[id][icoll] == '3'){ // The slow, general algorithm, not implemented yet
      coulomb_general_coll(adv[id], iran, id, icoll);
    }
  }					



  // Boost a random subset of particles to a different energy:
  void boost(int id, particle_dist &adv);
  if (boost_rate[id]>0) boost(id, adv[id]);

  // Let's time the amount of time in this routine:
  vadvance_time += times(&times_buf)-start_time;

}
