/* Advance the particle velocities, positions and the fields 
 to the next time steps allowing different time-steps for different species.
 This routine uses a 2nd order Leap-frog method
 Other methods can be substitted here. */

#include "eppic.h"
#include "eppic-mpi.h"
#include "eppic-efield.h"
#include "eppic-io.h"

// void leapadv_subcycle(particle_dist *pic, fluid *fspecie, field &Efield, 
// 		      FArrayND &rho, FArrayND_ranged &qden, 
// 		      FArrayND_ranged &divG,
// 		      particle_misc *misc, int it)
void leapadv_subcycle(particle_dist *pic, fluid *fspecie, field &Efield, 
		      FArrayND &rho, FArrayND_ranged &qden, 
		      INDICIES(FArrayND_ranged &Gx,
			       FArrayND_ranged &Gy,
			       FArrayND_ranged &Gz),
		      particle_misc *misc, int it)
{

  
  extern void vadvance(int, particle_dist *, particle_dist *, particle_dist *,
		        particle_misc *, fluid *, field &, FTYPE);
  extern void xadvance(int, particle_dist *, particle_dist *, particle_dist *,
                       particle_misc *, fluid *, FTYPE, FArrayND &);

  static int first_entry=TRUE;
  static field *Estore;
  static intAVec subcycle_repeat;
  static int iss=0;
  static FILE *fphistore[MAXDIST]={0};
  

  if (first_entry == TRUE) {
    first_entry=FALSE;
    /* For each subcycling (slowly cycling) distribution create a field 
       storage array */
    Estore=new field[ndist];
    subcycle_repeat=intAVec(ndist) = 0;
    
    for (int id=0; id<ndist; ++id) {
      if (subcycle[id] > 1) {
	/* check to see if this is the same as a previous subcycle rate */
	int repeat=0;
	for (int id2=id-1; id2>=0; --id2) 
	  if (subcycle[id] == subcycle[id2]) repeat=id2;
	if (repeat == 0) {
	  Estore[id]=Efield;
	  subcycle_repeat[id]=id;
	} else {
	  subcycle_repeat[id]=repeat;
	}
      }
    }
    // We will output Estore from this routine rather than output:
    FILE* openarrayfile(char* fileroot, int it, int subcycle);
    char nameroot[32];
    if (subdomain.rank == 0)
      for (int id=0; id<ndist; ++id) {
	if (phistore_out_subcycle[id] > 0 && subcycle[id] > 1) {
	  sprintf(nameroot,"phistore%d",id);
	  fphistore[id]=openarrayfile(nameroot,it-1,phistore_out_subcycle[id]);
	} 
      }
  }

  /* Advance the subcycle counter */
  iss++;
    
  /* Advance base cycle particle velocities to the next 1/2 timestep 
     from the previous 1/2 timestep: */
  for (int id=0; id<ndist; ++id) {
    if (method[id] >= 0)
      if (iss%subcycle[id]==0) 
	if (subcycle[id]<=1) {
	  vadvance(id,pic,pic,pic,misc,fspecie,Efield,(FTYPE)1.*subcycle[id]);
	} else { 
	  int id_store=subcycle_repeat[id];
	  vadvance(id,pic,pic,pic,misc,fspecie,Estore[id_store],
		   (FTYPE)1.*subcycle[id]);
	}
  }

  /* Push the particle positions to the next timestep: */
  for (int id=0; id<ndist; ++id) {
    if (method[id] >= 0)
      if (iss%subcycle[id]==0) {
	xadvance(id,pic,pic,pic,misc,fspecie,(FTYPE)(1. * subcycle[id]), rho);
      }
  } /* END   for (id=0; id<ndist; ++id) */
  
  // Fluid advance of distribution. This will also update particles
  // if there are any. Therefore, any particle-related updates need
  // to go in this routine, too.
  if (method.min() < 0) {
    if (efield_algorithm == 2) {
      // extern void fluid_den_meNon0_push(fluid *fspecie, field &Efield, 
      // 					particle_dist *pic, 
      // 					FArrayND_ranged &qden, 
      // 					FArrayND_ranged &divG,
      // 					int it);
      // fluid_den_meNon0_push(fspecie, Efield, pic, qden, divG, iss);
      extern void fluid_den_meNon0_push(fluid *fspecie, field &Efield, 
					particle_dist *pic, 
					FArrayND_ranged &qden, 
					INDICIES(FArrayND_ranged &Gx,
						 FArrayND_ranged &Gy,
						 FArrayND_ranged &Gz),
					int it);
      fluid_den_meNon0_push(fspecie, Efield, pic, qden, 
			    INDICIES(Gx,Gy,Gz), iss);
      // charges are easily calculated in fluid_den_meNon0_push
    } else {
      extern void fluid_den_meNon0_push(fluid *fspecie, field &Efield, 
					particle_dist *pic, 
					FArrayND &rho,
					int it);
      fluid_den_meNon0_push(fspecie, Efield, pic, rho, iss);
      // charges are easily calculated in fluid_den_meNon0_push
    }
  } else {
    /* Find the charges and currents on the grid at n+1 */
    if (efield_algorithm == 2) {
      // quasineutral_den_flux(qden, divG, pic, fspecie, it);
      quasineutral_den_flux(qden, INDICIES(Gx,Gy,Gz), pic, fspecie, it);
    }
    else if (efield_algorithm != 8) {
      charges(rho, pic, fspecie, it);
    }
  }


  // Added by Glenn 2018/01/18 To test twostream
  // Get mean of rho and subtract it to maintain neutrality
  /*  float rho_mean = 0;
  for (int i = 0; i<rho.size(); i++){
    rho_mean += rho(i);
  }
  rho_mean /= rho.size();
  rho -= rho_mean;
  */

  //  if (method.min() < 0) {
  //    // Adams Bashforth advance of inertialless fluid :
  //    extern void fluid_den_me0_push(fluid *fspecie, field &Efield, 
  //				   particle_dist *pic, FArrayND &rho);
  //    fluid_den_me0_push(fspecie, Efield, pic, rho);

  //  }

  //  Find the electric field on the grid at t=n:
  // efield_wrapper(Efield, rho, qden, divG, it);
  efield_wrapper(Efield, rho, qden, INDICIES(Gx,Gy,Gz), it);

  /* Store the appropriate E field values in Estore */
  for (int id=0; id<ndist; ++id) 
    if (subcycle[id] > 1 && subcycle_repeat[id] == id) {
      /* Zero the Estore arrays after using them */
      if (iss%subcycle[id]==0) {
	Estore[id].phi_rho = 0;
	/* Not using E* in this code
	Estore[id].Ex = 0.;
  	Estore[id].Ey = 0.; 
	Estore[id].Ez = 0.;
	*/
      }
      /* Sum Efield arrays and store in Estore */
      Estore[id].phi_rho += Efield.phi_rho;
      /* not using E* in this code
      Estore[id].Ex += Efield.Ex;
      Estore[id].Ey += Efield.Ey;
      Estore[id].Ez += Efield.Ez;
      */
      /* Convert the summed Efield arrays into an average Efield array*/
      if ((iss+1)%subcycle[id]==0) {
	Estore[id].phi_rho /= (FTYPE) subcycle[id];
	/* Not using E* in this code
	Estore[id].Ex /= (FTYPE) subcycle[id];
	Estore[id].Ey /= (FTYPE) subcycle[id];
	Estore[id].Ez /= (FTYPE) subcycle[id];
	*/

	// We will perform output of the Estore field data here because doing it
	// in output require extensive modifications
	int ngrid[] = {INDICIES(nx,ny,nz)};
	if ((iss+1)%(nout*phistore_out_subcycle[id])==0) 
	  if (subdomain.rank == 0) output_array(fphistore[id], Estore[id].phi_rho,ngrid,nout_avg); 

      }
    }

  return;

  
} /* advance */
