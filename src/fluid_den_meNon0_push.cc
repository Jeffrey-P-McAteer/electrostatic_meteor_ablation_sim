// WARNING: SUPERCYCLING has been disabled, to make this routine more efficient

// make sure I'm adding the fluid densities correctly for a pure fluid run
// I might be using den_new0 to advance den1, it should be den0

#include <stdio.h>
#include <math.h>
#include "eppic.h"
#include "eppic-efield.h"

// intermediate array; needs to be external for restart purposes
extern fluid *f_old;
fluid *f_old;

// For quasineutral runs, which need qden & divG
// void fluid_den_meNon0_push(fluid *fspecie,field &Efield,particle_dist *pic,
// 			   FArrayND_ranged &qden,FArrayND_ranged &divG,
// 			   int it)
// For quasineutral runs, which need qden & fluxes
void fluid_den_meNon0_push(fluid *fspecie,field &Efield,particle_dist *pic,
			   FArrayND_ranged &qden,
			   INDICIES(FArrayND_ranged &Gx,
				    FArrayND_ranged &Gy,
				    FArrayND_ranged &Gz),
			   int it)
{

  /* This routine uses an Adams-Bashforth method to advance a fluid species 
     assuming it is massless, magnetized, warm, and collisional.  
     It appears as if it applies to more than 1 fluid species,
     but that capability is untested and probably will not work. 

     There are two instabilities which can cause this algorithm to fail.  
     1) The courant condition must be met in all directions
     2) Along B: dt<~coll_rate/w_p^2 where w_p is the plasma freq for the species

     One can push these limits somewhat because the Adams-Bashforth scheme is 
     quasi-implicit but not too much.

  */

  // Initialize some things: 
  void fluid_den_meNon0_step(fluid &fspecie, field &Efield, fluid &f0, fluid &f1, fluid &d2,
			     fluid &f_old, bool store_f, FTYPE dt2, FTYPE dt3);

  // Stored intermediate arrays:
  //  static FArrayND den_p(INDICIES(nx,ny,nz));
  //  fluid *f_old;
  static bool first_entry=true;
  static fluid *tmp;
  
  // Declare Routines to calculate the fluid velocities:
  extern FTYPE fluid_ve_me0_advance(fluid &,field &);
  extern FTYPE fluid_ve_me0_advance_calcE(fluid &,field &);

  // Let's time the amount of time in this routine:
  clock_t start_time=times(&times_buf);

  if (first_entry) {
    // This part initializes the variables
    if (iread == 0) {
      try {
	f_old = new fluid [ndist];
      } catch (const std::bad_alloc &e) {
	std::cout << "Allocation failed (fluid): " << e.what() << "\n";
	MPI_Abort(MPI_COMM_WORLD,-1);
      }
    }
    for (int id=0; id<ndist; ++id) {
      if (method[id] == -2) {
	extern void init_fluid(fluid *&, particle_dist *&);
	init_fluid(tmp, pic); // Will store intermediate density
	f_old[id].den = fspecie[id].den; // Will store intermediate RHS
      }
      if (method[id] == -3) {
	extern void init_fluid(fluid *&, particle_dist *&);
	init_fluid(tmp, pic); // Will store intermediate density
	if (iread == 0) {
	  f_old[id].den = fspecie[id].den; // Will store intermediate RHS
	  f_old[id].vx = fspecie[id].vx;
	  f_old[id].vy = fspecie[id].vy;
	  if (vel_dim[id] == 3) f_old[id].vz = fspecie[id].vz;
	}
      }
    } 
  }

  // Store the PIC charge densities from the previous time step.
  // Do this by subtracting the existing fluid elements rho (passed in).

  FArrayND_ranged den_pic_old=qden;

  for (int id=0; id<ndist; ++id) {
    if (method[id] < 0 && method[id] != -4) {
		
#ifdef USE_DOMAINS
      qden = 0.; // this just forces the qden guard cells to be zero 
      void sub_array(FArrayND_ranged &in,FArrayND &out,
		     INDICIES(int startx, int starty, int startz),
		     INDICIES(int endbeforex, 
			      int endbeforey, 
			      int endbeforez));
      sub_array(fspecie[id].den,qden,
		INDICIES(0,0,0),INDICIES(nx,ny,nz));
      den_pic_old -= qden;
	
#else      
      qden = fspecie[id].den;
      den_pic_old -= qden;
#endif
    }
  }

  // SUPERCYCLING piece:
  /*
  FArrayND den_pic_slope = rho; // to set correct dimensions
  den_pic_slope = 0.; // just to clean up -- this all should be done better
  // Find the charges on the grid at n+1 
  //  charges(rho, pic, fspecie, 0);

  // Remove the current fluid densities (using den_pic_slope as workspace)
  for (int id=0; id<ndist; ++id) {
    if (method[id] < 0 ) {
#ifdef USE_DOMAINS
      void sub_array(FArrayND_ranged &in,FArrayND &out,
		     INDICIES(int startx, int starty, int startz),
		     INDICIES(int endbeforex, 
			      int endbeforey, 
			      int endbeforez));
      sub_array(fspecie[id].den,den_pic_slope,
		INDICIES(0,0,0),INDICIES(nx,ny,nz));
      den_pic_slope *= fspecie[id].q;
      rho -= den_pic_slope;
#else
      den_pic_slope = fspecie[id].den;
      den_pic_slope *= fspecie[id].q;
      rho -= den_pic_slope;
#endif
    }
  }

  // Find the interpolation slope for rho:
  for (int id=0; id<ndist; ++id) {
    if (method[id] < 0 ) {
      den_pic_slope = rho;
      den_pic_slope -= den_pic_old;
      den_pic_slope *= (1./supercycle[id]);
    }
  }
  */

  // Loop through the fluid distributions
  for (int id=0; id<ndist; ++id) {
    if ((efield_algorithm == 2)&&(id == electron_dist)) continue;
    if (method[id] == -2) {

#ifdef FLUID_METHOD_TWO            
      // I have not tested this for domains yet
      FTYPE dt2=0.5*dt/supercycle[id];

      FTYPE vgrid_max=dx;
      if (ndim >= 2) vgrid_max=max(vgrid_max,dy);
      if (ndim >= 3) vgrid_max=max(vgrid_max,dz);
      vgrid_max /= dt/supercycle[id];

      // Supercycle this distribution 
      for (int isuper=0; isuper<supercycle[id]; isuper++) {
	
	// Calculate the electric field for the current fluid densities:
	// not done right because outer cell not interpolated correctly
	for (int ix=0; ix < nx; ix++) {
	  for (int iy=0; iy < ny; iy++) {
	    for (int iz=0; iz < nz; iz++) {
	      rho(INDICIES(ix,iy,iz)) = den_pic_old(INDICIES(ix,iy,iz))
		+(den_pic_slope(INDICIES(ix,iy,iz))*isuper)/supercycle[id]
		+fspecie[id].den(INDICIES(ix,iy,iz))*fspecie[id].q;
	    }
	  }
	}
	
	// Find the electric field on the grid at the n+1 timestep: 
	// efield(Efield, rho);
	// efield_wrapper(Efield, rho, qden, divG, it);
	efield_wrapper(Efield, rho, qden, INDICIES(Gx,Gy,Gz), it);

	// Calculate the updated fluid velocities using the old Efield:
	FTYPE vmax;
	if (Efield.Ex.size() > 0) vmax=fluid_ve_me0_advance(fspecie[id],Efield);
	else vmax=fluid_ve_me0_advance_calcE(fspecie[id],Efield);
	
	if (vmax > 2*vgrid_max) 
	  terminate(-1,"Courant violation: fluid velocities too large");

	if (first_entry && iread == 0) {
	  // On the first time step do a simple (non-predictor-corrector) advance
	  fluid_den_me0_step(fspecie[id], fspecie[id].den, tmp[id].den, 
			     fspecie[id].den, f_old[id].den, 
			     TRUE, dt/supercycle[id], 0.);
	  
	} else { 
	  //On all but the first iteration, do this Adams-Bashforth 2. order :
	  /*       use second order predictor to calculate the right hand side 
		   based on the old value of fluid densities and velocities */
	  
	  fluid_den_me0_step(fspecie[id], fspecie[id].den, tmp[id].den, 
			     fspecie[id].den, f_old[id].den, 
			     TRUE, 3.*dt2, -1.*dt2);
	  
	}
	

	// Temporary output diagnostics:
	/*	particle_misc *misc;
		weight *w;
		extern void output(char *outdir, particle *, fluid *, field &, FArrayND &, 
		particle_misc *misc, weight *, int );
		char outdir = NULL;
		FArrayND rho_bck=rho;
	output (&outdir, pic, tmp, Efield, rho_bck, misc, w, isuper);
	*/
	// Calculate electric field for these new densities:
	// Interpolate to find the PIC densities to this time step (isuper+1):
	for (int ix=0; ix < nx; ix++) {
	  for (int iy=0; iy < ny; iy++) {
	    for (int iz=0; iz < nz; iz++) {
	      rho(INDICIES(ix,iy,iz)) = den_pic_old(INDICIES(ix,iy,iz))
		+(den_pic_slope(INDICIES(ix,iy,iz))*(isuper))/supercycle[id]
		+tmp[id].den(INDICIES(ix,iy,iz))*tmp[id].q;
	    }
	  }
	}
	
	// Find the electric field on the grid at the n+1 timestep: 
	// efield(Efield, rho);
	efield_wrapper(Efield, rho, qden, divG, it);
 
	// Calculate the updated fluid velocities:
	if (Efield.Ex.size() > 0) vmax=fluid_ve_me0_advance(tmp[id],Efield);
	else vmax=fluid_ve_me0_advance_calcE(tmp[id],Efield);
	if (vmax > 2*vgrid_max) 
	  terminate(-1,"Courant violation: fluid velocities too large");
	
	
	// Second order corrector: Adams-Moulton
	
	fluid_den_me0_step(tmp[id], tmp[id].den, fspecie[id].den, 
			   fspecie[id].den, f_old[id].den, FALSE, dt2, dt2);
	
	// Temporary output diagnostics:
	/*	rho_bck=rho;
		output (&outdir, pic, fspecie, Efield, rho_bck, misc, w, isuper);
	*/

      }
#endif
      // METHOD -3: Intertial 
    } else if (method[id] == -3) {

      FTYPE dt2=0.5*dt; // /supercycle[id];
      

      // NO!!!! SUPERCYCLE
      // Supercycle this distribution 
      //      for (int isuper=0; isuper<supercycle[id]; isuper++) {
	
	// Calculate the electric field for the current fluid densities:
	// 1) interpolate the pic values

      // supercycle uses slope
      // no supercycle = no slope
      /*
	for (int ix=0; ix < nx; ix++) {
	  for (int iy=0; iy < ny; iy++) {
	    for (int iz=0; iz < nz; iz++) {
	      rho(INDICIES(ix,iy,iz)) = den_pic_old(INDICIES(ix,iy,iz))
		+(den_pic_slope(INDICIES(ix,iy,iz))*isuper)/supercycle[id];
	    }}}
      */
      qden = den_pic_old;
      // 2) loop over each fluid and add its density to qden
      for (int id2=0; id2 < ndist; id2++) {
	if (method[id2] < 0) {
	  for (int ix=0; ix < nx; ix++) {
	    for (int iy=0; iy < ny; iy++) {
	      for (int iz=0; iz < nz; iz++) {
		qden(INDICIES(ix,iy,iz)) += 
		  fspecie[id2].den(INDICIES(ix,iy,iz));
	      }
	    }
	  }
	}
      }	

      // Find the electric field on the grid at the n+1 timestep: 
      FArrayND dummy_rho;
      // efield_wrapper(Efield, dummy_rho, qden, divG, it);
      efield_wrapper(Efield, dummy_rho, qden, INDICIES(Gx,Gy,Gz), it);

      if (first_entry && iread == 0) {
	  fluid_den_meNon0_step(fspecie[id], Efield, fspecie[id], tmp[id], 
				fspecie[id], f_old[id], 
				TRUE, dt/supercycle[id], 0.);
#ifdef USE_DOMAINS
	  // update guard cells
	  void pass_fluid_guard(fluid&,int);
	  pass_fluid_guard(tmp[id],id);
#endif
 
	} else { 
	  //On all but the first iteration, do this Adams-Bashforth 2. order :
	  /*       use second order predictor to calculate the right hand side 
		   based on the old value of fluid densities and velocities */
	  
	  fluid_den_meNon0_step(fspecie[id], Efield, fspecie[id], tmp[id], 
				fspecie[id], f_old[id], 
				TRUE, 3.*dt2, -1.*dt2);
#ifdef USE_DOMAINS
	  // update guard cells
	  void pass_fluid_guard(fluid&,int);
	  pass_fluid_guard(tmp[id],id);
#endif

      }
	
	// Calculate electric field for these new densities:
	// Interpolate to find the PIC densities to this time step (isuper+1):
	// Calculate the electric field for the current fluid densities:
	// 1) interpolate the pic values

      // again no super, so no slope
      /* 
	for (int ix=0; ix < nx; ix++) {
	  for (int iy=0; iy < ny; iy++) {
	    for (int iz=0; iz < nz; iz++) {
	      rho(INDICIES(ix,iy,iz)) = den_pic_old(INDICIES(ix,iy,iz))
		+(den_pic_slope(INDICIES(ix,iy,iz))*isuper)/supercycle[id];
	    }}}
      */
      qden = den_pic_old;
      // 2) loop over each fluid and add its density to qden
      for (int id2=0; id2 < ndist; id2++) {
	if (method[id2] < 0) {
	  if (id2 == id) {
	    for (int ix=0; ix < nx; ix++) {
	      for (int iy=0; iy < ny; iy++) {
		for (int iz=0; iz < nz; iz++) {
		  qden(INDICIES(ix,iy,iz)) += 
		    tmp[id2].den(INDICIES(ix,iy,iz));
		}
	      }
	    }
	  } else {
	    for (int ix=0; ix < nx; ix++) {
	      for (int iy=0; iy < ny; iy++) {
		for (int iz=0; iz < nz; iz++) {
		  qden(INDICIES(ix,iy,iz)) += 
		    fspecie[id2].den(INDICIES(ix,iy,iz));
		}
	      }
	    }
	  }
	}
      }

      // Find the electric field on the grid at the n+1 timestep: 
      // efield_wrapper(Efield, dummy_rho, qden, divG, it);
      efield_wrapper(Efield, dummy_rho, qden, INDICIES(Gx,Gy,Gz), it);
      
      // Second order corrector: Adams-Moulton
	
      fluid_den_meNon0_step(tmp[id], Efield, tmp[id], fspecie[id], 
			    fspecie[id], f_old[id], FALSE, dt2, dt2);

#ifdef USE_DOMAINS
	// update guard cells
	void pass_fluid_guard(fluid&,int);
	pass_fluid_guard(fspecie[id],id);
#endif
	
	// NO SUPERCYCLE
	//}
    }
  }
  first_entry=false;
  // done with method; reset rho to new value
  // this limits calls to charges() following this routine
  // (in leapadv_subcycle)
  // if we have ions
  // loop over distributions
  // if it's time to update, ie b/c of subcycling 
  // call charges
  
  // 2) loop over each fluid and add its density to rho
  if (method.max() >= 0) { // we have some ions
    for (int id=0; id < ndist; id++) {
      if (subcycle[id] <= 1 || it%subcycle[id] == 0) {
	// quasineutral_den_flux(qden, divG, pic, fspecie, it);
	quasineutral_den_flux(qden, INDICIES(Gx,Gy,Gz), pic, fspecie, it);
	fluid_time += times(&times_buf)-start_time;
	return;
      }
    }
  }
      
  // You get here if you are in between subcycles steps or there are 
  // no particles 
  qden = den_pic_old;
  for (int id=0; id < ndist; id++) {
    if (method[id] < 0) {
      for (int ix=0; ix < nx; ix++) {
	for (int iy=0; iy < ny; iy++) {
	  for (int iz=0; iz < nz; iz++) {
	    qden(INDICIES(ix,iy,iz)) += 
	      fspecie[id].den(INDICIES(ix,iy,iz));
	  }
	}
      }
    }
  }

  // Let's time the amount of time in this routine:
  fluid_time += times(&times_buf)-start_time;

}


// For non-quasineutral runs, which need rho
void fluid_den_meNon0_push(fluid *fspecie,field &Efield,particle_dist *pic,
			   FArrayND &rho,
			   int it)
{

  /* This routine uses an Adams-Bashforth method to advance a fluid species 
     assuming it is massless, magnetized, warm, and collisional.  
     It appears as if it applies to more than 1 fluid species,
     but that capability is untested and probably will not work. 

     There are two instabilities which can cause this algorithm to fail.  
     1) The courant condition must be met in all directions
     2) Along B: dt<~coll_rate/w_p^2 where w_p is the plasma freq for the species

     One can push these limits somewhat because the Adams-Bashforth scheme is 
     quasi-implicit but not too much.

  */

  // Initialize some things: 
  void fluid_den_meNon0_step(fluid &fspecie, field &Efield, fluid &f0, fluid &f1, fluid &d2,
			     fluid &f_old, bool store_f, FTYPE dt2, FTYPE dt3);

  // Stored intermediate arrays:
  //  static FArrayND den_p(INDICIES(nx,ny,nz));
  //  fluid *f_old;
  static bool first_entry=true;
  static fluid *tmp;
  
  // Declare Routines to calculate the fluid velocities:
  extern FTYPE fluid_ve_me0_advance(fluid &,field &);
  extern FTYPE fluid_ve_me0_advance_calcE(fluid &,field &);

  // Let's time the amount of time in this routine:
  clock_t start_time=times(&times_buf);

  if (first_entry) {
    // This part initializes the variables
    if (iread == 0) {
      try {
	f_old = new fluid [ndist];
      } catch (const std::bad_alloc &e) {
	std::cout << "Allocation failed (fluid): " << e.what() << "\n";
      }
    }
    for (int id=0; id<ndist; ++id) {
      if (method[id] == -2) {
	extern void init_fluid(fluid *&, particle_dist *&);
	init_fluid(tmp, pic); // Will store intermediate density
	f_old[id].den = fspecie[id].den; // Will store intermediate RHS
      }
      if (method[id] == -3) {
	extern void init_fluid(fluid *&, particle_dist *&);
	init_fluid(tmp, pic); // Will store intermediate density
	if (iread == 0) {
	  f_old[id].den = fspecie[id].den; // Will store intermediate RHS
	  f_old[id].vx = fspecie[id].vx;
	  f_old[id].vy = fspecie[id].vy;
	  if (vel_dim[id] == 3) f_old[id].vz = fspecie[id].vz;
	}
      }
    } 
  }

  // Store the PIC charge densities from the previous time step.
  // Do this by subtracting the existing fluid elements rho (passed in).

  FArrayND den_pic_old=rho;

  for (int id=0; id<ndist; ++id) {
    if (method[id] < 0 && method[id] != -4) {
		
#ifdef USE_DOMAINS
      rho = 0.; // this just forces the rho guard cell to be zero 
      void sub_array(FArrayND_ranged &in,FArrayND &out,
		     INDICIES(int startx, int starty, int startz),
		     INDICIES(int endbeforex, 
			      int endbeforey, 
			      int endbeforez));
      sub_array(fspecie[id].den,rho,
		INDICIES(0,0,0),INDICIES(nx,ny,nz));
      rho *= fspecie[id].q;
      den_pic_old -= rho;
	
#else      
      rho = fspecie[id].den;
      rho *= fspecie[id].q;
      den_pic_old -= rho;
#endif
    }
  }

  // SUPERCYCLING piece:
  /*
  FArrayND den_pic_slope = rho; // to set correct dimensions
  den_pic_slope = 0.; // just to clean up -- this all should be done better
  // Find the charges on the grid at n+1 
  //  charges(rho, pic, fspecie, 0);

  // Remove the current fluid densities (using den_pic_slope as workspace)
  for (int id=0; id<ndist; ++id) {
    if (method[id] < 0 ) {
#ifdef USE_DOMAINS
      void sub_array(FArrayND_ranged &in,FArrayND &out,
		     INDICIES(int startx, int starty, int startz),
		     INDICIES(int endbeforex, 
			      int endbeforey, 
			      int endbeforez));
      sub_array(fspecie[id].den,den_pic_slope,
		INDICIES(0,0,0),INDICIES(nx,ny,nz));
      den_pic_slope *= fspecie[id].q;
      rho -= den_pic_slope;
#else
      den_pic_slope = fspecie[id].den;
      den_pic_slope *= fspecie[id].q;
      rho -= den_pic_slope;
#endif
    }
  }

  // Find the interpolation slope for rho:
  for (int id=0; id<ndist; ++id) {
    if (method[id] < 0 ) {
      den_pic_slope = rho;
      den_pic_slope -= den_pic_old;
      den_pic_slope *= (1./supercycle[id]);
    }
  }
  */

  // Loop through the fluid distributions
  for (int id=0; id<ndist; ++id) {
    if ((efield_algorithm == 2)&&(id == electron_dist)) continue;
    if (method[id] == -2) {

#ifdef FLUID_METHOD_TWO            
      // I have not tested this for domains yet
      FTYPE dt2=0.5*dt/supercycle[id];

      FTYPE vgrid_max=dx;
      if (ndim >= 2) vgrid_max=max(vgrid_max,dy);
      if (ndim >= 3) vgrid_max=max(vgrid_max,dz);
      vgrid_max /= dt/supercycle[id];

      // Supercycle this distribution 
      for (int isuper=0; isuper<supercycle[id]; isuper++) {
	
	// Calculate the electric field for the current fluid densities:
	// not done right because outer cell not interpolated correctly
	for (int ix=0; ix < nx; ix++) {
	  for (int iy=0; iy < ny; iy++) {
	    for (int iz=0; iz < nz; iz++) {
	      rho(INDICIES(ix,iy,iz)) = den_pic_old(INDICIES(ix,iy,iz))
		+(den_pic_slope(INDICIES(ix,iy,iz))*isuper)/supercycle[id]
		+fspecie[id].den(INDICIES(ix,iy,iz))*fspecie[id].q;
	    }
	  }
	}
	
	// Find the electric field on the grid at the n+1 timestep: 
	// efield(Efield, rho);
	efield_wrapper(Efield, rho, qden, divG, it);

	// Calculate the updated fluid velocities using the old Efield:
	FTYPE vmax;
	if (Efield.Ex.size() > 0) vmax=fluid_ve_me0_advance(fspecie[id],Efield);
	else vmax=fluid_ve_me0_advance_calcE(fspecie[id],Efield);
	
	if (vmax > 2*vgrid_max) 
	  terminate(-1,"Courant violation: fluid velocities too large");

	if (first_entry && iread == 0) {
	  // On the first time step do a simple (non-predictor-corrector) advance
	  fluid_den_me0_step(fspecie[id], fspecie[id].den, tmp[id].den, 
			     fspecie[id].den, f_old[id].den, 
			     TRUE, dt/supercycle[id], 0.);
	  
	} else { 
	  //On all but the first iteration, do this Adams-Bashforth 2. order :
	  /*       use second order predictor to calculate the right hand side 
		   based on the old value of fluid densities and velocities */
	  
	  fluid_den_me0_step(fspecie[id], fspecie[id].den, tmp[id].den, 
			     fspecie[id].den, f_old[id].den, 
			     TRUE, 3.*dt2, -1.*dt2);
	  
	}
	

	// Temporary output diagnostics:
	/*	particle_misc *misc;
		weight *w;
		extern void output(char *outdir, particle *, fluid *, field &, FArrayND &, 
		particle_misc *misc, weight *, int );
		char outdir = NULL;
		FArrayND rho_bck=rho;
	output (&outdir, pic, tmp, Efield, rho_bck, misc, w, isuper);
	*/
	// Calculate electric field for these new densities:
	// Interpolate to find the PIC densities to this time step (isuper+1):
	for (int ix=0; ix < nx; ix++) {
	  for (int iy=0; iy < ny; iy++) {
	    for (int iz=0; iz < nz; iz++) {
	      rho(INDICIES(ix,iy,iz)) = den_pic_old(INDICIES(ix,iy,iz))
		+(den_pic_slope(INDICIES(ix,iy,iz))*(isuper))/supercycle[id]
		+tmp[id].den(INDICIES(ix,iy,iz))*tmp[id].q;
	    }
	  }
	}
	
	// Find the electric field on the grid at the n+1 timestep: 
	// efield(Efield, rho);
	efield_wrapper(Efield, rho, qden, divG, it);
 
	// Calculate the updated fluid velocities:
	if (Efield.Ex.size() > 0) vmax=fluid_ve_me0_advance(tmp[id],Efield);
	else vmax=fluid_ve_me0_advance_calcE(tmp[id],Efield);
	if (vmax > 2*vgrid_max) 
	  terminate(-1,"Courant violation: fluid velocities too large");
	
	
	// Second order corrector: Adams-Moulton
	
	fluid_den_me0_step(tmp[id], tmp[id].den, fspecie[id].den, 
			   fspecie[id].den, f_old[id].den, FALSE, dt2, dt2);
	
	// Temporary output diagnostics:
	/*	rho_bck=rho;
		output (&outdir, pic, fspecie, Efield, rho_bck, misc, w, isuper);
	*/

      }
#endif
      // METHOD -3: Intertial 
    } else if (method[id] == -3) {

      FTYPE dt2=0.5*dt; // /supercycle[id];
      

      // NO!!!! SUPERCYCLE
      // Supercycle this distribution 
      //      for (int isuper=0; isuper<supercycle[id]; isuper++) {
	
	// Calculate the electric field for the current fluid densities:
	// 1) interpolate the pic values

      // supercycle uses slope
      // no supercycle = no slope
      /*
	for (int ix=0; ix < nx; ix++) {
	  for (int iy=0; iy < ny; iy++) {
	    for (int iz=0; iz < nz; iz++) {
	      rho(INDICIES(ix,iy,iz)) = den_pic_old(INDICIES(ix,iy,iz))
		+(den_pic_slope(INDICIES(ix,iy,iz))*isuper)/supercycle[id];
	    }}}
      */
      rho = den_pic_old;
      // 2) loop over each fluid and add its density to rho
      for (int id2=0; id2 < ndist; id2++) {
	if (method[id2] < 0) {
	  for (int ix=0; ix < nx; ix++) {
	    for (int iy=0; iy < ny; iy++) {
	      for (int iz=0; iz < nz; iz++) {
		rho(INDICIES(ix,iy,iz)) += 
		  fspecie[id2].den(INDICIES(ix,iy,iz))*fspecie[id2].q;
	      }
	    }
	  }
	}
      }

      // Find the electric field on the grid at the n+1 timestep: 
      // FArrayND_ranged dummy_qden,dummy_divG;
      // efield_wrapper(Efield, rho, dummy_qden, dummy_divG, it);
      FArrayND_ranged dummy_qden,dummy_Gx,dummy_Gy,dummy_Gz;
      efield_wrapper(Efield, rho, dummy_qden, 
		     INDICIES(dummy_Gx,
			      dummy_Gy,
			      dummy_Gz),it);

      if (first_entry && iread == 0) {
	  fluid_den_meNon0_step(fspecie[id], Efield, fspecie[id], tmp[id], 
			     fspecie[id], f_old[id], 
			     TRUE, dt/supercycle[id], 0.);
#ifdef USE_DOMAINS
	  // update guard cells
	  void pass_fluid_guard(fluid&,int);
	  pass_fluid_guard(tmp[id],id);
#endif
 
	} else { 
	  //On all but the first iteration, do this Adams-Bashforth 2. order :
	  /*       use second order predictor to calculate the right hand side 
		   based on the old value of fluid densities and velocities */
	  
	  fluid_den_meNon0_step(fspecie[id], Efield, fspecie[id], tmp[id], 
			     fspecie[id], f_old[id], 
			     TRUE, 3.*dt2, -1.*dt2);
#ifdef USE_DOMAINS
	  // update guard cells
	  void pass_fluid_guard(fluid&,int);
	  pass_fluid_guard(tmp[id],id);
#endif

      }
	
	// Calculate electric field for these new densities:
	// Interpolate to find the PIC densities to this time step (isuper+1):
	// Calculate the electric field for the current fluid densities:
	// 1) interpolate the pic values

      // again no super, so no slope
      /* 
	for (int ix=0; ix < nx; ix++) {
	  for (int iy=0; iy < ny; iy++) {
	    for (int iz=0; iz < nz; iz++) {
	      rho(INDICIES(ix,iy,iz)) = den_pic_old(INDICIES(ix,iy,iz))
		+(den_pic_slope(INDICIES(ix,iy,iz))*isuper)/supercycle[id];
	    }}}
      */
      rho = den_pic_old;
      // 2) loop over each fluid and add its density to rho
      for (int id2=0; id2 < ndist; id2++) {
	if (method[id2] < 0) {
	  if (id2 == id) {
	    for (int ix=0; ix < nx; ix++) {
	      for (int iy=0; iy < ny; iy++) {
		for (int iz=0; iz < nz; iz++) {
		  rho(INDICIES(ix,iy,iz)) += 
		    tmp[id2].den(INDICIES(ix,iy,iz))*fspecie[id2].q;
		}
	      }
	    }
	  } else {
	    for (int ix=0; ix < nx; ix++) {
	      for (int iy=0; iy < ny; iy++) {
		for (int iz=0; iz < nz; iz++) {
		  rho(INDICIES(ix,iy,iz)) += 
		    fspecie[id2].den(INDICIES(ix,iy,iz))*fspecie[id2].q;
		}
	      }
	    }
	  }
	}
      }

      // Find the electric field on the grid at the n+1 timestep: 
      // efield_wrapper(Efield, rho, dummy_qden, dummy_divG, it);
      efield_wrapper(Efield, rho, dummy_qden, 
		     INDICIES(dummy_Gx,
			      dummy_Gy,
			      dummy_Gz),it);

	// Second order corrector: Adams-Moulton
	
	fluid_den_meNon0_step(tmp[id], Efield, tmp[id], fspecie[id], 
			   fspecie[id], f_old[id], FALSE, dt2, dt2);

#ifdef USE_DOMAINS
	// update guard cells
	void pass_fluid_guard(fluid&,int);
	pass_fluid_guard(fspecie[id],id);
#endif
	
	// NO SUPERCYCLE
	//}
    }
  }
  first_entry=false;
  // done with method; reset rho to new value
  // this limits calls to charges() following this routine
  // (in leapadv_subcycle)
  // if we have ions
  // loop over distributions
  // if it's time to update, ie b/c of subcycling 
  // call charges
  
  // 2) loop over each fluid and add its density to rho
  if (method.max() >= 0) { // we have some ions
    for (int id=0; id < ndist; id++) {
      if (subcycle[id] <= 1 || it%subcycle[id] == 0) {
	charges(rho, pic, fspecie, it);
	fluid_time += times(&times_buf)-start_time;
	return;
      }
    }
  }
      
  // You get here if you are in between subcycles steps or there are 
  // no particles 
  rho = den_pic_old;
  for (int id=0; id < ndist; id++) {
    if (method[id] < 0) {
      for (int ix=0; ix < nx; ix++) {
	for (int iy=0; iy < ny; iy++) {
	  for (int iz=0; iz < nz; iz++) {
	    rho(INDICIES(ix,iy,iz)) += 
	      fspecie[id].den(INDICIES(ix,iy,iz))*fspecie[id].q;
	  }
	}
      }
    }
  }

  // Let's time the amount of time in this routine:
  fluid_time += times(&times_buf)-start_time;

}

