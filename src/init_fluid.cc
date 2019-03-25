/* Initialize the fluid variables */
// this nees to be updated to cleanly deal with different CPP NDIM defs. 
#include "eppic.h"
#include "eppic-mpi.h"
#include <math.h>

#include "init_func_graddrift.h"

void init_fluid(fluid *&fspecie, particle_dist *&pic)
{
  if (method.min() >= 0) return; /* no fluid distributions */
  
  /* We have ndist distributions, not all are fluid */
  /* Create a set of fluid arrays, one for each distribution */
  try {
    fspecie = new fluid [ndist];
  } catch (const std::bad_alloc &e) {
    std::cout << "Allocation failed (fluid): " << e.what() << "\n";
    MPI_Abort(MPI_COMM_WORLD,-1);
  }

  // Change the random number seed for each processor 
  static int iran = -33 - subdomain.id_number;
  FTYPE ran3(int *);

  // use enum to give different initializations memorable names
  enum {
    random = -2,
    gaussian = -3,
    sine_guass = -4,
    graddrift = -5,
    constant = -6,
    copy = 0,
    evolve = 10
  };


  /* Fluid distribution setup */
  for (int id=0; id<ndist; ++id) {
    // if (method[id] < 0) {
    if (method[id] < 0 && method[id] != -4) {
      fspecie[id].m=md[id];
      fspecie[id].q=qd[id];
      fspecie[id].n0=pic[id].n0avg;
      fspecie[id].nu=coll_rate[id];
      fspecie[id].gamma=thermal_gamma[id];
      fspecie[id].T0=fspecie[id].m*(Sqr(vxthd[id])+Sqr(vythd[id]))/2.;
      fspecie[id].diffc=diffc[id];

#ifdef USE_DOMAINS
      /* create space for the density array and set to n0*/
      const int nsize_den[] =
	{INDICIES(nx+2*denx_guard_size,ny,nz)};
      int ndenstart[] =
	{INDICIES(-1*denx_guard_size,0,0)};
      fspecie[id].den = FArrayND_ranged(ndenstart, nsize_den);
      fspecie[id].den = fspecie[id].n0;

      /* create space for the velocity arrays and set v0*/
      const int nsize_vel[] = 
	{INDICIES(nx+2*velx_guard_size,ny,nz)};
      int nvelstart[] = 
	{INDICIES(-1*velx_guard_size,0,0)};
      fspecie[id].vx = FArrayND_ranged(nvelstart, nsize_vel);
      fspecie[id].vx = vx0d[id];
      if (vel_dim[id] >= 2) { 
	fspecie[id].vy = FArrayND_ranged(nvelstart, nsize_vel);
	fspecie[id].vy = vy0d[id];
      }
      if (vel_dim[id]==3) { 
	fspecie[id].vz = FArrayND_ranged(nvelstart, nsize_vel);
	fspecie[id].vz = vz0d[id];
      }
#else      
      /* create space for the density array and set to n0*/
      fspecie[id].den = FArrayND(INDICIES(nx, ny, nz));
      fspecie[id].den = fspecie[id].n0;
      
      /* create space for the velocity arrays and set v0*/
      fspecie[id].vx = FArrayND(INDICIES(nx, ny, nz)) = vx0d[id];
      if (vel_dim[id] >= 2) 
	fspecie[id].vy = FArrayND(INDICIES(nx, ny, nz)) = vy0d[id];
      if (vel_dim[id]==3) 
	fspecie[id].vz = FArrayND(INDICIES(nx, ny, nz)) = vz0d[id];
#endif    
      if (init_dist[id] == constant) { // Set the density equal to 
	// the constant value n0
	// printf("Using constant density for dist %d\n",id);
	for (int ix=0;
	     ix < (nsubdomains*nx);
	     ix++) {
	  if ((ix >= subdomain.id_number*nx) &&
	      (ix < (subdomain.id_number*nx+nx))) {
	    for (int iy=0; iy < ny; iy++) {
	      for (int iz=0; iz < nz; iz++) {
		fspecie[id].den(INDICIES(ix-subdomain.id_number*nx,iy,iz)) +=
		  fspecie[id].n0;
	      }
	    }
	  }
	}
      }
  
      if (init_dist[id] == random) { // Set the density equal to
	// a random noise plus a background
	// printf("Using random density for dist %d\n",id);
	for (int ix=0; 
	     ix < (nsubdomains*nx); 
	     ix++) {
	  if ((ix >= subdomain.id_number*nx) && 
	      (ix < (subdomain.id_number*nx+nx))) {
	    for (int iy=0; iy < ny; iy++) {
	      for (int iz=0; iz < nz; iz++) {
		fspecie[id].den(INDICIES(ix-subdomain.id_number*nx,iy,iz)) += 
		  fspecie[id].n0*param1[id]*(ran3(&iran)-0.5);
	      }
	    }
	  }
	}
      }

      if (init_dist[id] == gaussian) { // Set the density equal to 
	// a gaussian plus a background
	// printf("Using gaussian density for dist %d\n",id);
	for (int ix=0; 
	     ix < (nsubdomains*nx); 
	     ix++) {
	  if ((ix >= subdomain.id_number*nx) && 
	      (ix < (subdomain.id_number*nx+nx))) {
	    for (int iy=0; iy < ny; iy++) {
	      for (int iz=0; iz < nz; iz++) {
		fspecie[id].den(INDICIES(ix-subdomain.id_number*nx,iy,iz)) += 
		  param1[id]*fspecie[id].n0*
		  exp(-Sqr((ix-nx*nsubdomains/2.)*param2[id]))*
		  exp(-Sqr((iy-ny/2.)*param3[id]));
	      }
	    }
	  }
	}
      }

      if (init_dist[id] == sine_guass) { // Set the density equal to
	// a sinx wave plus a background
	// printf("Using sine_guass for density for dist %d\n",id);
	FTYPE kxval = 2*PI*int(param2[id])/(nx*nsubdomains);
	FTYPE kyval = 2*PI*int(param3[id])/(ny);
#if NDIM == 3
	FTYPE kzval = 2*PI*int(param4[id])/(nz);
#endif
	for (int ix=0; 
	     ix < (nsubdomains*nx); 
	     ix++) {
	  if ((ix >= subdomain.id_number*nx) && 
	      (ix <  (subdomain.id_number*nx+nx))) {
	    for (int iy=0; iy < ny; iy++) {
#if NDIM == 3
	      for (int iz=0; iz < nz; iz++) {
		fspecie[id].den(INDICIES(ix-subdomain.id_number*nx,iy,iz)) += 
		  fspecie[id].n0*param1[id]*
		  cos(kxval*ix)*
		  cos(kyval*iy)*
		  cos(kzval*iz);
#elif NDIM < 3
		fspecie[id].den(INDICIES(ix-subdomain.id_number*nx,iy,iz)) += 
		  fspecie[id].n0*param1[id]*
		  cos(kxval*ix)*
		  cos(kyval*iy);
#endif
#if NDIM == 3
	      }
#endif
	    }
	  }
	}
      }

      if (init_dist[id] == graddrift) {
	// printf("Using graddrift for density for dist %d\n",id);
	for(int ix=0;ix<nx;ix++) {
	  for( int iy=0; iy<ny; iy++) {
	    for( int iz=0; iz<nz; iz++) {
	      fspecie[id].den(INDICIES(ix,iy,iz))=
		init_func_graddrift(id,
				    INDICIES(PTYPE(ix),PTYPE(iy),PTYPE(iz)))
		*fspecie[id].n0;
	    }
	  }
	}
      }


      if (init_dist[id] >= copy) {
	// printf("Using copy or evolve for density for dist %d\n",id);
	int copy_dist = init_dist[id];
	if (init_dist[id] >= evolve) copy_dist -= 10;
	if (copy_dist >= id && method[init_dist[id]] < 0) {
	  char string[130];
	  sprintf(string,"%s%d%s%d%s\n","Error: init_dist for distribution ", 
		  id," is >= 0, but > ",id,
		  ",ie it has not been created yet!\n");
	  terminate(-1,string);
	}
	if (method[copy_dist] < 0) { // fluid
	  fspecie[id].den = fspecie[copy_dist].den;
	} else {
	  // Set the density equal to the density distribution of dist[id]
	  void density(FArrayND &, int, particle_dist *, fluid *, FTYPE);
	  FArrayND den,dentmp;
	  const int nsize_den[]={INDICIES(nx+xguard_size,ny,nz)};
	  den =  FArrayND(nsize_den);
	  dentmp =  FArrayND(nsize_den);
	  density(dentmp, copy_dist, pic, fspecie, 
		  //		fspecie[id].n0/pic[init_dist[id]].n0);
		  1.0);
#ifdef USE_MPI
	  // Sum the arrays around the networked processes 
	  void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);
	
	  // Use density as the target array
	  den = 0;
	  int mpi_err=MPI_Allreduce((void*) &(dentmp(INDICIES(0,0,0))), 
				    (void*) &(den(INDICIES(0,0,0))), 
				    dentmp.size(),
				    MPI_FTYPE, MPI_SUM,
				    subdomain.internal_comm);
	  if ( mpi_err != MPI_SUCCESS) 
	    mpi_error(mpi_rank, mpi_err,"Allreduce call in init_fluid");
	  den /= subdomain.np;
	
#ifdef USE_DOMAINS
	  void pass_sum_guard(FArrayND &, int nx);
	  pass_sum_guard(den, nx);
#endif
#endif
	  for (int ix=0; ix<nx; ix++) 
	    {
#if NDIM >=2
	      for (int iy=0; iy<ny; iy++) 
		{
#endif
#if NDIM == 3
		  for (int iz=0; iz<nz; iz++)
		    {
#endif
		      fspecie[id].den(INDICIES(ix,iy,iz)) = 
			den(INDICIES(ix,iy,iz));
#if NDIM >=2
		    }
#endif
#if NDIM == 3
		}
#endif
	    }
	  fspecie[id].n0 = pic[copy_dist].n0avg;
	  //	  n0peak[id] = fspecie[id].n0;
	  if (mpi_rank == 0) {
	    // printf("; WARNING: N0 for fluid adjusted!\n");
	    printf("\tn0avgd%1d = %g\n",id,fspecie[id].n0);
	  }
	    
	} // else pic
	
      } // if copy

#ifdef USE_DOMAINS
      // initialize guard cells
      void pass_fluid_guard(fluid&,int);
      pass_fluid_guard(fspecie[id],id);
#endif

      if (init_dist[id] >= evolve) {
	// take one step

	// initialize test

	// loop while test unsatisfied
      }
    }
  } /* END for (id=0; id<ndist; ++id) */
  
}


