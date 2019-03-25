/* This subroutine calculates the divergence of the current, useful for 
   determining the envelope field in a two time-scale approximation */

// still periodic, needs to be updated for domains

#include "eppic.h"
#include "eppic-mpi.h"

void divj(FArrayND &divj, particle_dist *pic, fluid *fspecie)
{

  int id, ix, iy, iz;
  FArrayND j=divj;

  void wgather(FArrayND &den, INDICIES(PTYPEAVec &x, PTYPEAVec &y, PTYPEAVec &z), 
	       PTYPEAVec &wt, FTYPE n0);

  // Clear output array:
  divj=0;
  
  /* Calculate Jx (current in the x direction) */
  for (id=0; id<ndist; ++id) {

    /* The standard PIC gather: */
    if (method[id] >= 0) {
      wgather(j, INDICIES(pic[id].x, pic[id].y, pic[id].z), pic[id].vx, 
	      pic[id].n0avg*pic[id].q*dx/dt);
      
    } else {//Fluid approach:
#ifdef USE_DOMAINS
      /*  Does not work yet.. not very important so just commented it out
      void sub_array(FArrayND_ranged &in,FArrayND &out,
		     INDICIES(int startx, int starty, int startz),
		     INDICIES(int endbeforex, int endbeforey, int endbeforez));
      void sub_array_mult(FArrayND_ranged &in,FArrayND &out,
			  INDICIES(int startx, int starty, int startz),
			  INDICIES(int endbeforex, int endbeforey, int endbeforez));
      sub_array(fspecie[id].vx,j,
		INDICIES(0,0,0),
		INDICIES(nx,ny,nz)
		);
      sub_array_mult(fspecie[id].den,j,
		INDICIES(0,0,0),
		INDICIES(nx,ny,nz)
		);
      j *= fspecie[id].q;
      */
#else
      j = fspecie[id].vx;
      j *= fspecie[id].den;
      j *= fspecie[id].q;
#endif
    }

    /* Calculate the derivative of Jx w.r.t. x */
    for (ix=0; ix<nx; ix++)
      for (iy=0; iy<ny; iy++) 
	for (iz=0; iz<nz; iz++) 
	  divj(INDICIES(ix,iy,iz)) += ( j(INDICIES((ix+1)%nx,iy,iz)) -
					j(INDICIES((ix+nx-1)%nx,iy,iz)) )/(2*dx);

  }
#if NDIM > 1  
  for (id=0; id<ndist; ++id) {
    /* Calculate Jy (current in the y direction) */

    /* The standard weighted PIC gather: */
    if (method[id] >= 0) {
      wgather(j, INDICIES(pic[id].x, pic[id].y, pic[id].z), pic[id].vy, 
	      pic[id].n0avg*pic[id].q*dx/dt);

    } else {//Fluid approach:
#ifdef USE_DOMAINS
      void sub_array(FArrayND_ranged &in,FArrayND &out,
		     INDICIES(int startx, int starty, int startz),
		     INDICIES(int endbeforex, int endbeforey, int endbeforez));
      void sub_array_mult(FArrayND_ranged &in,FArrayND &out,
			  INDICIES(int startx, int starty, int startz),
			  INDICIES(int endbeforex, int endbeforey, int endbeforez));

      sub_array(fspecie[id].vy,j,
		INDICIES(0,0,0),
		INDICIES(nx,ny,nz)
		);
      sub_array_mult(fspecie[id].den,j,
		INDICIES(0,0,0),
		INDICIES(nx,ny,nz)
		);
      j *= fspecie[id].q;
#else
      j = fspecie[id].vy;
      j *= fspecie[id].den;
      j *= fspecie[id].q;
#endif
    }
    
    /* Calculate the derivative of Jy w.r.t. y */
    for (ix=0; ix<nx; ix++)
      for (iy=0; iy<ny; iy++) 
	for (iz=0; iz<nz; iz++) 
	  divj(INDICIES(ix,iy,iz)) += (j(INDICIES(ix,(iy+1)%ny,iz))-
				       j(INDICIES(ix,(iy+ny-1)%ny,iz)) )/(2*dy);
  }
#endif
  
#if NDIM == 3
  for (id=0; id<ndist; ++id) {
    /* Calculate Jz (current in the z direction) */
    /* The standard weighted PIC gather: */
    if (method[id] >= 0) {
      wgather(j, pic[id].x, pic[id].y, pic[id].z, pic[id].vz, 
	      pic[id].n0avg*pic[id].q*dx/dt);
    } else {//Fluid approach:
#ifdef USE_DOMAINS
      void sub_array(FArrayND_ranged &in,FArrayND &out,
		     INDICIES(int startx, int starty, int startz),
		     INDICIES(int endbeforex, int endbeforey, int endbeforez));
      void sub_array_mult(FArrayND_ranged &in,FArrayND &out,
			  INDICIES(int startx, int starty, int startz),
			  INDICIES(int endbeforex, int endbeforey, int endbeforez));

      sub_array(fspecie[id].vz,j,
		INDICIES(0,0,0),
		INDICIES(nx,ny,nz)
		);
      sub_array_mult(fspecie[id].den,j,
		INDICIES(0,0,0),
		INDICIES(nx,ny,nz)
		);
      j *= fspecie[id].q;
#else
      j = fspecie[id].vz;
      j *= fspecie[id].den;
      j *= fspecie[id].q;
#endif
    }
    /* Calculate the derivative of Jz w.r.t. z */
    for (ix=0; ix<nx; ix++)
      for (iy=0; iy<ny; iy++) {
	for (iz=0; iz<nz; iz++) 
	  divj(INDICIES(ix,iy,iz)) += (j(ix,iy,(iz+1)%nz)-
				       j(ix,iy,(iz+nz-1)%nz))/(2*dz);
      }
  }
#endif

  // Collect from all processors:  
#ifdef USE_MPI
  j=0;
  int mpi_err=MPI_Reduce((void*) &(divj(INDICIES(0,0,0))),
			 (void*) &(j(INDICIES(0,0,0))), divj.size(),
			 MPI_FTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if ( mpi_err != MPI_SUCCESS) 
    mpi_error(mpi_rank, mpi_err,"Reduce divj call in divj");
  j /= mpi_np;
  divj = j;
#endif    

}

	
