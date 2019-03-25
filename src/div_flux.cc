#include <math.h>
#include <stdio.h>
#include "eppic.h"

void div_flux(FArrayND &divFlux,int id,particle_dist &pic,fluid &fspecie, 
	      FTYPE scaler)
{

  enum {XCOMP,YCOMP,ZCOMP};
  int ix,iy,iz;
  divFlux = 0.0;

  void pic_flux(particle_dist &pic, int dim, FArrayND &fluxDim, FTYPE scaler, int n_avg=1);

  static FArrayND flux;
  const int nsize[]={INDICIES(nx+1,ny,nz)}; 
  flux = FArrayND(nsize);


  // -->Are loops necessary for fluid flux, or can
  //    we just do, e.g., flux = fspecie.den*fspecie.vx
  //    if the arrays are the same size?
  
  /* Calculate X flux */
  if (method[id]>=0) {
    pic_flux(pic,XCOMP,flux,1.0);
  } else {
    for (ix=0;ix<nx;ix++) {
      for (iy=0;iy<ny;iy++) {
	for (iz=0;iz<nz;iz++) {
	  flux(INDICIES(ix,iy,iz)) = fspecie.den(INDICIES(ix,iy,iz))*
	    fspecie.vx(INDICIES(ix,iy,iz));
	}
      }
    }
  }

  /* Calculate the derivative w.r.t X (PERIODIC) */
  for (ix=0; ix<nx; ix++)
    for (iy=0; iy<ny; iy++) 
      for (iz=0; iz<nz; iz++) 
	divFlux(INDICIES(ix,iy,iz)) += (flux(INDICIES((ix+1)%nx,iy,iz)) -
					flux(INDICIES((ix+nx-1)%nx,iy,iz)))/(2*dx);

#if NDIM > 1  
  /* Calculate Y flux */
  if (method[id]>=0) {
    pic_flux(pic,YCOMP,flux,1.0);
  } else {
    for (ix=0;ix<nx;ix++) {
      for (iy=0;iy<ny;iy++) {
	for (iz=0;iz<nz;iz++) {
	  flux(INDICIES(ix,iy,iz)) = fspecie.den(INDICIES(ix,iy,iz))*
	    fspecie.vy(INDICIES(ix,iy,iz));
	}
      }
    }
  }

  /* Calculate the derivative w.r.t Y (PERIODIC) */
  for (ix=0; ix<nx; ix++)
    for (iy=0; iy<ny; iy++) 
      for (iz=0; iz<nz; iz++) 
	divFlux(INDICIES(ix,iy,iz)) += (flux(INDICIES(ix,(iy+1)%ny,iz))-
					flux(INDICIES(ix,(iy+ny-1)%ny,iz)))/(2*dy);

#endif
  
#if NDIM == 3

  /* Calculate Z flux */
  if (method[id]>=0) {
    pic_flux(pic,ZCOMP,flux,1.0);
  } else {
    for (ix=0;ix<nx;ix++) {
      for (iy=0;iy<ny;iy++) {
	for (iz=0;iz<nz;iz++) {
	  flux(INDICIES(ix,iy,iz)) = fspecie.den(INDICIES(ix,iy,iz))*
	    fspecie.vz(INDICIES(ix,iy,iz));
	}
      }
    }
  }

  /* Calculate the derivative w.r.t Z (PERIODIC) */
  for (ix=0; ix<nx; ix++)
    for (iy=0; iy<ny; iy++)
      for (iz=0; iz<nz; iz++) 
	divFlux(INDICIES(ix,iy,iz)) += (flux(ix,iy,(iz+1)%nz)-
					flux(ix,iy,(iz+nz-1)%nz))/(2*dz);

#endif

#ifdef USE_DOMAINS
  // if (nsubdomains>1) {
  //   void pass_sum_guard(FArrayND &rho, int nx_local);
  //   pass_sum_guard(divFlux, nx);
  // }
  // void pass_sum_guard(FArrayND &rho, int nx_local);
  // pass_sum_guard(divFlux, nx);
#endif

}
