/* Calculate quasineutral density and flux divergence for a given
   species.

   This rountine is designed to calculate quantities as efficiently 
   as possible, at the expense of temporary array storage. */

#if USE_QN
#include <math.h>
#include <stdio.h>
#include "eppic.h"

void qden_divG(FArrayND &qden,FArrayND &divG,int id,
	       particle_dist *pic,fluid *fspecie,FTYPE scaler=1.0)
{

  /* Local function declarations */
  void gather_den_flux(FArrayND &den, 
  		       INDICIES(FArrayND &xflux,FArrayND &yflux,FArrayND &zflux), 
  		       INDICIES(PTYPEAVec &x,PTYPEAVec &y,PTYPEAVec &z),
  		       INDICIES(PTYPEAVec &xwt,PTYPEAVec &ywt,PTYPEAVec &zwt),
		       const int np,
  		       const int nxmin);
  void write_local_bin(FArrayND_ranged array, const char *name, int *nx_ghost);
  void write_local_bin(FArrayND array, const char *name);


  /* Zero in/out arrays */
  // qden = 0.0;
  // divG = 0.0;

  /* Local variables */
  int ix,iy,iz;
  int tmp_size[]={INDICIES(nx+xguard_size,ny,nz)};
  int flux_size[]={INDICIES(nx+qnx_guard_size[0]+qnx_guard_size[1],ny,nz)};
  int flux_start[]={INDICIES(0-qnx_guard_size[0],0,0)};
  FArrayND xflux_tmp;
  FArrayND_ranged xflux;
  xflux_tmp = FArrayND(tmp_size);
  xflux_tmp = 0.0;
  xflux = FArrayND_ranged(flux_start,flux_size);
#if NDIM > 1
  FArrayND yflux_tmp;
  FArrayND_ranged yflux;
  yflux_tmp = FArrayND(tmp_size);
  yflux_tmp = 0.0;
  yflux = FArrayND_ranged(flux_start,flux_size);
#if NDIM > 2
  FArrayND zflux_tmp;
  FArrayND_ranged zflux;
  zflux_tmp = FArrayND(tmp_size);
  zflux_tmp = 0.0;
  zflux = FArrayND_ranged(flux_start,flux_size);
#endif
#endif

  FTYPE den_scale = scaler*pic[id].n0avg*(nx*nsubdomains)*ny*nz/
    (npd[id]*nsubdomains);
  int nxmin = 0;
  if (boundary_type[0] == INJECT) nxmin = -1;

  /* Gather PIC density and flux, or copy them from fluid */
  if (method[id] == 0) { // PIC
    gather_den_flux(qden,
		    INDICIES(xflux_tmp,yflux_tmp,zflux_tmp),
		    INDICIES(pic[id].x,pic[id].y,pic[id].z),
    		    INDICIES(pic[id].vx,pic[id].vy,pic[id].vz),
		    pic[id].np,
    		    nxmin);
    qden *= den_scale;
    xflux_tmp *= den_scale*dx/dt;
#if NDIM > 1
    yflux_tmp *= den_scale*dy/dt;
#if NDIM > 2
    zflux_tmp *= den_scale*dz/dt;
#endif
#endif

    // OPEN: UNTESTED
    if ((boundary_type[0] == OPEN)) {
      if (subdomain.id_number==0) {
    	for (iy=0; iy<ny; iy++) 
    	  for (iz=0; iz<nz; iz++) 
    	    qden(INDICIES(0,iy,iz)) += scaler*n0lhsd[id]/2;
      } 
      if (subdomain.id_number==nsubdomains-1) {
    	for (iy=0; iy<ny; iy++) 
    	  for (iz=0; iz<nz; iz++) 
    	    qden(INDICIES(nx,iy,iz)) += scaler*n0rhsd[id]/2;
      }
    }

    // INJECT: there are no particles < 0, so density at 0 is half what it should be
    if ((boundary_type[0] == INJECT)) {
      if (subdomain.id_number==0 && nsubdomains > 1) {
    	for (iy=0; iy<ny; iy++) 
    	  for (iz=0; iz<nz; iz++) 
    	    qden(INDICIES(0,iy,iz)) += scaler*n0lhsd[id]/2;
      } 
      if (subdomain.id_number==nsubdomains-1) {
    	for (iy=0; iy<ny; iy++) 
    	  for (iz=0; iz<nz; iz++) 
    	    qden(INDICIES(nx,iy,iz)) += scaler*n0rhsd[id]/2;
      }
    }
    // if ((boundary_type[0] == INJECT)) {
    //   if (subdomain.id_number==0) {
    // 	for (iy=0; iy<ny; iy++) 
    // 	  for (iz=0; iz<nz; iz++) 
    // 	    qden(INDICIES(0,iy,iz))+=n0lhsd[id]/2;
    //   } 
    //   if (subdomain.id_number==nsubdomains-1) {
    // 	for (iy=0; iy<ny; iy++) 
    // 	  for (iz=0; iz<nz; iz++) 
    // 	    qden(INDICIES(nx,iy,iz))+=n0rhsd[id]/2;
    //   }
    // }

  } else if (method[id] == -2) { // inertialess fluid
    FTYPE den_value=0.0;
    for (ix=0;ix<nx;ix++) {
      for (iy=0;iy<ny;iy++) {
	for (iz=0;iz<nz;iz++) {
	  den_value = fspecie[id].den(INDICIES(ix,iy,iz)); // Only index array once
	  qden(INDICIES(ix,iy,iz)) = den_value;
	  xflux_tmp(INDICIES(ix,iy,iz)) = den_value*fspecie[id].vx(INDICIES(ix,iy,iz));
#if NDIM > 1
	  yflux_tmp(INDICIES(ix,iy,iz)) = den_value*fspecie[id].vy(INDICIES(ix,iy,iz));
#if NDIM > 2
	  zflux_tmp(INDICIES(ix,iy,iz)) = den_value*fspecie[id].vz(INDICIES(ix,iy,iz));	  
#endif
#endif
	}
      }
    }
  } else if (method[id] == -3) { // inertial fluid (TO DO: see density.cc)
    char message[130];
    sprintf(message,"Method=%d not implemented",method[id]);
    terminate(-1,message);
  } else if (method[id] == -4) { // static, inertialess, quasineutral fluid
  } else {
    char message[130];
    sprintf(message,"Method=%d not implemented",method[id]);
    terminate(-1,message);
  }

#ifdef USE_DOMAINS
  void pass_sum_guard(FArrayND &, int nx);
  pass_sum_guard(xflux_tmp,nx);
#if NDIM > 1
  pass_sum_guard(yflux_tmp,nx);
#if NDIM > 2
  pass_sum_guard(zflux_tmp,nx);
#endif
#endif
#endif

  // write_local_bin(qden,"qden_tmp");
  // write_local_bin(xflux_tmp,"xflux_tmp");
  // write_local_bin(yflux_tmp,"yflux_tmp");
  
  /* Copy values from non-ghosted arrays into ghosted arrays */
  xflux = 0.0;
#if NDIM > 1
  yflux = 0.0;
#if NDIM > 2
  zflux = 0.0;
#endif
#endif
  for (ix=0;ix<nx+xguard_size;ix++) {
    for (iy=0;iy<ny;iy++) {
      for (iz=0;iz<nz;iz++) {
  	xflux(INDICIES(ix,iy,iz)) = xflux_tmp(INDICIES(ix,iy,iz));
#if NDIM > 1
  	yflux(INDICIES(ix,iy,iz)) = yflux_tmp(INDICIES(ix,iy,iz));
#if NDIM > 2
	zflux(INDICIES(ix,iy,iz)) = zflux_tmp(INDICIES(ix,iy,iz));
#endif
#endif
      }
    }
  }  

  /* Pass flux guard cell values around before taking derivative */
#ifdef USE_DOMAINS
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);
  pass_guards(xflux,qnx_guard_size);
#if NDIM > 1
  pass_guards(yflux,qnx_guard_size);
#if NDIM > 2
  pass_guards(zflux,qnx_guard_size);
#endif
#endif
#endif

  // write_local_bin(xflux,"xflux");
  // write_local_bin(yflux,"yflux");
    
  /* Calculate local Div[flux] */
  FTYPE dxInv=1.0,dyInv=1.0,dzInv=1.0;
  switch(ndim)
    {
    case 3: dzInv = 1.0/dz;
    case 2: dyInv = 1.0/dy;
    case 1: dxInv = 1.0/dx;
    }
  FTYPE divGx=0.0,divGy=0.0,divGz=0.0;
  for (ix=0;ix<nx;ix++) {
    for (iy=0;iy<ny;iy++) {
      for (iz=0;iz<nz;iz++) {
 	divGx = 0.5*dxInv*(xflux(INDICIES(ix+1,iy,iz))
			   -xflux(INDICIES(ix-1,iy,iz)));
#if NDIM > 1
	divGy = 0.5*dyInv*(yflux(INDICIES(ix,(iy+1)%ny,iz))
			   -yflux(INDICIES(ix,(iy+ny-1)%ny,iz)));
#if NDIM > 2
	divGz = 0.5*dzInv*(zflux(INDICIES(ix,iy,(iz+1)%nz))
			   -zflux(INDICIES(ix,iy,(iz+nz-1)%nz)));
#endif
#endif
	divG(INDICIES(ix,iy,iz)) += divGx+divGy+divGz;
      }
    }
  }

}
#endif
