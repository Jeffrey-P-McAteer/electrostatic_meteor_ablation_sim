/* Calculate quasineutral density and flux for a given species.

   This rountine is designed to calculate quantities as efficiently 
   as possible, at the expense of temporary array storage. */

#if USE_QN
#include <math.h>
#include <stdio.h>
#include "eppic.h"

void qden_flux(FArrayND &qden,
	       INDICIES(FArrayND &Gx,FArrayND &Gy,FArrayND &Gz),
	       int id,
	       particle_dist *pic,fluid *fspecie,FTYPE scaler=1.0)
{

  /* Local function declarations */
  void gather_den_flux(FArrayND &den, 
  		       INDICIES(FArrayND &xflux,FArrayND &yflux,FArrayND &zflux), 
  		       INDICIES(PTYPEAVec &x,PTYPEAVec &y,PTYPEAVec &z),
  		       INDICIES(PTYPEAVec &xwt,PTYPEAVec &ywt,PTYPEAVec &zwt),
		       const int np,
  		       const int nxmin);
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);
  void pass_sum_guard(FArrayND &, int nx);
  void write_local_bin(FArrayND_ranged array, const char *name, int *nx_ghost);
  void write_local_bin(FArrayND array, const char *name);


  /* Local variables */
  int ix,iy,iz;
  FTYPE den_scale = scaler*pic[id].n0avg*(nx*nsubdomains)*ny*nz/
    (npd[id]*nsubdomains);
  int nxmin = 0;
  if (boundary_type[0] == INJECT) nxmin = -1;

  /* Gather PIC density and flux, or copy them from fluid */
  if (method[id] == 0) { // PIC
    gather_den_flux(qden,
		    INDICIES(Gx,Gy,Gz),
		    INDICIES(pic[id].x,pic[id].y,pic[id].z),
    		    INDICIES(pic[id].vx,pic[id].vy,pic[id].vz),
		    pic[id].np,
    		    nxmin);
    qden *= den_scale;
    Gx *= den_scale*dx/dt;
#if NDIM > 1
    Gy *= den_scale*dy/dt;
#if NDIM > 2
    Gz *= den_scale*dz/dt;
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
	  // Only index den array once
	  den_value = fspecie[id].den(INDICIES(ix,iy,iz));
	  qden(INDICIES(ix,iy,iz)) = den_value;
	  Gx(INDICIES(ix,iy,iz)) = den_value*
	    fspecie[id].vx(INDICIES(ix,iy,iz));
#if NDIM > 1
	  Gy(INDICIES(ix,iy,iz)) = den_value*
	    fspecie[id].vy(INDICIES(ix,iy,iz));
#if NDIM > 2
	  Gz(INDICIES(ix,iy,iz)) = den_value*
	    fspecie[id].vz(INDICIES(ix,iy,iz));	  
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

  if (hybrid_diagnostic_subcycle > 0 && 
      it%nout*hybrid_diagnostic_subcycle == 0) {
    char wlb_name[32];
    sprintf(wlb_name,"qden_tmp-%06d-",it);
    write_local_bin(qden,wlb_name);
    sprintf(wlb_name,"Gx_tmp-%06d-",it);
    write_local_bin(Gx,wlb_name);
    sprintf(wlb_name,"Gy_tmp-%06d-",it);
    write_local_bin(Gy,wlb_name);
  }
  
}
#endif
