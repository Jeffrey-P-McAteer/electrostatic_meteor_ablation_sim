#include <stdio.h>
#include <math.h>
#include <sys/times.h> //Clock related functions
#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif

#ifdef USE_DOMAINS
#include "classes/realfftwmpi_transpose.h"
typedef realfftwmpi<FTYPE,NDIM> FArrayND_fftw;
#else 
#include <complex>
#include "classes/realfftw2.h"
typedef realfftw2<FTYPE,NDIM> FArrayND_fftw;
#endif

void periodic_filter(FArrayND &array)
{
  /* This subroutine filters array using the fsteep fwidth variables, */
  clock_t start_time=times(&times_buf);


  if (fsteep != 0 && fwidth != 0) { 


    /* Local variables */
    
    const int nsize[]={INDICIES(nx*nsubdomains,ny,nz)};
#ifdef USE_DOMAINS
    static FArrayND_fftw array_trans(subdomain.neighbor_comm, nsize);
#else
    static FArrayND_fftw array_trans(nsize);
#endif
    
    
    // Copy phi into phi_trans, setting all non-overlapping values to zero
    array_trans.data=0.0;
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) {
	  array_trans(INDICIES(ix,iy,iz))=array(INDICIES(ix,iy,iz));
	}

    // Transform into k space
    array_trans.transform();

    int nx_global=nx*nsubdomains;
    FTYPE ksqr, damp, kx, ky, kz; 
    FTYPE kxbase=fwidth/(nx_global/2.), 
      kybase=fwidth/(ny/2.), kzbase=fwidth/(nz/2.);
    
    for (int ikx=0; ikx<nx_global; ikx++) {
      int ikx_global=ikx;
      if (ikx_global <= nx_global/2) kx=kxbase*ikx_global;
      else kx=kxbase*(-(nx_global-ikx_global));
      
      for (int iky=0; iky<array_trans.ny_transpose; iky++) {
	int iky_global=iky+array_trans.y_start_transpose;
	if (iky_global <= ny/2) ky=kybase*iky_global;
	else ky=kybase*(-(ny-iky_global));
	for (int ikz=0; ikz<nz; ikz++) {
	  if (ikz <= nz/2) kz=kzbase*ikz;
	  else kz=kzbase*(-(nz-ikz));
	  ksqr=kx*kx+ky*ky+kz*kz;
	  damp=exp(-1*pow(ksqr,fsteep));
	  if (ndim == 1) {
	    if (ikx_global <= nx_global/2) {
	      array_trans.cdata(INDICIES(0,ikx,0)) *= damp;
	      }
	  }
	  else if (ndim == 2) {
	    if (iky_global <= ny/2) {
	      array_trans.cdata(INDICIES(iky,ikx,0)) *= damp;
	      }
	  } else {
	    if (ikz <= nz/2) {
		array_trans.cdata(INDICIES(iky,ikx,ikz)) *= damp;
	    }
	  }
	}
      }
    }
    

    // Inverse fft transform:
    array_trans.invtransform();

    //Copy the rest of the phi_trans array to phi
    for (int ix = phix_guard_size[1]; ix < nx-phix_guard_size[0]; ++ix)
      for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) {
	  array(INDICIES(ix,iy,iz))=array_trans(INDICIES(ix,iy,iz));
	}

      
  } /* end if fsteep,fwidth!=0 */
}

void periodic_filter(FArrayND_ranged &array,const int nx_guard[])
{
  /* This subroutine filters array using the fsteep fwidth variables, */
  clock_t start_time=times(&times_buf);


  if (fsteep != 0 && fwidth != 0) { 


    /* Local variables */
    
    const int nsize[]={INDICIES(nx*nsubdomains,ny,nz)};
#ifdef USE_DOMAINS
    static FArrayND_fftw array_trans(subdomain.neighbor_comm, nsize);
#else
    static FArrayND_fftw array_trans(nsize);
#endif
    
    
    // Copy phi into phi_trans, setting all non-overlapping values to zero
    array_trans.data=0.0;
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) {
	  array_trans(INDICIES(ix,iy,iz))=array(INDICIES(ix,iy,iz));
	}

    // Transform into k space
    array_trans.transform();

    int nx_global=nx*nsubdomains;
    FTYPE ksqr, damp, kx, ky, kz; 
    FTYPE kxbase=fwidth/(nx_global/2.), 
      kybase=fwidth/(ny/2.), kzbase=fwidth/(nz/2.);
    
    for (int ikx=0; ikx<nx_global; ikx++) {
      int ikx_global=ikx;
      if (ikx_global <= nx_global/2) kx=kxbase*ikx_global;
      else kx=kxbase*(-(nx_global-ikx_global));
      
      for (int iky=0; iky<array_trans.ny_transpose; iky++) {
	int iky_global=iky+array_trans.y_start_transpose;
	if (iky_global <= ny/2) ky=kybase*iky_global;
	else ky=kybase*(-(ny-iky_global));
	for (int ikz=0; ikz<nz; ikz++) {
	  if (ikz <= nz/2) kz=kzbase*ikz;
	  else kz=kzbase*(-(nz-ikz));
	  ksqr=kx*kx+ky*ky+kz*kz;
	  damp=exp(-1*pow(ksqr,fsteep));
	  if (ndim == 1) {
	    if (ikx_global <= nx_global/2) {
	      array_trans.cdata(INDICIES(0,ikx,0)) *= damp;
	      }
	  }
	  else if (ndim == 2) {
	    if (iky_global <= ny/2) {
	      array_trans.cdata(INDICIES(iky,ikx,0)) *= damp;
	      }
	  } else {
	    if (ikz <= nz/2) {
		array_trans.cdata(INDICIES(iky,ikx,ikz)) *= damp;
	    }
	  }
	}
      }
    }
    

    // Inverse fft transform:
    array_trans.invtransform();

    //Copy the rest of the phi_trans array to phi
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) {
	  array(INDICIES(ix,iy,iz))=array_trans(INDICIES(ix,iy,iz));
	}

  } /* end if fsteep,fwidth!=0 */
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);
  pass_guards(array,nx_guard);

}
