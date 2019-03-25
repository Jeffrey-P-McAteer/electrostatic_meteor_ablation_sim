

#include "efield_inject_tridiag.h" /* for field, types */
#include "eppic-mpi.h" /* for subdomain stuff */
#include "eppic-math.h" /* for ceil function */

/** Passes nx by ny_subdomain array to other procs in subdomain and assigns to
    correct portion of rho. 
    
    This is used if there are multiple processors per subdomain. In this case, 
    the efield_inject_parallel routines divide the work calculating each iy 
    row amongst the subdomain processors, and this info needs to be 
    interchanged before the calculation finishes.

    \param rho a FArrayND to be modified (input and output)
*/

inline void set_real(std::complex<FTYPE> &rhs, FTYPE lhs) {
  rhs = std::complex<FTYPE>(lhs, rhs.imag());
}

inline void set_imag(std::complex<FTYPE> &rhs, FTYPE lhs) {
  rhs = std::complex<FTYPE>(rhs.real(), lhs);
}


void  efield_inject_tridiag_sync_subdomain(rfftw_many &working_array,
					   rfftw_many &boundary_working_array,
					   field &Efield)
{

  if (subdomain.np==1) return;
  MPI_Barrier(MPI_COMM_WORLD);
  int nx = Efield.nx;
  int ny = Efield.ny;
  int nz = Efield.nz;

#if NDIM==2
  int dny = static_cast<int>(ceil(static_cast<FTYPE>((ny/2)/subdomain.np)));
  int nz_end=nz;

#else
  int dny = ny/subdomain.np;
  int nz_end=nz/2+1;
#endif

  int ny_rho_tmp = 2*dny;

  // add extra columns for phi-1 and phiNx, ie rho_boundary 
  // MAY NOT WORK FOR 1 PROC!!
  int nx_extra=0;
  if (subdomain.id_number==0) nx_extra++;
  if (subdomain.id_number==nsubdomains-1) {
    nx_extra++;
  }
  FArrayND rho_tmp = FArrayND(INDICIES(nx+nx_extra,ny_rho_tmp,nz_end));

  for(int isubdomain=0;isubdomain<subdomain.np;isubdomain++){
  int iy_start = dny*isubdomain;
  int iy_end = iy_start+dny;
#if NDIM==2
  if (iy_end>ny/2) iy_end = ny/2;
#else
  if (iy_end>ny) iy_end = ny;
#endif
    if (subdomain.rank == isubdomain) {
      rho_tmp=0;
      if (iy_start==0) {//special!
	for (int iz=0;iz<nz_end;iz++)  
	  for (int ix=0;ix<nx;ix++) {
	    rho_tmp(INDICIES(ix,0,iz)) = 
	      working_array.cdata(RFFTW_C_INDICIES(ix,0,iz)).real();
#if NDIM==2
	    rho_tmp(INDICIES(ix,ny_rho_tmp-1,iz)) = 
	      working_array.cdata(RFFTW_C_INDICIES(ix,ny/2,iz)).real();
#else
	    rho_tmp(INDICIES(ix,ny_rho_tmp-1,iz)) = 
	      working_array.cdata(RFFTW_C_INDICIES(ix,0,nz/2)).real();
#endif
	  }
	if (subdomain.id_number==0) {
	  for (int iz=0;iz<nz_end;iz++) {
	    rho_tmp(INDICIES(nx,0,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(0,0,iz)).real();
#if NDIM==2
	    rho_tmp(INDICIES(nx,ny_rho_tmp-1,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(0,ny/2,iz)).real();
#else
	    rho_tmp(INDICIES(nx,ny_rho_tmp-1,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(0,0,nz/2)).real();
#endif
	  }
	}
	if (subdomain.id_number==nsubdomains-1) {
	  int xpos=nx;
	  if (subdomain.id_number==0) xpos++;
	  for (int iz=0;iz<nz_end;iz++) {
	    rho_tmp(INDICIES(xpos,0,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,iz)).real();
#if NDIM==2
	    rho_tmp(INDICIES(xpos,ny_rho_tmp-1,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(1,ny/2,iz)).real();
#else
	    rho_tmp(INDICIES(xpos,ny_rho_tmp-1,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,nz/2)).real();
#endif
	  }
	  for (int iz=0;iz<nz_end;iz++) {
	    rho_tmp(INDICIES(xpos,0,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,iz)).real();
#if NDIM==2
	    rho_tmp(INDICIES(xpos,ny_rho_tmp-1,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(1,ny/2,iz)).real();
#else
	    rho_tmp(INDICIES(xpos,ny_rho_tmp-1,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,nz/2)).real();
#endif
	  }

	}
	//iy_start++;
      }
      for (int iz=0;iz<nz_end;iz++) 
	for (int iy=iy_start;iy<iy_end;iy++) {
	  if (iy==0) continue;
	  for (int ix=0;ix<nx;ix++) {
	    
	    rho_tmp(INDICIES(ix,iy-iy_start,iz)) = 
	      working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)).real();
	    rho_tmp(INDICIES(ix,ny_rho_tmp-(iy-iy_start)-1,iz)) = 
	      working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)).imag();
	  }
	}
      if (subdomain.id_number==0) {
	for (int iz=0;iz<nz_end;iz++) 
	  for (int iy=iy_start;iy<iy_end;iy++) {
	    if (iy==0) continue;
	    rho_tmp(INDICIES(nx,iy-iy_start,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(0,iy,iz)).real();
	    rho_tmp(INDICIES(nx,ny_rho_tmp-(iy-iy_start)-1,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(0,iy,iz)).imag();
	  }
      }
      if (subdomain.id_number==nsubdomains-1) {
	int xpos=nx-1+nx_extra;
	for (int iz=0;iz<nz_end;iz++) 
	  for (int iy=iy_start;iy<iy_end;iy++) {
	    if (iy==0) continue;
	    rho_tmp(INDICIES(xpos,iy-iy_start,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(1,iy,iz)).real();
	    rho_tmp(INDICIES(xpos,ny_rho_tmp-(iy-iy_start)-1,iz)) = 
	      boundary_working_array.cdata(RFFTW_C_INDICIES(1,iy,iz)).imag();
	  }
      }
    }
    MPI_Bcast(
	      &(rho_tmp(INDICIES(0,0,0))),
	      rho_tmp.size(), 
	      MPI_FTYPE,
	      isubdomain,
	      subdomain.internal_comm);
    if (subdomain.rank != isubdomain) {
      if (iy_start==0) {// special
	for (int iz=0;iz<nz_end;iz++)  
	  for (int ix=0;ix<nx;ix++) {
	    set_real(working_array.cdata(RFFTW_C_INDICIES(ix,0,iz)), 
	      rho_tmp(INDICIES(ix,0,iz)));
#if NDIM==2
	    set_real(working_array.cdata(RFFTW_C_INDICIES(ix,ny/2,iz)), 
	      rho_tmp(INDICIES(ix,ny_rho_tmp-1,iz)));
#else
	    set_real(working_array.cdata(RFFTW_C_INDICIES(ix,0,nz/2)), 
	      rho_tmp(INDICIES(ix,ny_rho_tmp-1,iz)));
#endif
	  }
	if (subdomain.id_number==0) {
	  for (int iz=0;iz<nz_end;iz++) {
	    set_real(boundary_working_array.cdata(RFFTW_C_INDICIES(0,0,iz)), 
	      rho_tmp(INDICIES(nx,0,iz)));
#if NDIM==2
	    set_real(boundary_working_array.cdata(RFFTW_C_INDICIES(0,ny/2,iz)), 
	      rho_tmp(INDICIES(nx,ny_rho_tmp-1,iz)));
#else
	    set_real(boundary_working_array.cdata(RFFTW_C_INDICIES(0,0,nz/2)), 
	      rho_tmp(INDICIES(nx,ny_rho_tmp-1,iz)));
#endif
	  }
	}
	if (subdomain.id_number==nsubdomains-1) {
	  int xpos=nx-1+nx_extra;
	  for (int iz=0;iz<nz_end;iz++) {
	    set_real(boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,iz)), 
	      rho_tmp(INDICIES(xpos,0,iz)));
#if NDIM==2
	    set_imag(boundary_working_array.cdata(RFFTW_C_INDICIES(1,ny/2,iz)), 
	      rho_tmp(INDICIES(xpos,ny_rho_tmp-1,iz)));
#else
	    set_imag(boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,nz/2)), 
	      rho_tmp(INDICIES(xpos,ny_rho_tmp-1,iz)));
#endif
	  }
	}
	
	//iy_start++;
      }
      for (int iz=0;iz<nz_end;iz++) 
	for (int iy=iy_start;iy<iy_end;iy++) {
	  if (iy==0) continue;
	  for (int ix=0;ix<nx;ix++) {
	    set_real(working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)), 
	      rho_tmp(INDICIES(ix,iy-iy_start,iz)));
	    set_imag(working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)), 
	      rho_tmp(INDICIES(ix,ny_rho_tmp-(iy-iy_start)-1,iz)));
	  }
	}
      if (subdomain.id_number==0) {
	for (int iz=0;iz<nz_end;iz++) 
	  for (int iy=iy_start;iy<iy_end;iy++) {
	  if (iy==0) continue;

	    set_real(boundary_working_array.cdata(RFFTW_C_INDICIES(0,iy,iz)), 
	      rho_tmp(INDICIES(nx,iy-iy_start,iz)));
	    set_imag(boundary_working_array.cdata(RFFTW_C_INDICIES(0,iy,iz)), 
	      rho_tmp(INDICIES(nx,ny_rho_tmp-(iy-iy_start)-1,iz)));
	  }
      }
      if (subdomain.id_number==nsubdomains-1) {
	int xpos=nx-1+nx_extra;
	for (int iz=0;iz<nz_end;iz++) 
	  for (int iy=iy_start;iy<iy_end;iy++) {
	    if (iy==0) continue;
	    set_real(boundary_working_array.cdata(RFFTW_C_INDICIES(1,iy,iz)), 
	      rho_tmp(INDICIES(xpos,iy-iy_start,iz)));
	    set_imag(boundary_working_array.cdata(RFFTW_C_INDICIES(1,iy,iz)), 
	      rho_tmp(INDICIES(xpos,ny_rho_tmp-(iy-iy_start)-1,iz)));
	  }
      }
    }
  }
}
