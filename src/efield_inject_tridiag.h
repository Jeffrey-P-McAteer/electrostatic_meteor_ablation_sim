

#ifndef EPPIC_INJECT_TRIDIAG_H
#define EPPIC_INJECT_TRIDIAG_H

#include "eppic-efield.h" /* for definition of field, other types */


#include "realfftwmpi_many.h"
typedef realfftwmpi_many<FTYPE,NDIM-1> rfftw_many;

// indicies for data, real array
#if NDIM==2
#define RFFTW_INDICIES(x,y,z) y, x
#elif NDIM==3
#define RFFTW_INDICIES(x,y,z) y, z, x
#endif


// indicies for cdata, ie complex array
#if NDIM==2
#define RFFTW_C_INDICIES(x,y,z) y, x
#elif NDIM==3
#define RFFTW_C_INDICIES(x,y,z) z, y, x
#endif

#include "eppic_fftw.h"

// -------------------------------------------------------------------------
// Main routine: efield_inject_tridiag
// ------------------------------------------------------------------------- 

void efield_inject_tridiag(field &Efield , FArrayND &rho);


// -------------------------------------------------------------------------
// helper routines used in the main routine: 
// -------------------------------------------------------------------------

void efield_inject_tridiag_init(FArrayND &rho,
				field &Efield,
				rfftw_many &working_array);

void efield_inject_tridiag_solve(rfftw_many &working_array,
				  rfftw_many &boundary_working_array,
				  field &Efield);

/// _solve_batch is being replaced with a series of routines for each case

void efield_inject_tridiag_solve_periodic(rfftw_many &working_array,
				  rfftw_many &boundary_working_array,
				  field &Efield);

void efield_inject_tridiag_solve_open(rfftw_many &working_array,
				  rfftw_many &boundary_working_array,
				  field &Efield);

void  efield_inject_tridiag_solve_mixed(rfftw_many &working_array,
					rfftw_many &boundary_working_array,
					field &Efield);

void  efield_inject_tridiag_solve_dirichlet(rfftw_many &working_array,
					    rfftw_many &boundary_working_array,
					    field &Efield);

void efield_inject_tridiag_solve_batch(rfftw_many &working_array,
				  rfftw_many &boundary_working_array,
				  field &Efield);

void efield_inject_tridiag_sync_subdomain(rfftw_many &working_array,
					  rfftw_many &boundary_working_array,
					  field &Efield);

void efield_inject_tridiag_sync_subdomain_batch(rfftw_many &working_array,
					  rfftw_many &boundary_working_array,
					  field &Efield);

void efield_inject_tridiag_tophi(rfftw_many &working_array,
				 rfftw_many &boundary_working_array,
				 field &Efield);

inline int isolve_to_iz(int isolve, 
			int local_ny, int starty, 
			int local_nz, int startz) {
  return int(isolve/2/local_ny)+startz;
};

inline int isolve_to_iy(int isolve, 
			int local_ny, int starty, 
			int local_nz, int startz) {
  return isolve/2 - int(isolve/2/local_ny)*local_ny +starty;
};
  
#endif
