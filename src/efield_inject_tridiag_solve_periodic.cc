

/// Sets up the main equations for the many tridiagonal solves, calls the 
/// batch solver, and stores the solutions back in the working arrays.
/// 
/// This is set up for both 2D an 3D cases, and for both cases where there
/// are multiple processors per domain. The equations and boundary condition
/// definitions are taken from Birdsall. 
///
///  \param rfftw_many working_array -- solutions on the grid
///  \param rfftw_many boundary_working_array -- solutions at the boundary
///  \param field      Efield -- used for frequnetly used values
///

#include "efield_inject_tridiag.h" /* for field, types and fftw stuff */
#include "eppic-mpi.h" /* for subdomain stuff */
#include "eppic-math.h" /* for ceil function */
#include "tridiag.h"
#include <iostream>

void  efield_inject_tridiag_solve_periodic(rfftw_many &working_array,
				  rfftw_many &boundary_working_array,
				  field &Efield)
{

  // initialize dimensions
  int nx = Efield.nx;
  int ny = Efield.ny;
  int nz = Efield.nz;

  // divide workload in slowest varying direction
  // in 2D, this is y, in 3D this is z
  // division along slowest leads to most efficient MPI communications
  // we do real and imaginary component in pairs 
  bool handlesKy0Eqn = false;
  int nNz=nz/2+1;  
#if NDIM==2  
  int nNy=ny/2+1;
  int dny = nNy/subdomain.np;
  int ystart = dny*subdomain.rank;
  int yend = ystart+dny;
  // this next line is really for the right most processor
  // basically dny is calculated by rounding down, so in odd cases, the last
  // processor might need to do more equations.
  if (yend +dny > nNy) yend = nNy; 
  if (ystart == 0) handlesKy0Eqn = true;
  int zstart = 0;
  int zend = 1;
  //'*2' applies to real and imag terms
  int n_solve = (yend-ystart)*nNz*2; 
#else
  int nNy=ny;
  int ystart = 0;
  int yend = ny;
  int zstart = working_array.y_start_transpose;
  int zend = zstart+working_array.ny_transpose;
  if (zstart==0) handlesKy0Eqn = true;
  //'*2' applies to real and imag terms
  int n_solve = (zend-zstart)*nNy*2; 
#endif

  int n_eqn = n_solve;
  if ((n_eqn > MAX_TRIDIAG_SOLVE) && MAX_TRIDIAG_SOLVE > 0) {
    n_eqn = MAX_TRIDIAG_SOLVE;
  }

  // first entry, includes:
  // set up b(y,z) terms and rinv terms (for open cases)
  static bool first_entry=true;
  static ArrayNd<FTYPE,2> byz,rinv;

  if(first_entry) {
    first_entry=false;

    index_type constants_size[] = {nNy,nNz};
    byz = ArrayNd<FTYPE,2>(constants_size);
    rinv = ArrayNd<FTYPE,2>(constants_size);

    FTYPE dx = Efield.dx;
    FTYPE dy = Efield.dy;
    FTYPE dz = Efield.dz;

    for(int iy=0;iy<nNy;iy++) 
      for(int iz=0;iz<nNz;iz++) {
      FTYPE dm=1.
	+2.*Sqr((dx/dy)*sin(M_PI* (FTYPE)iy /ny))
	+2.*Sqr((dx/dz)*sin(M_PI* (FTYPE)iz /nz));
      FTYPE r=dm+sqrt(Sqr(dm)-1.);
      rinv(iy)=1./r;// used only for open cases
      byz(iy,iz)=-2*dm;
      }

  }

   
 // setup values common to all boundary conditions
  if (n_solve == 0) 
    return;
  ArrayNd<FTYPE,2> a,b,c,lhs,rhs;
  a=ArrayNd<FTYPE,2>(n_eqn,nx)= 1;
  b=ArrayNd<FTYPE,2>(n_eqn,nx);
  c=ArrayNd<FTYPE,2>(n_eqn,nx) = 1;
  ArrayNd<FTYPE,1> alpha,beta;
  alpha = ArrayNd<FTYPE,1>(n_eqn) = 1;
  rhs=ArrayNd<FTYPE,2>(n_eqn,nx);
  lhs=ArrayNd<FTYPE,2>(n_eqn,nx);
  FTYPE factor = -Efield.dx*Efield.dx/Efield.eps;
  int local_nNy=yend-ystart;
  int local_nNz=zend-zstart;
  int isolve=0;

  // some where about here I need to loop over sub nsolve
  for (int curr_eqn=0; curr_eqn < n_solve; curr_eqn+=n_eqn) {
    
    if (Efield.boundary[0]==periodic) {
      beta = alpha = 1;
      a=1;
      b=1;
      c=1;
      if (handlesKy0Eqn && (curr_eqn == 0)) {
	// even though alpha,beta apply to boundaries, all domains solve
	// a global problem involving these terms. So if this processor solves
	// for cdata(0).real and cdata(1).imag (ie isolve = 0 and 1) then
	// they need the alpha and beta terms adjusted (not cyclic in this case)
	alpha(0) = 0;
	alpha(1) = 0;
	beta(0) = 0;
	beta(1) = 0;
      }
    }

    for (isolve=0; isolve<n_eqn; isolve+=2) {
      int iz=isolve_to_iz(isolve+curr_eqn, local_nNy, ystart, local_nNz, zstart);
      int iy=isolve_to_iy(isolve+curr_eqn, local_nNy, ystart, local_nNz, zstart);
      // this handles the cases where n_solve is not a multiple of n_eqn
      // basically the last equation is repeatedly solved to runtime errors
      if (iz>=zend) {iz=zend-1;};
      if (iy>=yend) {iy=yend-1;};

      for (int ix=0;ix<nx;ix++) {
	rhs(isolve,ix) = 
	  working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)).real()*factor;
	b(isolve,ix) = byz(iy,iz);
      }
      for (int ix=0;ix<nx;ix++) {
	rhs(isolve+1,ix) = 
	  working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)).imag()*factor;
	b(isolve+1,ix) = byz(iy,iz); // same for both real,imag
      }
    }

    if (curr_eqn == 0) {
      // set up left side boundaries
      if (subdomain.id_number==0) {
	// nothing to do for non kyz==0 terms
	// kyz=0 term is special
	isolve=0;
	if (handlesKy0Eqn) {
	  // first term
	  a(isolve,0) = 0;
	  b(isolve,0) = 1;
	  c(isolve,0) = 0;
	  rhs(isolve,0) = 0;
	  alpha(isolve) = 0;
	  beta(isolve) = 0;
	  isolve++;
	  a(isolve,0) = 0;
	  b(isolve,0) = 1;
	  c(isolve,0) = 0;
	  rhs(isolve,0) = 0;
	  alpha(isolve) = 0;
	  beta(isolve) = 0;
	  // second term;
	  // nothing special for this term, done correctly above!
	}
      }
      // set up right side boundaries
      if (subdomain.id_number==nsubdomains-1) {
	// nothing to do for non kyz==0 terms
	// kyz=0 term is not special
	isolve=0;
	if (handlesKy0Eqn) {
	  // first term
	  alpha(isolve) = 0;
	  beta(isolve) = 0;
	  isolve++;
	  alpha(isolve) = 0;
	  beta(isolve) = 0;
	  // second term;
	  // nothing special for this term, done correctly above!
	}
      }
    }

    // batch solve
    batch_tridiag_cyclic_solver_parallel(a,b,c,rhs,lhs,alpha,beta,
					 subdomain.neighbor_comm);
    
    // record solutions
    for (isolve=0; isolve<n_eqn; isolve+=2) {
      int iz=isolve_to_iz(isolve+curr_eqn, local_nNy, ystart, local_nNz, zstart);
      int iy=isolve_to_iy(isolve+curr_eqn, local_nNy, ystart, local_nNz, zstart);
      if (iz>=zend) {continue;};
      if (iy>=yend) {continue;};
      for (int ix=0;ix<nx;ix++) {
	working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)) = 
	  std::complex<FTYPE>(lhs(isolve,ix),lhs(isolve+1,ix));
      }
    }
  }
    
}
