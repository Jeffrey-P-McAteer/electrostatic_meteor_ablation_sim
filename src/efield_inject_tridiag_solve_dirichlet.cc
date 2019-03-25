

/// Sets up the main equations for the many tridiagonal solves, calls the 
/// batch solver, and stores the solutions back in the working arrays.
/// 
/// This is set up for both 2D an 3D cases, and for both cases where there
/// are multiple processors per domain. The equations are taken from Birdsall. 
/// And the boundary conditions are set to Dirichlet on the left hand side
/// (x=0) and Neumann on the right (x=Nx-1/2). 
///
///  \param rfftw_many working_array -- solutions on the grid
///  \param rfftw_many boundary_working_array -- solutions at the boundary
///  \param field      Efield -- used for frequnetly used values
///
/// This solver does full dirichlet on both boundaries, and the potential 
/// difference is set by Ex0_external in the input file.


#include "efield_inject_tridiag.h" /* for field, types and fftw stuff */
#include "eppic-mpi.h" /* for subdomain stuff */
#include "eppic-math.h" /* for ceil function */
#include "tridiag.h"

void  efield_inject_tridiag_solve_dirichlet(rfftw_many &working_array,
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

  // first entry, includes:
  // set up b(y,z) terms 
  static bool first_entry=true;
  static ArrayNd<FTYPE,2> byz;

  FTYPE dx = Efield.dx;
  if(first_entry) {
    first_entry=false;

    index_type constants_size[] = {nNy,nNz};
    byz = ArrayNd<FTYPE,2>(constants_size);

    FTYPE dy = Efield.dy;
    FTYPE dz = Efield.dz;

    for(int iy=0;iy<nNy;iy++) 
      for(int iz=0;iz<nNz;iz++) {
      FTYPE dm=1.
	+2.*Sqr((dx/dy)*sin(M_PI* (FTYPE)iy /ny))
	+2.*Sqr((dx/dz)*sin(M_PI* (FTYPE)iz /nz));
      byz(iy,iz)=-2*dm;
      }
  }

  // setup values common to all boundary conditions
  if (n_solve == 0) 
    return;
  ArrayNd<FTYPE,2> a,b,c,lhs,rhs;
  a=ArrayNd<FTYPE,2>(n_solve,nx)= 1;
  b=ArrayNd<FTYPE,2>(n_solve,nx);
  c=ArrayNd<FTYPE,2>(n_solve,nx) = 1;
  rhs=ArrayNd<FTYPE,2>(n_solve,nx);
  FTYPE factor = -Efield.dx*Efield.dx/Efield.eps;
  int isolve=0;
  for (int iz=zstart;iz<zend;iz++) 
    for (int iy=ystart;iy<yend;iy++) {
      for (int ix=0;ix<nx;ix++) {
	rhs(isolve,ix) = 
	  working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)).real()*factor;
	b(isolve,ix) = byz(iy,iz);
      }
      isolve++;
      for (int ix=0;ix<nx;ix++) {
	rhs(isolve,ix) = 
	  working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)).imag()*factor;
	b(isolve,ix) = byz(iy,iz); // same for both real,imag
      }
      isolve++;
    }
  
  lhs=ArrayNd<FTYPE,2>(n_solve,nx);

  
  // set up left side boundaries
  if (subdomain.id_number==0) {
    isolve=0;
      for (int iz=zstart;iz<zend;iz++) 
	for (int iy=ystart;iy<yend;iy++) {
	  // iz=0,iy=0 special, but for speed done wrong here anyway
	  // corrected below
	  // first term
	  a(isolve,0) = 0;
	  b(isolve,0) = 1;
	  c(isolve,0) = 0;
	  rhs(isolve,0) = 0;
	  isolve++;

	  // second term
	  a(isolve,0) = 0;
	  b(isolve,0) = 1;
	  c(isolve,0) = 0;
	  rhs(isolve,0) = 0;
	  isolve++;
	}
      // kyz=0 term is special
      isolve=0;
      if (handlesKy0Eqn) {
	  // first term
	a(isolve,0) = 0.0;
	b(isolve,0) = 1;
	c(isolve,0) = 0;
	rhs(isolve,0) = Efield.boundary_values[0];  
	isolve++;
	  // this next eqn is meaningless, but if not treated properly,
	  // then the tridiagonal solve fails
	a(isolve,0) = 0.0;
	b(isolve,0) = 1;
	c(isolve,0) = 0;
	rhs(isolve,0) = Efield.boundary_values[0];  
	isolve++;
      }

  }

  // set up right side boundaries
  if (subdomain.id_number==nsubdomains-1) {
      isolve=0;
      for (int iz=zstart;iz<zend;iz++) 
	for (int iy=ystart;iy<yend;iy++) {
	  // iz=0,iy=0 special, but for speed done wrong here anyway
	  // corrected below
	  // first term
	  c(isolve,nx-1) = 0;
	  isolve++;
	  // second term
	  c(isolve,nx-1) = 0;
	  isolve++;
	}
      // kyz=0 term is special
      isolve=0;
      if (handlesKy0Eqn) {

	c(isolve,nx-1) = 0;
	isolve++;

	c(isolve,nx-1) = 0;
	  isolve++;

      }
  }

  // batch solve
  batch_tridiag_solver_parallel(a,b,c,rhs,lhs,subdomain.neighbor_comm);

  // - left side solution
  // do it first because need working array (ie rho) before it's overwritten
  isolve=0;
  if (subdomain.id_number==0) {
    for (int iz=zstart;iz<zend;iz++) 
      for (int iy=ystart;iy<yend;iy++) {
	// assumes phi(0) = 0
	boundary_working_array.cdata(RFFTW_C_INDICIES(0,iy,iz))= 
	  // -rho(0)*dx^2/eps:
	  working_array.cdata(RFFTW_C_INDICIES(0,iy,iz))*factor
	  // -phi(1):
	  -std::complex<FTYPE>(lhs(isolve,1),lhs(isolve+1,1));
	isolve+=2;
      }
    if (handlesKy0Eqn) {
	boundary_working_array.cdata(RFFTW_C_INDICIES(0,0,0))= 
	  // -rho(0)*dx^2/eps:
	  working_array.cdata(RFFTW_C_INDICIES(0,0,0))*factor	  
	  // -phi(1)
	  -std::complex<FTYPE>(lhs(0,1),lhs(1,1))
	  // +2*phi(0):
	  +2*Efield.boundary_values[0];
	  ;
    }
  }


  // record solutions
  isolve=0;
  for (int iz=zstart;iz<zend;iz++) 
    for (int iy=ystart;iy<yend;iy++) {
      for (int ix=0;ix<nx;ix++) 
	working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)) = 
	  std::complex<FTYPE>(lhs(isolve,ix),lhs(isolve+1,ix));
      isolve+=2;
    }


  // - right side soluion
  isolve=0;
  if (subdomain.id_number==nsubdomains-1) {
    for (int iz=zstart;iz<zend;iz++) 
      for (int iy=ystart;iy<yend;iy++) {
	// first term
	boundary_working_array.cdata(RFFTW_C_INDICIES(1,iy,iz)) = 0;
	// following assumes lhs(isolve,nx) == 0, ie previous line
	// also assumes rho(isolve,nx) == 0
	boundary_working_array.cdata(RFFTW_C_INDICIES(2,iy,iz)) = 
	  -std::complex<FTYPE>(lhs(isolve,nx-1),lhs(isolve+1,nx-1));
	isolve+=2;
      }
    if (handlesKy0Eqn) {
      boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,0)) = 
	std::complex<FTYPE>(0,0);
      // following assumes lhs(isolve,nx) == 0, ie previous line
      // also assumes rho(isolve,nx) == 0
      boundary_working_array.cdata(RFFTW_C_INDICIES(2,0,0)) = 
	-std::complex<FTYPE>(lhs(0,nx-1),lhs(1,nx-1));
    }
    
  }
  
	
}
