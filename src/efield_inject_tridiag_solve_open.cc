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
/// This is the solver for rho=0 outside the simulation;
/// It is currently not working.

#include "efield_inject_tridiag.h" /* for field, types and fftw stuff */
#include "eppic-mpi.h" /* for subdomain stuff */
#include "eppic-math.h" /* for ceil function */
#include "tridiag.h"

void  efield_inject_tridiag_solve_open(rfftw_many &working_array,
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


  // if either side open, setup ky0 stuff
  ArrayNd<FTYPE,1> ky0_boundary=ArrayNd<FTYPE,1>(4)=0; // 0,1-left//2,3-right

  if (handlesKy0Eqn) {
    ArrayNd<FTYPE,1> ky0_boundary_buffer=ArrayNd<FTYPE,1>(4)=0;
    FTYPE xshift=0;
      
    // xshift for domains and 'j'=-1, see formula cited in birdsall book
    xshift=subdomain.id_number*nx; 
    FTYPE nxTotal = nx*nsubdomains;
    int ixstart=0;
    // this next line forces a skip of ixGlobal=0, which assumes that 
    // rho(ixGlobal=0) is 0, so adding tha term is not needed. 
    if (subdomain.id_number==0) ixstart=0;
    
    for (int ix=ixstart; ix<nx; ix++) {
      // corresponds to j=-1
      //               (ix_global-j)
      ky0_boundary(0)+=(ix+xshift+1)*
	working_array.cdata(RFFTW_C_INDICIES(ix,0,0)).real(); 
      // corresponds to j=0
      //               (ix_global-j)
      ky0_boundary(1)+=(ix+xshift)*
	working_array.cdata(RFFTW_C_INDICIES(ix,0,0)).real();
    }
    
    /* shift for 'j'=nx_total */
    ixstart=0;
    if (subdomain.id_number==0) ixstart=0;
    for (int ix=ixstart; ix<nx; ix++) {
      // corresponds to j=nx
      //               (j      -ixGlobal )
      ky0_boundary(2)+=(nxTotal-xshift-ix)*
	working_array.cdata(RFFTW_C_INDICIES(ix,0,0)).real();
      // corresponds to j=nx+1
      //               (j      -ixGlobal )
      ky0_boundary(3)+=(nxTotal+1-xshift-ix)*
	working_array.cdata(RFFTW_C_INDICIES(ix,0,0)).real();
    }
    if (nsubdomains>1) {
      ky0_boundary_buffer=ky0_boundary;
      MPI_Reduce(ky0_boundary_buffer.address(0),
		 ky0_boundary.address(0),2,MPI_FTYPE,MPI_SUM,0,
		 subdomain.neighbor_comm);
      MPI_Reduce(ky0_boundary_buffer.address(2),
		 ky0_boundary.address(2),2,MPI_FTYPE,MPI_SUM,nsubdomains-1,
		 subdomain.neighbor_comm);
    }
    


    ky0_boundary(0)-=ky0_boundary(2);
    cout 
      << " ky 0 = " << ky0_boundary(0)
      << " ky 1 = " << ky0_boundary(1)
      << " ky 2 = " << ky0_boundary(2)
      << endl;
    ky0_boundary*=(-Efield.dx*Efield.dx/Efield.eps/2.);
    
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

  cout 
    << "rho 0 = " << rhs(0,0)
    << "rho nx-1 = " << rhs(0,nx-1)
    << endl;

  
  // set up left side boundaries
  if (subdomain.id_number==0) {
    isolve=0;
      for (int iz=zstart;iz<zend;iz++) 
	for (int iy=ystart;iy<yend;iy++) {
	  // iz=0,iy=0 special, but for speed done wrong here anyway
	  // corrected below
	  // first term
// 	  a(isolve,0) = 0;
// 	  b(isolve,0) = 1;
// 	  c(isolve,0) = 0;
// 	  rhs(isolve,0) = 1; // will redefine after the solve
// 	  a(isolve,1) = 0;
// 	  b(isolve,1) = byz(iy,iz)+rinv(iy,iz);
// 	  isolve++;
	  a(isolve,0) = 0;
	  b(isolve,0) = byz(iy,iz)+rinv(iy,iz);
	  isolve++;

	  // second term
// 	  a(isolve,0) = 0;
// 	  b(isolve,0) = 1;
// 	  c(isolve,0) = 0;
// 	  rhs(isolve,0) = 1; // will redefine after the solve
// 	  a(isolve,1) = 0;
// 	  b(isolve,1) = byz(iy,iz)+rinv(iy,iz);
// 	  isolve++;
	  a(isolve,0) = 0;
	  b(isolve,0) = byz(iy,iz)+rinv(iy,iz);
	  isolve++;
	}
      // kyz=0 term is special
      isolve=0;
      if (handlesKy0Eqn) {
	// equations setup for x=0 ignored:
// 	// first term
// 	a(isolve,0) = 0.0;
// 	b(isolve,0) = 1.0;
// 	c(isolve,0) = 0.0;
// 	rhs(isolve,0) = ky0_boundary(1);
// 	b(isolve,1) = byz(0,0);
// 	a(isolve,1) = 0.0;
// 	rhs(isolve,1) -= ky0_boundary(1);
// 	isolve++;

// 	// this next eqn is meaningless, but if not treated properly,
// 	// then the tridiagonal solve fails
// 	a(isolve,0) = 0.0;
// 	b(isolve,0) = 1.0;
// 	c(isolve,0) = 0.0;
// 	rhs(isolve,0) = ky0_boundary(1);
// 	b(isolve,1) = byz(0,0);
// 	a(isolve,1) = 0.0;
// 	rhs(isolve,1) -= ky0_boundary(1);
// 	isolve++;

	// equations setup for x=0 included
	b(isolve,0) = byz(0,0);
	a(isolve,0) = 0.0;
	rhs(isolve,0) -= ky0_boundary(0);
	isolve++;

	// this next eqn is meaningless, but if not treated properly,
	// then the tridiagonal solve fails
	b(isolve,0) = byz(0,0);
	a(isolve,0) = 0.0;
	rhs(isolve,0) -= ky0_boundary(0);
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
	  b(isolve,nx-1)+=rinv(iy,iz);
	  isolve++;
	  // second term
	  c(isolve,nx-1) = 0;
	  b(isolve,nx-1)+=rinv(iy,iz);
	  isolve++;
	}
      // kyz=0 term is special
      isolve=0;
      if (handlesKy0Eqn) {
	// first term
	rhs(isolve,nx-1)-=ky0_boundary(2);
	c(isolve,nx-1) = 0;
	b(isolve,nx-1) = byz(0,0);
	isolve++;
	rhs(isolve,nx-1)-=ky0_boundary(2);
	c(isolve,nx-1) = 0;
	b(isolve,nx-1) = byz(0,0);

	isolve++;


      }
  }

  // batch solve
  batch_tridiag_solver_parallel(a,b,c,rhs,lhs,subdomain.neighbor_comm);

  // record solutions
  isolve=0;
  for (int iz=zstart;iz<zend;iz++) 
    for (int iy=ystart;iy<yend;iy++) {
      for (int ix=0;ix<nx;ix++) 
	working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)) = 
	  std::complex<FTYPE>(lhs(isolve,ix),lhs(isolve+1,ix));
      isolve+=2;
    }

  // - left side solution

  if (subdomain.id_number==0) {
    isolve=0;
    for (int iz=zstart;iz<zend;iz++) 
      for (int iy=ystart;iy<yend;iy++) {
	
// 	// first term
// 	working_array.cdata(RFFTW_C_INDICIES(0,iy,iz)) = 
// 	  std::complex<FTYPE>(lhs(isolve,1)*rinv(iy,iz),
// 			      lhs(isolve+1,1)*rinv(iy,iz));
	boundary_working_array.cdata(RFFTW_C_INDICIES(0,iy,iz))= 
	  std::complex<FTYPE>(lhs(isolve,1)*rinv(iy,iz)*rinv(iy,iz),
			      lhs(isolve,1)*rinv(iy,iz)*rinv(iy,iz));
	isolve++;
	isolve++;
      }
    isolve=0;
    if (handlesKy0Eqn) {
      // first term
//       working_array.cdata(RFFTW_C_INDICIES(0,0,0)) = 
// 	std::complex<FTYPE>(ky0_boundary(1),ky0_boundary(1));
      boundary_working_array.cdata(RFFTW_C_INDICIES(0,0,0)) = 
	std::complex<FTYPE>(ky0_boundary(0),ky0_boundary(0));

      isolve++;
    }
  }

  // right side soluion
  if (subdomain.id_number==nsubdomains-1) {
    isolve=0;
    for (int iz=zstart;iz<zend;iz++) 
      for (int iy=ystart;iy<yend;iy++) {
	// first term
	boundary_working_array.cdata(RFFTW_C_INDICIES(1,iy,iz)) = 
	  std::complex<FTYPE>(lhs(isolve,nx-1)*rinv(iy,iz),
			      lhs(isolve+1,nx-1)*rinv(iy,iz));
	boundary_working_array.cdata(RFFTW_C_INDICIES(2,iy,iz)) = 
	  std::complex<FTYPE>(lhs(isolve,nx-1)*rinv(iy,iz)*rinv(iy,iz),
			      lhs(isolve+1,nx-1)*rinv(iy,iz)*rinv(iy,iz));
	isolve++;
	isolve++;
      }
    if (handlesKy0Eqn) {
      // first term
      boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,0)) = 
	std::complex<FTYPE>(ky0_boundary(2),ky0_boundary(2));
      boundary_working_array.cdata(RFFTW_C_INDICIES(2,0,0)) = 
	std::complex<FTYPE>(ky0_boundary(3),ky0_boundary(3));
    }
    
  }
  
	
}
