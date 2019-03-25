

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
///  \todo This routine is very long, and could be broken up into several 
///        subroutine calls, eg setup and post processing, perhaps each with
///        their own subroutine calls. 
///
///  \todo Verify equations for Neumann case
///
///  \todo Add start and end variables for Y and Z, remove current awkward 
///        setup. In 2D, the most slowly changing index is Y, in 3D the most 
///        slowly changing index is Z. The slowest index is used to 
///        subdomain-decompose. Currently there are CPP if statements all 
///        over the place to decide between 2D and 3D, but this should be done
///        once at the top along with start and end variables for Y and Z. 
///  
///  \todo Update documentation somewhere that defines the equations being 
///        solved. 

#include "efield_inject_tridiag.h" /* for field, types and fftw stuff */
#include "eppic-mpi.h" /* for subdomain stuff */
#include "eppic-math.h" /* for ceil function */
#include "tridiag.h"

void  efield_inject_tridiag_solve_batch(rfftw_many &working_array,
				  rfftw_many &boundary_working_array,
				  field &Efield)
{

  // initialize dimensions
  int nx = Efield.nx;
  int ny = Efield.ny;
  int nz = Efield.nz;

#if NDIM==2
  int nNy=ny/2+1;
#else
  int nNy=ny;
#endif
  int nNz=nz/2+1;

  // divide workload in slowest varying direction
  // in 2D, this is y, in 3D this is z
  // division along slowest leads to most efficient MPI communications
  // we do real and imaginary component in pairs 
#if NDIM==2  
  int dnslow = nNy/subdomain.np;
  // old way: static_cast<int>(ceil((nNy)/subdomain.np));
  int slow_start = dnslow*subdomain.rank;
  int slow_end = slow_start+dnslow;
  if (slow_end +dnslow > nNy) slow_end = nNy; 
#else
  int slow_start = working_array.y_start_transpose;
  int slow_end = slow_start+working_array.ny_transpose;
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

#if NDIM==2
  int n_solve = (slow_end-slow_start)*nNz*2; 
  //'*2' applies to real and imag terms
#else
  int n_solve = (slow_end-slow_start)*nNy*2; 
  //'*2' applies to real and imag terms
#endif
  if (slow_start==0) {// ie only relevent if ky0 is done by this proc
    ArrayNd<FTYPE,1> ky0_boundary_buffer=ArrayNd<FTYPE,1>(4)=0;

    if ((Efield.boundary[0]==open)||(Efield.boundary[1]==open)) {
      
      
      FTYPE xshift=0;
      
      if (Efield.boundary[0]==open) {
	// xshift for domains and 'j'=-1, see formula cited in birdsall book
	xshift=subdomain.id_number*nx+1; 
	int ixstart=0;
	if (subdomain.id_number==0) ixstart=1;
	
	for (int ix=ixstart; ix<nx; ix++) {
	  // corresponds to x=-1
	  ky0_boundary(0)+=(1+ix+xshift)*
	    working_array.cdata(RFFTW_C_INDICIES(ix,0,0)).real(); 
	  // corresponds to x=0
	  ky0_boundary(1)+=(ix+xshift)*
	    working_array.cdata(RFFTW_C_INDICIES(ix,0,0)).real();
	}
      }
      
      if (Efield.boundary[1]==open) {
	/* shift for 'j'=nx_total */
	xshift=nx*nsubdomains-nx*subdomain.id_number; 
	int ixstart=0;
	if (subdomain.id_number==0) ixstart=1;
	for (int ix=ixstart; ix<nx; ix++) {
	  // corresponds to x=nx
	  ky0_boundary(2)+=(xshift-ix)*
	    working_array.cdata(RFFTW_C_INDICIES(ix,0,0)).real();
	// corresponds to x=nx+1
	  ky0_boundary(3)+=(1+xshift-ix)*
	    working_array.cdata(RFFTW_C_INDICIES(ix,0,0)).real();
	}
      }
    
    
      if (nsubdomains>1) {
	ky0_boundary_buffer=ky0_boundary;
	if (Efield.boundary[0]==open)
	  MPI_Reduce(ky0_boundary_buffer.address(0),
		     ky0_boundary.address(0),2,MPI_FTYPE,MPI_SUM,0,
		     subdomain.neighbor_comm);
	if (Efield.boundary[1]==open)
	  MPI_Reduce(ky0_boundary_buffer.address(2),
		     ky0_boundary.address(2),2,MPI_FTYPE,MPI_SUM,nsubdomains-1,
		     subdomain.neighbor_comm);
	
      }
    }
  }

  // setup values common to all boundary conditions
  if (n_solve == 0) 
    return;
  ArrayNd<FTYPE,2> a,b,c,lhs,rhs;
  a=ArrayNd<FTYPE,2>(n_solve,nx)= 1;
  b=ArrayNd<FTYPE,2>(n_solve,nx);
  c=ArrayNd<FTYPE,2>(n_solve,nx) = 1;
  ArrayNd<FTYPE,1> alpha,beta;
  if (Efield.boundary[0]==periodic) {
    alpha = ArrayNd<FTYPE,1>(n_solve) = 1;
    beta = alpha;
    if (slow_start==0) {
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
  rhs=ArrayNd<FTYPE,2>(n_solve,nx);
  FTYPE factor = -Efield.dx*Efield.dx/Efield.eps;
  int isolve=0;
#if NDIM==2  
  for (int iz=0;iz<1;iz++) 
    for (int iy=slow_start;iy<slow_end;iy++) {
#else
  for (int iz=slow_start;iz<slow_end;iz++) 
    for (int iy=0;iy<ny;iy++) {
#endif
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
    switch(Efield.boundary[0]) {
    case open:
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	  // iz=0,iy=0 special, but for speed done wrong here anyway
	  // corrected below
	  // first term
	  a(isolve,0) = 0;
	  b(isolve,0) = 1;
	  c(isolve,0) = 0;
	  rhs(isolve,0) = 1; // will redefine after the solve
	  a(isolve,1) = 0;
	  b(isolve,1) = byz(iy,iz)+rinv(iy,iz);
	  isolve++;

	  // second term
	  a(isolve,0) = 0;
	  b(isolve,0) = 1;
	  c(isolve,0) = 0;
	  rhs(isolve,0) = 1; // will redefine after the solve
	  a(isolve,1) = 0;
	  b(isolve,1) = byz(iy,iz)+rinv(iy,iz);
	  isolve++;
	}
      // kyz=0 term is special
      isolve=0;
      if (slow_start==0) {
	// first term
	a(isolve,0) = 0.0;
	b(isolve,0) = 1.0;
	c(isolve,0) = 0.0;
	rhs(isolve,0) = ky0_boundary(1);
	b(isolve,1) = byz(0,0);
	a(isolve,1) = 0.0;
	rhs(isolve,1) -= ky0_boundary(1);
	isolve++;
	// this next eqn is meaningless, but if not treated properly,
	// then the tridiagonal solve fails
	a(isolve,0) = 0.0;
	b(isolve,0) = 1.0;
	c(isolve,0) = 0.0;
	rhs(isolve,0) = ky0_boundary(1);
	b(isolve,1) = byz(0,0);
	a(isolve,1) = 0.0;
	rhs(isolve,1) -= ky0_boundary(1);
	isolve++;
      }
      break;

    case dirichlet:
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	  // iz=0,iy=0 special, but for speed done wrong here anyway
	  // corrected below
	  // first term
	  a(isolve,0) = 0;
	  b(isolve,0) = 1;
	  c(isolve,0) = 0;
	  rhs(isolve,0) = 0;
	  a(isolve,1) = 0;
	  isolve++;
	  // second term
	  a(isolve,0) = 0;
	  b(isolve,0) = 1;
	  c(isolve,0) = 0;
	  rhs(isolve,0) = 0;
	  a(isolve,1) = 0;
	  isolve++;
	}
      // kyz=0 term is special
      isolve=0;
      if (slow_start==0) {
	// first term
	a(isolve,0) = 0.0;
	b(isolve,0) = 1.0;
	c(isolve,0) = 0.0;
	rhs(isolve,0) = Efield.boundary_values[1];
	a(isolve,1) = 0.0;
	rhs(isolve,1) -= Efield.boundary_values[1];
	isolve++;
	a(isolve,0) = 0.0;
	b(isolve,0) = 1.0;
	c(isolve,0) = 0.0;
	rhs(isolve,0) = Efield.boundary_values[1];
	a(isolve,1) = 0.0;
	rhs(isolve,1) -= Efield.boundary_values[1];
	isolve++;
	// second term;
	// nothing to do, correct above
	isolve++;
      }
      break;
    case neumann: // very likely not correct!
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	  // iz=0,iy=0 special, but for speed done wrong here anyway
	  // corrected below
	  // first term
	  a(isolve,0) = 0;
	  b(isolve,0)++;
	  isolve++;
	  // second term
	  a(isolve,0) = 0;
	  b(isolve,0)++;
	  isolve++;
	}
      // kyz=0 term is NOT special for neumann (apparently!)
      break;
    case periodic:
      // nothing to do for non kyz==0 terms
      // kyz=0 term is special
      isolve=0;
      if (slow_start==0) {
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
	isolve++;
	// second term;
	// nothing special for this term, done correctly above!
      }
    break;
    }
  }

  // set up right side boundaries
  if (subdomain.id_number==nsubdomains-1) {
    switch(Efield.boundary[0]) {
    case open:
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
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
      if (slow_start==0) {
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
      
      break;
    case dirichlet:
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
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
      if (slow_start==0) {
	// first term
	c(isolve,nx-1) = 0;
	rhs(isolve,nx-1) -= Efield.boundary_values[2];
	isolve++;
	c(isolve,nx-1) = 0;
	rhs(isolve,nx-1) -= Efield.boundary_values[2];
	isolve++;
	// second term;
	// same as typical case, done already
      }

      break;
    case neumann: // very likely not correct!
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	  // iz=0,iy=0 special, but for speed done wrong here anyway
	  // corrected below
	  // first term
	  c(isolve,nx-1) = 0;
	  b(isolve,nx-1)++;
	  isolve++;
	  // second term
	  c(isolve,nx-1) = 0;
	  b(isolve,nx-1)++;
	  isolve++;
	}
      // kyz=0 term is NOT special for neumann (apparently!)
      break;
    case periodic:
      // nothing to do for non kyz==0 terms
      // kyz=0 term is not special
      isolve=0;
      if (slow_start==0) {
	// first term
	alpha(isolve) = 0;
	beta(isolve) = 0;
	isolve++;
	alpha(isolve) = 0;
	beta(isolve) = 0;
	// second term;
	// nothing special for this term, done correctly above!
      }
      break;
    }
  }

  // batch solve
  if (Efield.boundary[0] == periodic) // then [1] term must be periodic too
    batch_tridiag_cyclic_solver_parallel(a,b,c,rhs,lhs,alpha,beta,
					 subdomain.neighbor_comm);
  else
    batch_tridiag_solver_parallel(a,b,c,rhs,lhs,subdomain.neighbor_comm);

  // record solutions
  // - generical solution record
  isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	for (int ix=0;ix<nx;ix++) 
	  working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)) = 
	    std::complex<FTYPE>(lhs(isolve,ix),lhs(isolve+1,ix));
	//isolve++;
	//for (int ix=0;ix<nx;ix++) 
	//working_array.cdata(RFFTW_C_INDICIES(ix,iy,iz)).imag() = 
	//lhs(isolve,ix);
	isolve++;
	isolve++;
    }
  
  
  // - left side solution
  if (subdomain.id_number==0) {
    switch(Efield.boundary[0]) {
    case open:
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	  // first term
	  working_array.cdata(RFFTW_C_INDICIES(0,iy,iz)) = 
	    std::complex<FTYPE>(lhs(isolve,1)*rinv(iy,iz),
				lhs(isolve+1,1)*rinv(iy,iz));
	  boundary_working_array.cdata(RFFTW_C_INDICIES(0,iy,iz))= 
	    std::complex<FTYPE>(lhs(isolve,1)*rinv(iy,iz)*rinv(iy,iz),
				lhs(isolve,1)*rinv(iy,iz)*rinv(iy,iz));
	  isolve++;
	  isolve++;
	}
      isolve=0;
      if (slow_start==0) {
	// first term
	working_array.cdata(RFFTW_C_INDICIES(0,0,0)) = 
	  std::complex<FTYPE>(ky0_boundary(1),ky0_boundary(0));
	isolve++;
      }

      break;
    case dirichlet:
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	  // first term
	  boundary_working_array.cdata(RFFTW_C_INDICIES(0,iy,iz)) = 
	    std::complex<FTYPE>(0,0);
	isolve++;
	isolve++;
	}
      if (slow_start==0) {
	// first term
	boundary_working_array.cdata(RFFTW_C_INDICIES(0,0,0)) = 
	  std::complex<FTYPE>(Efield.boundary_values[0],0);
	isolve++;
	isolve++;
      }
      
      break;
    case neumann: // very likely not correct!
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	  // first term
	  boundary_working_array.cdata(RFFTW_C_INDICIES(0,iy,iz)) = 
	    std::complex<FTYPE>(lhs(isolve,0),lhs(isolve+1,0));
	  isolve++;
	  isolve++;
	}
      break;
    case periodic:
      break;
    }
  }

  // right side soluion
  if (subdomain.id_number==nsubdomains-1) {
    switch(Efield.boundary[0]) {
    case open:
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
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
      if (slow_start==0) {
	// first term
	boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,0)) = 
	  std::complex<FTYPE>(ky0_boundary(2),0);
	boundary_working_array.cdata(RFFTW_C_INDICIES(2,0,0)) = 
	  std::complex<FTYPE>(ky0_boundary(3),0);
      }

      break;
    case dirichlet:
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	  //first term
	  boundary_working_array.cdata(RFFTW_C_INDICIES(1,iy,iz)) =
	    std::complex<FTYPE>(0,0);
	}
      if (slow_start==0) {
	// first term
	boundary_working_array.cdata(RFFTW_C_INDICIES(1,0,0)) = 
	  std::complex<FTYPE>(Efield.boundary_values[2],0);
	boundary_working_array.cdata(RFFTW_C_INDICIES(2,0,0)) = 
	  std::complex<FTYPE>(Efield.boundary_values[3],0);
	isolve++;
	isolve++;
      }

      break;
    case neumann: // very likely not correct!
      isolve=0;
#if NDIM==2
      for (int iz=0;iz<1;iz++) 
	for (int iy=slow_start;iy<slow_end;iy++) {
#else 
      for (int iz=slow_start;iz<slow_end;iz++) 
	for (int iy=0;iy<ny;iy++) {
#endif
	  // first term
	  boundary_working_array.cdata(RFFTW_C_INDICIES(1,iy,iz)) = 
	    std::complex<FTYPE>(lhs(isolve,nx-1),lhs(isolve+1,nx-1));
	  isolve++;
	  isolve++;
	}
      break;
    case periodic:
      break;
    }
  }
}
