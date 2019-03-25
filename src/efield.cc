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

//#define write_dbg_efield // Turn this on to write out copious diagnostics

void efield(field &Efield, FArrayND &rho)
{

  /* This subroutine solves for the electric field, Ex, */
  /* phi, given the charge density, ne. */
  
  /* Local variables */

  FArrayND_ranged &phi = Efield.phi_rho;  // This is done for convenience only

  // phi  and ksqrinv must have some padding in the last dimension to 
  // allow fftw to perform an inplace transform

  /*
  FTYPE rho_dc;

  if (boundary_type==INJECT){ //between codes if INJECT FD, if PERIODIC FFT

    if (NDIM==1) {
      //void Tridiag(FArrayND &phi,FTYPE &rho_dc,FArrayND &rho,int b);
      //Tridiag(phi,rho_dc,rho,1);
      terminate(-1,"1D solver not yet implemented\n");
    } else if (NDIM==2) {
      void efield_inject(FArrayND_ranged &phi,FTYPE &rho_dc,FArrayND &rho,int b);
      void efield_inject_parallel(FArrayND_ranged &phi,FTYPE &rho_dc,FArrayND &rho,int b);
      if (nsubdomains>1) efield_inject_parallel(phi,rho_dc,rho,1);
      else efield_inject(phi,rho_dc,rho,1);
    }

  }else{ //between codes
  */
    const int nsize[]={INDICIES(nx*nsubdomains,ny,nz)};
#ifdef USE_DOMAINS
    static FArrayND_fftw phi_trans(subdomain.neighbor_comm, nsize);
    if (phi_trans.local_nx != nx) 
      terminate(-1,"Error: fftw_mpi is not assigning  Efield.phi to have nx=local_nx  as required\n");
    if (phi_trans.workspace == 0) phi_trans.make_workspace(nsize,0.);
    static FArrayND_fftw ksqrinv(subdomain.neighbor_comm, nsize);
    MPI_Request request1, request2;
    MPI_Status  status1,  status2;
#else
    static FArrayND phi_trans(nsize);
    static FArrayND_fftw ksqrinv(nsize);
#endif

#ifdef write_dbg_efield
    static ofstream diag;
#endif

    static bool first_entry=true;
    if (first_entry) {
      first_entry=false;

#ifdef write_dbg_efield
      //output this for diagnostic purposes (use 1 file per rank)
      std::stringstream fname;
      fname << "ksqrinv_mpi" <<mpi_rank<<".txt";
      string cname = fname.str();
      diag.open(cname.c_str());
#endif

      // The operator we need is just the fft of the finite difference opperator
      if (ksqrinv.x_start + nx == nsize[0] && nx > 0) 
	// Only the last processor needs this last index:
	ksqrinv(INDICIES(nx-1,0,0))+= 1/Sqr(dx);
      if (ksqrinv.x_start == 0 && nx > 0) {  // Only the first processor
	ksqrinv(INDICIES(0,0,0))   +=-2/Sqr(dx);
	ksqrinv(INDICIES(1%nx,0,0))+= 1/Sqr(dx);
	ksqrinv(INDICIES(0,0,0))   +=-2/Sqr(dy);
	ksqrinv(INDICIES(0,1%ny,0))+= 1/Sqr(dy);
	ksqrinv(INDICIES(0,ny-1,0))+= 1/Sqr(dy);
#if NDIM == 3
	ksqrinv(INDICIES(0,0,0))   +=-2/Sqr(dz);
	ksqrinv(INDICIES(0,0,1%nz))+= 1/Sqr(dz);
	ksqrinv(INDICIES(0,0,nz-1))+= 1/Sqr(dz);
#endif
      }
#ifdef write_dbg_efield
      for (int ix = 0; ix < ksqrinv.data.size(0); ++ix)
	for (int iy = 0; iy < ksqrinv.data.size(1); ++iy)
	  for (int iz = 0; iz < ksqrinv.data.size(2); iz++) {
	    diag <<"ksqr(pr="<< mpi_rank <<", "<< ix <<", "<< iy <<", "<< iz 
	      <<") = "<< ksqrinv(INDICIES(ix,iy,iz)) <<"\n";
	  }
#endif

      ksqrinv.transform();

      // ksqrinv is now a k^2 like operator. Invert, element by element
      // and make the complex = real part for rapid multiplication later.
#if NDIM==3
      for (int iy = 0; iy < ksqrinv.cdata.size(0); ++iy)
	for (int ix = 0; ix < ksqrinv.cdata.size(1); ++ix)
	  for (int iz = 0; iz < ksqrinv.cdata.size(2); iz++) {
#ifdef write_dbg_efield
	    diag <<"ksqrinv0(pr="<< mpi_rank <<", "<< ix <<", "<< iy <<", "<< iz 
	      <<") = "<< ksqrinv.cdata(INDICIES(iy,ix,iz)) <<"\n";
#endif
	    ksqrinv.cdata(iy,ix,iz)=1./ksqrinv.cdata(iy,ix,iz);
	  }
#elif NDIM==2
      for (int iy = 0; iy < ksqrinv.cdata.size(0); iy++) {
	for (int ix = 0; ix < ksqrinv.cdata.size(1); ++ix)
	  ksqrinv.cdata(iy,ix)=1./ksqrinv.cdata(iy,ix);
      }
#else
      for (int ix = 0; ix < nx/2+1; ++ix) {
	ksqrinv.cdata(ix)=1./ksqrinv.cdata(ix);
      }
#endif
    
      // Set DC component to zero:
      if (subdomain.id_number==0) 
	ksqrinv.cdata(INDICIES(0,0,0))=0.0;


      int nx_global=nx*nsubdomains;

      // Incorporate the epsilon coefficient and the n_elements factor into
      // ksqrinv:

      ksqrinv.data_transposed /= -eps* (nx*nsubdomains)*ny*nz;

      //To avoid taking an fft of ne twice, I'll do the high-frequency filtering
      // in this routine:

      // Calculate the damping coefficients

      if (fsteep != 0 && fwidth != 0) {
	FTYPE ksqr, damp, kx, ky, kz; 
	FTYPE kxbase=fwidth/(nx_global/2.), 
	  kybase=fwidth/(ny/2.), kzbase=fwidth/(nz/2.);

	for (int ikx=0; ikx<nx_global; ikx++) {
	  int ikx_global=ikx;
	  if (ikx_global <= nx_global/2) kx=kxbase*ikx_global;
	  else kx=kxbase*(-(nx_global-ikx_global));

	  for (int iky=0; iky<ksqrinv.ny_transpose; iky++) {
	    int iky_global=iky+ksqrinv.y_start_transpose;
	    if (iky_global <= ny/2) ky=kybase*iky_global;
	    else ky=kybase*(-(ny-iky_global));
	    for (int ikz=0; ikz<nz; ikz++) {
	      if (ikz <= nz/2) kz=kzbase*ikz;
	      else kz=kzbase*(-(nz-ikz));
	      ksqr=kx*kx+ky*ky+kz*kz;
	      damp=exp(-1*pow(ksqr,fsteep));
	      if (ndim == 1) {
		if (ikx_global <= nx_global/2) {
		  ksqrinv.cdata(INDICIES(ikx,0,0)) *= damp;
		}
	      }
	      else if (ndim == 2) {
		if (iky_global <= ny/2) {
		  ksqrinv.cdata(INDICIES(iky,ikx,0)) *= damp;
		}
	      } else {
		if (ikz <= nz/2) {
		  ksqrinv.cdata(INDICIES(iky,ikx,ikz)) *= damp;
		}
	      }
	    }
	  }
	}
      }

      // For some applications it is useful to eliminate 
      // parallel electric fields:
      if (no_parallel_efields && ksqrinv.x_start==0) {
	int ikx=0;
	for (int iky=0; iky<ksqrinv.ny_transpose; iky++) 
	  for (int ikz=0; ikz<nz; ikz++) 
	    ksqrinv.cdata(INDICIES(iky,ikx,ikz))=0.0;
      }
    
      // Eliminate modes with an angle closer than kill_modes_ang deg to z
      if (ndim == 3 && kill_modes_ang > 0) {
	FTYPE kill_modes_ang_rad=kill_modes_ang*M_PI/180.;
	FTYPE ksqr, kx, ky, kz; 
	FTYPE kxbase=2*M_PI/(nx_global*dx);
	FTYPE kybase=2*M_PI/(ny*dy);
	FTYPE kzbase=2*M_PI/(nz*dz);
	for (int ikx=0; ikx<nx_global; ikx++) {
	  int ikx_global=ikx;
	  if (ikx_global <= nx_global/2) kx=kxbase*ikx_global;
	  else kx=kxbase*(-(nx_global-ikx_global));
	  for (int iky=0; iky<ksqrinv.ny_transpose; iky++) {
	    int iky_global=iky+ksqrinv.y_start_transpose;
	    if (iky_global <= ny/2) ky=kybase*iky_global;
	    else ky=kybase*(-(ny-iky_global));
	    for (int ikz=0; ikz<nz; ikz++) {
	      if (ikz <= nz/2) kz=kzbase*ikz;
	      else kz=kzbase*(-(nz-ikz));
	      ksqr=kx*kx+ky*ky+kz*kz;
	      if (ksqr > 0.) {
		if (kill_modes_ang_rad > acos(kz/sqrt(ksqr)) ) {
		  if (ikz <= nz/2) {
		    ksqrinv.cdata(INDICIES(iky,ikx,ikz)) = 0.;
		  }
		}
	      }
	    }
	  }
	}
      }


      // Since ksqrinv is a real number, and in order to enable fast multiplication,
      // we will make the real and imag parts of ksqrinv be equal.
      int nksqrinv[]={INDICIES(ksqrinv.cdata.size(0),ksqrinv.cdata.size(1),
			       ksqrinv.cdata.size(2)),1,1};
      for (int ix = 0; ix < nksqrinv[0]; ++ix)
	for (int iy = 0; iy < nksqrinv[1]; ++iy)
	  for (int iz = 0; iz < nksqrinv[2]; iz++) {
	    ksqrinv.cdata(INDICIES(ix,iy,iz))=
	      std::complex<FTYPE>(ksqrinv.cdata(INDICIES(ix,iy,iz)).real(),
				  ksqrinv.cdata(INDICIES(ix,iy,iz)).real());
	    //output this for diagnostic purposes
#ifdef write_dbg_efield
	    diag <<"ksqrinv(pr="<< mpi_rank <<", "<< ix <<", "<< iy <<", "<< iz 
	      <<") = "<< ksqrinv.cdata(INDICIES(ix,iy,iz)) <<"\n";
#endif
	  }

      //      diag.close();

    } //END if (first_entry) 


    // Calculate the rhs of poisson's equation in MKS units: 

    /* TEST ROUTINE: initialize data to a simple wave function: */
    /*  for (int ix = 0; ix < nx; ++ix)
	for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) {
	phi(ix, iy, iz)
	= sin((ix)/(nx/(2.*M_PI)));
	cout <<"0 "<< mpi_rank <<" "<< ix <<" "<< iy <<" "<< iz <<" "<< 
	phi_trans(ix,iy,iz) <<"\n";
	} */

    // Copy phi into phi_trans, setting all non-overlapping values to zero
    phi_trans.data=0.0;

    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) {
	  phi_trans(INDICIES(ix,iy,iz))=rho(INDICIES(ix,iy,iz));
#ifdef write_dbg_efield
	    diag <<"rho(pr="<< mpi_rank <<", "<< ix <<", "<< iy <<", "<< iz 
		 <<") = "<< phi_trans(INDICIES(ix,iy,iz)) << "\n";
#endif
	}


    // Transform into k space
    phi_trans.transform();
  
    // Store the DC component: Only correct on the first subdomain
    //  FTYPE phi_trans_dc=0;
    //  if (subdomain.id_number == 0) phi_trans_dc=phi_trans(INDICIES(0,0,0))/(nx*ny*nz);
  
      int nksqrinv[]={INDICIES(ksqrinv.cdata.size(0),ksqrinv.cdata.size(1),
			       ksqrinv.cdata.size(2)),1,1};
#ifdef write_dbg_efield
      for (int ix = 0; ix < nksqrinv[0]; ++ix)
	for (int iy = 0; iy < nksqrinv[1]; ++iy)
	  for (int iz = 0; iz < nksqrinv[2]; iz++) {
	    diag <<"phi(pr="<< mpi_rank <<", "<< ix <<", "<< iy <<", "<< iz 
		 <<") = "<< phi_trans.cdata(INDICIES(ix,iy,iz)) << " * " << ksqrinv.cdata(INDICIES(ix,iy,iz)) <<"\n";
	  }
#endif

    // Calculate phi_trans.transform/ (-k^2 * eps*nx*ny*nz)
    phi_trans.data_transposed *= ksqrinv.data_transposed;
  
#ifdef write_dbg_efield
      for (int ix = 0; ix < nksqrinv[0]; ++ix)
	for (int iy = 0; iy < nksqrinv[1]; ++iy)
	  for (int iz = 0; iz < nksqrinv[2]; iz++) {
	    diag <<"phiFFT(pr="<< mpi_rank <<", "<< ix <<", "<< iy <<", "<< iz 
		 <<") = "<< phi_trans.cdata(INDICIES(ix,iy,iz)) <<"\n";
	  }
#endif

    // Inverse fft transform:
    phi_trans.invtransform();

#ifdef write_dbg_efield
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) {
	    diag <<"phi_f(pr="<< mpi_rank <<", "<< ix <<", "<< iy <<", "<< iz 
		 <<") = "<< phi_trans(INDICIES(ix,iy,iz)) << "\n";
	}
#endif

    // Test the solution:
    /* Needs implementing:
       static int test_done=FALSE;
       if (!test_done && !(fsteep != 0 && fwidth != 0)) {
       test_done=TRUE;
       if (phi_trans.address() != phi.address()) {
       FTYPE rho2, rho_start;
       int ixm1=nx-1;
       for (int ix=0; ix<nx; ix++) {
       int iym1=ny-1;
       for (int iy=0; iy<ny; iy++) {
       int izm1=nz-1;
       for (int iz=0; iz<nz; iz++) {
       rho2=0;
       rho2 += phi_trans(INDICIES(ix,iy,iz))*(-2/Sqr(dx));
       rho2 += phi_trans(INDICIES((ix+1)%nx,iy,iz))/Sqr(dx);
       rho2 += phi_trans(INDICIES(ixm1,iy,iz))/Sqr(dx);
	    
       #if NDIM > 1
       rho2 += phi_trans(INDICIES(ix,iy,iz))*(-2/Sqr(dy));
       rho2 += phi_trans(INDICIES(ix,(iy+1)%ny,iz))/Sqr(dy);
       rho2 += phi_trans(INDICIES(ix,iym1,iz))/Sqr(dy);
       #endif
	  
       #if NDIM == 3
       rho2 += phi_trans(INDICIES(ix,iy,iz))*(-2/Sqr(dz));
       rho2 += phi_trans(INDICIES(ix,iy,(iz+1)%nz))/Sqr(dz);
       rho2 += phi_trans(INDICIES(ix,iy,izm1))/Sqr(dz);
       #endif	  
       rho2 *= -1.0*eps;
       rho_start=phi(INDICIES(ix,iy,iz));
       if (fabs(rho2-(rho_start-rho_dc))/sqrt(Sqr(rho2)+Sqr(rho_start-rho_dc))>1e-5){
       char message[1024];
       sprintf(message,
       "efield.c incorrect: rho(%d,%d,%d)=%g  -del^2(phi_trans)*eps=%g\n",
       ix,iy,iz,rho_start-rho_dc,rho2);
       printf(message);
       }
       izm1=iz;
       }
       iym1=iy;
       }
       ixm1=ix;
       }
       }
       }

    */

    // Copy solution into phi

#ifdef USE_DOMAINS
    // This section passes guard cell(s) from an array on one processor 
    // to another on another processor.  It only passes the last dimension(s).
    // nx_guard[0] indicates how many guard cells to be passed on the LHS
    // nx_guard[1] indicates how many guard cells to be passed on the RHS
    // This can be done before copying the rest of the array over.

    //  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const unsigned int nx_guard[]);
    //  pass_guards(Efield.phi.data, nxguard_size);

    int size_send=phi.length()/phi.size(0);

    //Copy the RHS components which will be sent into the phi array (with no padding)
    for (int ix = nx - phix_guard_size[0]; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) 
	  phi(INDICIES(ix,iy,iz))=phi_trans(INDICIES(ix,iy,iz));
    }
    // Pass RHS to the East:

    const int index0[]={INDICIES(nx-phix_guard_size[0],0,0)};
    MPI_Isend(phi.address(index0), size_send*phix_guard_size[0], 
	      MPI_FTYPE, proc_east, 2, MPI_COMM_WORLD, &request2);

    // Recieve this data from the West and place in the LHS:
    const int index1[]={INDICIES(-((int)phix_guard_size[0]),0,0)};
    MPI_Recv(phi.address(index1), size_send*phix_guard_size[0], 
	     MPI_FTYPE, proc_west, 2, MPI_COMM_WORLD, &status2);
  
    //Copy the LHS components which will be sent into phi with no padding
    for (int ix = 0; ix < phix_guard_size[1]; ++ix)
      for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) 
	  phi(INDICIES(ix,iy,iz))=phi_trans(INDICIES(ix,iy,iz));
    // Pass LHS to the West:

    const int index2[]={INDICIES(0,0,0)};
    MPI_Isend(phi.address(index2), size_send*phix_guard_size[1], 
	      MPI_FTYPE, proc_west, 1, MPI_COMM_WORLD, &request1);

    // Recieve this data from the East and place in the RHS:
    const int index3[]={INDICIES(nx,0,0)};
    MPI_Recv(phi.address(index3), size_send*phix_guard_size[1], 
	     MPI_FTYPE, proc_east, 1, MPI_COMM_WORLD, &status1);

#endif

    //Copy the rest of the phi_trans array to phi
    for (int ix = phix_guard_size[1]; ix < nx-phix_guard_size[0]; ++ix)
      for (int iy = 0; iy < ny; ++iy)
	for (int iz = 0; iz < nz; ++iz) {
	  phi(INDICIES(ix,iy,iz))=phi_trans(INDICIES(ix,iy,iz));
	}

#ifdef USE_DOMAINS

    MPI_Wait(&request1, &status1);
    MPI_Wait(&request2, &status2);
#endif

    //  We don't need Ex because vpush is faster working with phi
    if (Efield.Ex.size() == Efield.phi_rho.size()) {
      void gradient(FArrayND &phi, 
		    INDICIES(FArrayND &Ex, FArrayND &Ey, FArrayND &Ez), 
		    INDICIES(FTYPE dx, FTYPE dy, FTYPE dz), FTYPE scaler);
      gradient(phi, INDICIES(Efield.Ex, Efield.Ey, Efield.Ez), 
	       INDICIES(dx, dy, dz), (FTYPE) -1.); 
      // Impose any External or driving E-field 
      //  (changes to phi must be periodic!)
      void externale(field &Efield);
      externale(Efield);

    }

#ifdef write_dbg_efield
    std::stringstream fname2;
    fname2 << "phi0_it"<< it << "_mpi" <<mpi_rank<<".txt";
    string cname2 = fname2.str();
    phi.output(cname2.c_str());
#endif

    // Clear workspace when this is an output timestep to make space for output work!
    if (it%nout == 0) phi_trans.del_workspace();

    // Let's time the amount of time in this routine:
  
    /*  
} //end "between codes" block
*/

}
