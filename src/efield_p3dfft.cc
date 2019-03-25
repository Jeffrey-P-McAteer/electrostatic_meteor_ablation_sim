  /* This subroutine solves for the electric potential phi, given the charge density, ne. 
  
     It uses the p3dfft library which requires a 2D decomposition of the 3D data.  
     It only works for 3D data. 
  */

#include <stdio.h>
#include <math.h>
#include <sys/times.h> //Clock related functions
#include "eppic.h"
#include "eppic-mpi.h"

#include <complex>

//#define write_dbg_p3d

#if NDIM != 3
// I tried to get this to work in 2D but it would not do the decomposition in a useful way
void efield_p3dfft(field &Efield, FArrayND &rho){
  terminate(-1,"Error: NDIM = 1, efield_p3dfft not applicable in 1D or 2D");
};
#else

#ifndef USE_P3DFFT
void efield_p3dfft(field &Efield, FArrayND &rho){
  terminate(-1,"Error: efield_p3dfft not compiled in");
};
#else

extern "C" {
#include "p3dfft.h"
}

int in_range(const int *ipoint, const int *imin, const int *imax) {
  int in_range=TRUE;
  for (int id=0; id<NDIM;id++) {
    if ( !(ipoint[id]>=imin[id] && ipoint[id]<=imax[id]) ) {
      in_range=FALSE;
      break;
    }
  }
  return in_range;
};

template <typename TYPE>
void add_if_in_range(TYPE &array, INDICIES(int nx, int ny, int nz) , const int *istart, const int *iend, FTYPE amount){
  int global_index[NDIM]={INDICIES(nx,ny,nz)};
  if (in_range(global_index,istart,iend)) {
    int local_index[]={INDICIES(nx-istart[0],ny-istart[1],nz-istart[2])};
    array(local_index) += amount;
  }
};

template void add_if_in_range (ArrayNd<FTYPE,NDIM> &array,  INDICIES(int nx, int ny, int nz), const int *istart, const int *iend, FTYPE amount);


void efield_p3dfft(field &Efield, FArrayND &rho)
{
  // Added by Glenn for debugging
  //  printf("Proc %d: rho(%d,%d,%d) = %f\n", mpi_rank,
  //     0,0,0, rho(INDICIES(0,0,0)));


  /* Local variables */

  FArrayND_ranged &phi = Efield.phi_rho;  // This is done for convenience only

  // phi  and ksqrinv must have some padding in the last dimension to 
  // allow fftw to perform an inplace transform

  int nsize[3]={nx*nsubdomains,ny,nz}; // This has to be 3D regardless and nz=1 if 2D

  // P3DFFT specific variables
  static int memory_dims[3], istart[3], iend[3], isize[3], fstart[3], fend[3], fsize[3], conf;
  // I need to break ny into pieces

  static const int ny2=(ny/subdomain.np);
  if (ny2 < 1) terminate(-1, "efield_p3dfft error: ny/nsubdomains less than 1.");
  static int p3dim[2]={subdomain.np,nsubdomains} ;
  static int local_size[3];

  // Define Arrays for the decomposed 3D transform data d1 and the operator array, k2inv.
  static ArrayNd<FTYPE,NDIM>  k2inv;
  static unsigned char op_f[]="fft", op_b[]="tff";

#ifdef write_dbg_p3d
  static ofstream diag2;
#endif

  int nx_global=nx*nsubdomains;
  
  static bool first_entry=true;
  if (first_entry) {
    first_entry=false;

#ifdef write_dbg_p3d
    //output this for diagnostic purposes (use 1 file per rank)
    std::stringstream fname2;
    fname2 << "k2_mpi" <<mpi_rank<<".txt";
    string cname2 = fname2.str();
    diag2.open(cname2.c_str());
#endif
    
    if (SUBDOMAINS_ADJACENT != false)
      terminate(-1,"efield_p3dfft error: Only works for SUBDOMAINS_ADJACENT=false\n");

    if ( ny%subdomain.np != 0 ) 
      terminate(-1,"efield_p3dfft error: ny not evenly divisible by number of processors per subdomain\n");


    // Set up p3dfft transforms
    

    /* Initialize P3DFFT */

    // Note that we are calling routines that assume a fortran ordering of data A(ix,iy,iz),
    // where ix is contiguous in memory.
    Cp3dfft_setup(p3dim, nsize[2], nsize[1], nsize[0],  MPI_Comm_c2f(MPI_COMM_WORLD),
		 nsize[2], nsize[1], nsize[0], 1, local_size);

#ifdef write_dbg_p3d
    diag2 << "p3dfft_setup memsize output:" << 
      local_size[0] << " " << local_size[1] << " " << local_size[2] << "\n";
#endif

    conf = 1; // Indicates a forward transform
    Cp3dfft_get_dims(istart,iend,isize,conf);
#ifdef write_dbg_p3d
    diag2 << " Cp3dfft_get_dims memsize output (mpi="<<mpi_rank<<"):" ;
    for (int i=0; i<NDIM;i++) diag2 <<" ("<<istart[i]<<" "<<iend[i]<<" "<<isize[i] <<"), ";
    diag2 << " " << conf << "\n";
#endif
    
    /* Get dimensions for output array - complex numbers, Z-pencil shape.
       Stride-1 dimension could be X or Z, depending on how the library 
       was compiled (stride1 option) */
    conf = 2;
    Cp3dfft_get_dims(fstart,fend,fsize,conf);
#ifdef write_dbg_p3d
    diag2 << " Cp3dfft_get_dims memsize output transformed (mpi="<<mpi_rank<<"):" ;
    for (int i=0; i<NDIM;i++) diag2 <<" ("<<fstart[i]<<" "<<fend[i]<<" "<<fsize[i]<<"), ";
    diag2 <<" " << conf << "\n";
#endif

    // Now we must declare arrays that have the appropriate size:
    // Since we're communicating with Fortran routines, we still want z continuous in memory
    // but p3dfft will interprit k2inv(x,y,z) as k2inv(z,y,x).  We'll call this k2inv_p
    
    k2inv=ArrayNd<FTYPE,NDIM> (INDICIES(local_size[2],local_size[1],local_size[0]));
    // Define a complex array containing the data (to make it easier to operate on
    // Assuming stride1 defined, the array will be x,y,z
    const int local_size_c[3]={fsize[2],fsize[1],fsize[0]};
    ArrayNd<std::complex<FTYPE>,NDIM> k2inv_c(INDICIES(local_size_c[0],local_size_c[1],local_size_c[2]),
					      (std::complex<FTYPE>*) k2inv.address());

    // The operator we need is just the fft of the finite difference opperator
    const int is2[3] = {istart[2]-1,istart[1]-1,istart[0]-1};  //Handle fortran indicies which start at 1
    const int ie2[3] = {iend[2]-1,iend[1]-1,iend[0]-1};  //Handle fortran indicies which start at 1

    add_if_in_range(k2inv,INDICIES(0,0,0),is2,ie2, -2/Sqr(dx));
    add_if_in_range(k2inv,INDICIES(1,0,0),is2,ie2, 1/Sqr(dx));
    add_if_in_range(k2inv,INDICIES(nsize[0]-1,0,0),is2,ie2, 1/Sqr(dx));
    
    add_if_in_range(k2inv,INDICIES(0,nsize[1]-1,0),is2,ie2, 1/Sqr(dy));
    add_if_in_range(k2inv,INDICIES(0,0,0),is2,ie2, -2/Sqr(dy));
    add_if_in_range(k2inv,INDICIES(0,1,0),is2,ie2, 1/Sqr(dy));

#if NDIM == 3
    add_if_in_range(k2inv,INDICIES(0,0,nsize[2]-1),is2,ie2, 1/Sqr(dz));
    add_if_in_range(k2inv,INDICIES(0,0,0),is2,ie2, -2/Sqr(dz));
    add_if_in_range(k2inv,INDICIES(0,0,1),is2,ie2, 1/Sqr(dz));
#endif

#ifdef write_dbg_p3d
    for (int ix = 0; ix < isize[2]; ++ix) 
      for (int iy = 0; iy < isize[1]; ++iy)
	for (int iz = 0; iz < isize[0]; iz++) {
	  diag2 <<"ksqr(pr="<< mpi_rank <<", "<< ix <<", "<< iy <<", "<< iz 
		 <<") = "<< k2inv(INDICIES(ix,iy,iz)) 
		 << " globalF: ("<< ix+istart[2] << "," << iy+istart[1] <<","<< iz+istart[0]<<")\n";
	  }
#endif

    // Transform into K space
    Cp3dfft_ftran_r2c(k2inv.address(),k2inv.address(),op_f);

    /*Note that the stride1 defined algorithm has the transpose so
      In Fortran it's K2invTF(nx,ny/M2,(nz+2)/2M1) which in C will appear as 
      K2invT((nz+2)/2M1,ny/M2,nx) where nx is contiguous in memory. */

#ifdef write_dbg_p3d
    for (int ix = 0; ix < isize[2]; ++ix) 
      for (int iy = 0; iy < isize[1]; ++iy)
	for (int iz = 0; iz < isize[0]; iz++) {
	  diag2 <<"ksqrinv(pr="<< mpi_rank <<", "<< ix <<", "<< iy <<", "<< iz 
		 <<") = "<< k2inv(INDICIES(ix,iy,iz)) 
		 << " globalF: ("<< ix+istart[2] << "," << iy+istart[1] <<","<< iz+istart[0]<<")\n";
	  }

    diag2 << "addresses:" << k2inv.address() <<", " <<  k2inv_c.address() <<"\n";
    for (int iz = 0; iz < fsize[2]; iz++) 
      for (int iy = 0; iy < fsize[1]; ++iy)
	for (int ix = 0; ix < fsize[0]; ix++) {
	  diag2 <<"ksqrinvC(pr="<< mpi_rank <<", "<< iz <<", "<< iy <<", "<< ix 
		<<") = ("<< k2inv_c(INDICIES(iz,iy,ix)) 
		 << ") globalFC: ("<< iz+istart[0] << "," << iy+istart[1] <<","<< ix+istart[2]<<")\n";
	  }
#endif
  
    // k2inv is now a k^2 like operator. Invert, element by element
    // and make the complex = real part for rapid multiplication later.
    for (int iz = 0; iz < fsize[2]; iz++) 
      for (int iy = 0; iy < fsize[1]; ++iy)
	for (int ix = 0; ix < fsize[0]; ix++) {
	  k2inv_c(INDICIES(iz,iy,ix))=1.0/k2inv_c(INDICIES(iz,iy,ix));
	}

    // Set DC component to zero:
    if (mpi_rank==0) k2inv_c(INDICIES(0,0,0))=std::complex<FTYPE>(0.0,0.0);

    // Incorporate the epsilon coefficient and the n_elements factor into k2inv:
    
    k2inv /= (-eps*nx_global)*ny*nz;

    //To avoid taking an fft of ne twice, I'll do the high-frequency filtering
    // in this routine:

    // Calculate the damping coefficients

      if (fsteep != 0 && fwidth != 0) {
	FTYPE ksqr, damp, kx, ky, kz; 
	FTYPE kxbase=fwidth/(nx_global/2.), 
	  kybase=fwidth/(ny/2.), kzbase=fwidth/(nz/2.);

	for (int ikx=0; ikx<fsize[2]; ikx++) {
	  int ikx_global=ikx+fstart[2]-1; // The -1 is because p3dfft has indices starting at 1 but we use indices starting at 0
	  if (ikx_global <= nx_global/2) kx=kxbase*ikx_global;
	  else kx=kxbase*(-(nx_global-ikx_global));

	  for (int iky=0; iky<fsize[1]; iky++) {
	    int iky_global=iky+fstart[1]-1; // The -1 is because p3dfft has indices starting at 1 but we use indices starting at 0
	    if (iky_global <= ny/2) ky=kybase*iky_global;
	    else ky=kybase*(-(ny-iky_global));

	    for (int ikz=0; ikz<fsize[0]; ikz++) {
	      int ikz_global=ikz+fstart[0]-1; // The -1 is because p3dfft has indices starting at 1 but we use indices starting at 0
	      if (ikz_global <= nz/2) kz=kzbase*ikz_global;
	      else kz=kzbase*(-(nz-ikz_global));

	      ksqr=kx*kx+ky*ky+kz*kz;
	      damp=exp(-1*pow(ksqr,fsteep));
	      if (ndim == 1) {
		if (ikx_global <= nx_global/2) {
		  k2inv_c(INDICIES(ikx,0,0)) *= damp;// Note the function of INDICIES
		}
	      }
	      else if (ndim == 2) {
		if (iky_global <= ny/2) {
		  k2inv_c(INDICIES(ikx,iky,0)) *= damp;// Note the function of INDICIES
		}
	      } else {
		if (ikz <= nz/2) {
		  k2inv_c(INDICIES(ikx,iky,ikz)) *= damp;
		}
	      }
	    }
	  }
	}
      }

    /* Add in later...
      // For some applications it is useful to eliminate 
      // parallel electric fields:
      if (no_parallel_efields && k2inv.x_start==0) {
	int ikx=0;
	for (int iky=0; iky<k2inv.ny_transpose; iky++) 
	  for (int ikz=0; ikz<nz; ikz++) 
	    k2inv_c(INDICIES(iky,ikx,ikz))=0.0;
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
	  for (int iky=0; iky<k2inv.ny_transpose; iky++) {
	    int iky_global=iky+k2inv.y_start_transpose;
	    if (iky_global <= ny/2) ky=kybase*iky_global;
	    else ky=kybase*(-(ny-iky_global));
	    for (int ikz=0; ikz<nz; ikz++) {
	      if (ikz <= nz/2) kz=kzbase*ikz;
	      else kz=kzbase*(-(nz-ikz));
	      ksqr=kx*kx+ky*ky+kz*kz;
	      if (ksqr > 0.) {
		if (kill_modes_ang_rad > acos(kz/sqrt(ksqr)) ) {
		  if (ikz <= nz/2) {
		    k2inv.cdata(INDICIES(iky,ikx,ikz)) = 0.;
		  }
		}
	      }
	    }
	  }
	}
      }

    */

#ifdef write_dbg_p3d
    //output this for diagnostic purposes (use 1 file per rank)
    ofstream diag;
    std::stringstream fname;
    fname << "k2inv_mpi" <<mpi_rank<<".txt";
    string cname = fname.str();
    diag.open(cname.c_str());

#endif
      // Since k2inv is a real number, and in order to enable fast multiplication,
      // we will make the real and imag parts of k2inv be equal.

      for (int iz = 0; iz < fsize[2]; iz++) 
      	for (int iy = 0; iy < fsize[1]; ++iy)
	  for (int ix = 0; ix < fsize[0]; ++ix) {
	    k2inv_c(INDICIES(iz,iy,ix))=
	      std::complex<FTYPE>(k2inv_c(INDICIES(iz,iy,ix)).real(),
				  k2inv_c(INDICIES(iz,iy,ix)).real());
	    //output this for diagnostic purposes
#ifdef write_dbg_p3d
	    diag <<"ksqrinv(pr="<< mpi_rank <<", "<< iz <<", "<< iy <<", "<< ix 
		 <<") = "<< k2inv_c(INDICIES(iz,iy,ix)) 
		 << "global: ("<< iz+fstart[2] << "," << iy+fstart[1] <<","<< ix+fstart[0]<<")\n";

#endif
	  }

#ifdef write_dbg_p3d
      diag.close();
#endif

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
  static FArrayND phi_trans(INDICIES(local_size[2],local_size[1],local_size[0]));

  // Define a complex array containing the data (to make it easier to operate on
  // Assuming stride1 defined, the array will be x,y,z
  const int local_size_c[3]={fsize[2],fsize[1],fsize[0]};
  static ArrayNd<std::complex<FTYPE>,NDIM> phi_trans_c(INDICIES(local_size_c[0],local_size_c[1],local_size_c[2]),
						(std::complex<FTYPE>*) phi_trans.address());

    // phi_trans(ix,iy,iz)  is organized like this to have z contiguous in memory,
    //In fortran it appears as phi_trans(iz,iy,ix)  but that means that 
    // the same array in my C++ code is k2inv(ix,iy,iz) 
    phi_trans=0.0;
    // We have to map the array from eppic's simpler 1D decomposition to the current 2D decomposition
    // So for nx=64, nsubdomains=2, ny= 64, nz=32 we will have on all processors
    // the local fortran organized array isize={z=32,y=8,x=64} (C reverses x and z)
    // Processor 0 will have (C) ix 0->63, iy 0->7, iz 0->31
    // Processor 1 will have (C) ix 0->63, iy 8->15, iz 0->31)

    for (int ix = 0; ix < isize[2]; ++ix) {
      int ix_global=ix+(istart[2]-1);  // This is the global ix value
      int ix2=ix_global-subdomain.id_number*nx;  // Hopefully, it only hits this in the appropriate domain.
      for (int iy = 0; iy < isize[1]; ++iy) {
	int iy2=iy+(istart[1]-1);  // This is the global and local rho iy value
	for (int iz = 0; iz < isize[0]; iz++) {
#ifdef write_dbg_p3d
	  diag2 <<"rho:"<< ix << "," << iy << "," << iz << ",  " << ix2 << "," << iy2;
#endif
	  phi_trans(INDICIES(ix,iy,iz))=rho(INDICIES(ix2,iy2,iz));
#ifdef write_dbg_p3d
	  diag2 << ", " << rho(INDICIES(ix2,iy2,iz)) << "\n";
#endif
	}}}

#ifdef write_dbg_p3d
    diag2 <<" Copied into phi_trans... \n";
#endif

    // Transform into k space
    Cp3dfft_ftran_r2c(phi_trans.address(),phi_trans.address(),op_f);
#ifdef write_dbg_p3d
    diag2 <<" Transformed... \n";
#endif

    ArrayNd<FTYPE,NDIM>  k2inv_copy=k2inv;
    ArrayNd<std::complex<FTYPE>,NDIM> 
      k2inv_c(INDICIES(local_size_c[0],local_size_c[1],local_size_c[2]),
	      (std::complex<FTYPE>*) k2inv_copy.address());
    k2inv_copy=k2inv;
#ifdef write_dbg_p3d
    for (int ix = 0; ix < fsize[2]; ++ix) 
      for (int iy = 0; iy < fsize[1]; ++iy)
	for (int iz = 0; iz < fsize[0]; iz++) {
	  diag2 <<"phi_trans:"<< ix << "," << iy << "," << iz ;
	  diag2 << ", " << phi_trans_c(INDICIES(ix,iy,iz)) << " x " << k2inv_c(INDICIES(ix,iy,iz)) << "\n";
	}
#endif

    // Calculate phi_trans.transform/ (-k^2 * eps*nx*ny*nz)
    // Test  k2inv=1./(nx_global*ny*nz);
    phi_trans *= k2inv;

#ifdef write_dbg_p3d
    diag2 <<" Multiplied...REVERSED X,Y for comparison with efield below:\n";
    for (int iy = 0; iy < fsize[1]; ++iy)
      for (int ix = 0; ix < fsize[2]; ++ix) 
	for (int iz = 0; iz < fsize[0]; iz++) {
	  diag2 <<"phiFFT(pr="<< mpi_rank<< ", "  << iy << ", " << ix << ", "<< iz ;
	  diag2 << ") = " << phi_trans_c(INDICIES(ix,iy,iz)) << "\n";
	}
#endif

    // Inverse fft transform:
    Cp3dfft_btran_c2r(phi_trans.address(),phi_trans.address(),op_b);

#ifdef write_dbg_p3d
    diag2 <<" Backwards Transformed... \n";
    for (int ix = 0; ix < isize[2]; ++ix) {
      for (int iy = 0; iy < isize[1]; ++iy) {
	for (int iz = 0; iz < isize[0]; iz++) {
	  diag2 <<"phi("<< ix << "," << iy << "," << iz ;
	  diag2 << ")=" << phi_trans(INDICIES(ix,iy,iz)) << "\n";
	}}}
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

    // Added by Glenn
    /*
    printf("Proc %d before ALLGATHER:\n", mpi_rank);
    for (int ix =0; ix<2; ix++){
      for (int iy=0; iy<2; iy++){
        for (int iz=0; iz<2; iz++){
          printf("Proc %d: phi(%d,%d,%d) = %f, phi_trans(%d,%d,%d) = %f\n", 
                 mpi_rank,ix,iy,iz,phi(INDICIES(ix,iy,iz)),
                 ix, iy, iz, phi_trans(INDICIES(ix,iy,iz)));
        }
      }
    }
    */

    // This needs to be done accross processors
    int phi_size_noguards=nx*ny;
    int phi_trans_size=isize[1];  // isize[1] is the y size
#if NDIM == 3
    phi_size_noguards *= nz;
    phi_trans_size *= isize[0]; // isize[0] is the z size
#endif
    // We have to do this one x row at a a time because otherwise, 
    // it get's block-transposed in memory.  There may be a faster way...
    for (int ix = 0; ix < isize[2]; ++ix) {
      int ifirst[]={INDICIES(ix,0,0)};
      int mpi_err = 
	MPI_Allgather(phi_trans.address(ifirst),phi_trans_size,MPI_FTYPE,
		      phi.address(ifirst),phi_trans_size,MPI_FTYPE,subdomain.internal_comm);
      if ( mpi_err != MPI_SUCCESS) 
	mpi_error(mpi_rank, mpi_err,"MPI_Allgather call in efield_p3dfft");
    }


    /*    for (int ix = 0; ix < isize[2]; ++ix) {
      int ix_global=ix+(istart[2]-1);  // This is the global ix value
      int ix2=ix_global-subdomain.id_number*nx;  // Hopefully, it only hits this in the appropriate domain.
      for (int iy = 0; iy < isize[1]; ++iy) {
	int iy2=iy+(istart[1]-1);  // This is the global and local rho iy value
	for (int iz = 0; iz < isize[0]; iz++) {
	  diag2 <<"phi:"<< ix << "," << iy << "," << iz << ",  " << ix2 << "," << iy2;
	  phi(INDICIES(ix2,iy2,iz))= phi_trans(INDICIES(ix,iy,iz));
	  diag2 << ", " << phi(INDICIES(ix2,iy2,iz)) << "\n";
	  }}}*/

#ifdef write_dbg_p3d
    diag2 <<" Copied to phi... \n";
#endif
    // This section passes guard cell(s) from an array on one processor 
    // to another on another processor.  It only passes the last dimension(s).
    // nx_guard[0] indicates how many guard cells to be passed on the LHS
    // nx_guard[1] indicates how many guard cells to be passed on the RHS
    // This can be done before copying the rest of the array over.

    void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);
    pass_guards(phi, phix_guard_size);

    // ADDED BY GLENN FOR DEBUGGING
    /*
    FTYPE phi_sum = 0.0;
    for (int ix=0;ix<local_size[2];ix++){
      for (int iy=0;iy<local_size[1]; iy++){
        for (int iz=0; iz<local_size[0]; iz++){
          phi_sum += phi_trans(INDICIES(ix,iy,iz));
        }
      }
    }
    printf("Proc %d phi_trans_sum = %f\n", mpi_rank, phi_sum);
    printf("Proc %d after after pass_gaurds:\n", mpi_rank);
    for (int ix =0; ix<2; ix++){
      for (int iy=0; iy<2; iy++){
        for (int iz=0; iz<2; iz++){
          printf("Proc %d: phi(%d,%d,%d) = %f\n", 
                 mpi_rank,ix,iy,iz,phi(INDICIES(ix,iy,iz)));
        }
      }
    }
    */

#ifdef write_dbg_p3d
    diag2 <<" Passed Guards... \n";

    for (int ix = phi.start(0); ix < phi.size(0)+phi.start(0); ++ix) {
      for (int iy = phi.start(1); iy < phi.size(1)+phi.start(1); ++iy) {
	for (int iz = phi.start(2); iz < phi.size(2)+phi.start(2); iz++) {
	  diag2 <<"phi_all("<< ix << "," << iy << "," << iz ;
	  diag2 << ")=" << phi(INDICIES(ix,iy,iz)) << "\n";
	}}}
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

#ifdef write_dbg_p3d
    std::stringstream fname2;
    fname2 << "phi_it" << it << "_mpi" <<mpi_rank<<".txt";
    string cname2 = fname2.str();
    phi.output(cname2.c_str());
#endif

    // Let's time the amount of time in this routine:
    
    /*  
} //end "between codes" block
*/

    // ADDED BY GLENN FOR DEBUGGING
    /*
      phi_sum = 0.0;
    for (int ix=0;ix<local_size[2];ix++){
      for (int iy=0;iy<local_size[1]; iy++){
        for (int iz=0; iz<local_size[0]; iz++){
          phi_sum += phi_trans(INDICIES(ix,iy,iz));
        }
      }
    }
    printf("Proc %d at the end of efield_p3dfft:\n", mpi_rank);
    printf("Proc %d phi_trans_sum = %f\n", mpi_rank, phi_sum);
    for (int ix =0; ix<2; ix++){
      for (int iy=0; iy<2; iy++){
        for (int iz=0; iz<2; iz++){
          printf("Proc %d: phi(%d,%d,%d) = %f\n", 
                 mpi_rank,ix,iy,iz,phi(INDICIES(ix,iy,iz)));
        }
      }
    }
    */



}
#endif
#endif
