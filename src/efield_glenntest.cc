/*
This is a play routine to test various functions.
 */
#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif

//#include "classes/realfftwmpi_transpose.h"
//typedef realfftwmpi<FTYPE,NDIM-1> FArrayND_fftw;
#include "realfftwmpi_many.h"
typedef realfftwmpi_many<FTYPE,NDIM-1> rfftw_many;
#include "tridiag.h"

// indicies for data, real array
// Flip the x and z dimensions used in efield_inject_tridiag_solve
#if NDIM==2
#define RFFTW_INDICIES(x,y,z) y, x
#define RFFTW_C_INDICIES(x,y,z) y, x
#elif NDIM==3
//#define RFFTW_INDICIES(x,y,z) y, z, x
//#define RFFTW_C_INDICIES(x,y,z) z, y, x
#define RFFTW_INDICIES(x,y,z) x, y, z
#define RFFTW_C_INDICIES(x,y,z) y, x, z
#endif

// Tridiagonal solver (z direction), periodic in x, y

void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);

void efield_glenntest(field &Efield, FArrayND &rho)
{
  //if (mpi_rank == 0) printf("Entering efield_glenntest.cc.\n");
  int domain_xmin = nx*subdomain.id_number;
  int nx_global = nx*nsubdomains;
  int nzpg = nz+zguard_size;
  int nypg = ny+yguard_size;
  static bool first_solve = TRUE;
  index_type array_size[]={nx*nsubdomains, nypg};
  static rfftw_many working_array = rfftw_many(subdomain.neighbor_comm,
                                               array_size, nzpg);
  
  // Initialize vectors to be used for tridiagonal solver
  static ArrayNd<FTYPE,2> d_m;

  int iy_start = working_array.x_start;
  int iy_end = iy_start+working_array.local_nx;
  int ystart = 0;
  int yend = nypg;
  int xstart = working_array.y_start_transpose;
  int xend = xstart+working_array.ny_transpose;
  
  int n_solve = nx_global*working_array.ny_transpose;
  FArrayND_ranged &phi = Efield.phi_rho;  // This is done for convenience only

  // If this is the first time entering this function, set up the working arrays
  if (first_solve){
    first_solve = FALSE;

    // Set up the d_m array (from eqs 2 and 3 in Birdsall pg 319)
    // use max(1,working_array_ny_transpose) in case ny_transpose = 0
    d_m = ArrayNd<FTYPE,2>({nx_global, max(1,working_array.ny_transpose)});
    for (int ix=0; ix<nx_global; ix++){
      for (int iy=working_array.y_start_transpose; 
           iy<working_array.y_start_transpose+working_array.ny_transpose; iy++){
        d_m(ix,iy-working_array.y_start_transpose) = -2.0 
          -4.0*Sqr((dz/dx)*sin(M_PI* ((FTYPE)ix)/nx_global))
          -4.0*Sqr((dz/dy)*sin(M_PI* (FTYPE)iy/nypg));
      }
    }
  } // END first_solve


  // Zero the DC if you want
  FTYPE rho_dc=0;
  if (Efield.zero_dc) {
    for(int ix=0;ix<nx;ix++) 
      for(int iy=0;iy<nypg;iy++) 
	for(int iz=0;iz<nzpg;iz++) 
	  rho_dc+=rho(INDICIES(ix,iy,iz));
    rho_dc/=(nx*nypg*nzpg);

    if (nsubdomains>1)  {
      FTYPE rho_dc_tmp=rho_dc;
      int mpi_err=MPI_Allreduce(&rho_dc_tmp,&rho_dc,1,MPI_FTYPE,MPI_SUM,
                                subdomain.neighbor_comm);
      if ( mpi_err != MPI_SUCCESS)
        mpi_error(mpi_rank, mpi_err, "Allreduce call failed in efield_glenntest");
      rho_dc/=nsubdomains;
    }
  }

  // Create the working_array to do multiple 2D FFTs at once
  FTYPE factor = -dz*dz/Efield.eps;
  for (int ix=0; ix<nx; ix++){
    for (int iy=0; iy<nypg; iy++){
      for (int iz=0; iz<nzpg; iz++){
        working_array.data(RFFTW_INDICIES(working_array.x_start+ix,iy,iz)) = 
          (rho(INDICIES(ix,iy,iz)) - rho_dc)*factor;
        //Testing how the fft works...
        //working_array.data(RFFTW_INDICIES(working_array.x_start+ix,iy,iz)) = 
        //sin(iz*2*M_PI/nz)+cos(iy*2*M_PI/ny)+
        //sin((ix+working_array.x_start)*2*M_PI/(nx*nsubdomains));
        //working_array.data(RFFTW_INDICIES(working_array.x_start+ix,iy,iz)) = 
        //working_array.x_start+ix+iy+iz - (nz+nx_global+ny-3.0)/2.0;
      }
    }
  }
  // Transform the working array (2D FFT of xy planes)
  working_array.transform();

  // Step 2: Solve n_solve tridiagonal systems
  // Use the Thomas Algorithm assuming a=c=1, b=constant=d_m(ix,iy)
  // Setup the tridiagonal system of equations
  // Do if statements first for speed
  if (field_boundary_type[2][0]==open && 
      field_boundary_type[2][1]==open){
    // Open-open solver...
    // Use boundary solution in Birdsall pg 320
    // Phi at end = r^-1*phi_nz-1
    for (int ix=0; ix<nx_global; ix++){
      for (int iy=0; iy<working_array.ny_transpose; iy++){
        // Set the first elements of the sweep
        // See if we are in the ix=0 and iy=0 situation
        FTYPE b = d_m(ix,iy);
        FTYPE r = -b + sqrt(b*b-1);
        FTYPE rinv = 1/r;
        tridiag_solver_complex_open_open(d_m(ix,iy), rinv,
                                         &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                         &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                         nzpg);

        /*
        if (ix==0 && iy+working_array.y_start_transpose==0){
          // Solve for phi_-1 and phi_0 using equation 19 on page 322 of Birdsall Langdon
          FTYPE phi_im1=0;
          FTYPE phi_i=0;
          for (int iz=0; iz<nzpg; iz++){
            phi_im1 += (iz+1)*working_array.cdata.data[iz].real();
            phi_i   += iz*working_array.cdata.data[iz].real();
          }
          phi_im1 *= factor/2;
          phi_i *= factor/2;
          // Sweep through all z indicies and solve for phi
          working_array.cdata.data[0] = phi_i;
          for (int iz=1; iz<nzpg; iz++){
            working_array.cdata.data[iz] += -d_m(0,0)*phi_i-phi_im1;
            phi_im1 = phi_i;
            phi_i = working_array.cdata.data[iz].real();
          }
        } // end if ix=iy+y_start_transpose==0
        else{
          // Use the usual equation for kx!=0, ky!=0
          tridiag_solver_complex_open_open(d_m(ix,iy), rinv,
                                           &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                           &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                           nzpg);
        }
        */
      }
    }
  } // open open boundary
  else if (field_boundary_type[2][0]==periodic){
    // Periodic solver...
    // Use a cyclic tridiagonal solver as in Numerical Recipes pg 79
    // alpha = 1, beta = 1
    FTYPE alpha = 1.0;
    FTYPE beta = 1.0;
    for (int ix=0; ix<nx_global; ix++){
      for (int iy=0; iy<working_array.ny_transpose; iy++){
        if (ix==0 && iy+working_array.y_start_transpose==0){
          // The usual periodic matrix will be singular.
          // To deal with this, eliminate the last row and set phi=0 at that point
          tridiag_solver_complex(d_m(ix,iy),
                                 &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                 &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                 nzpg-1);
          working_array.cdata.data[(ix+iy*nx_global+1)*nzpg-1] = 0.0;
        }
        else{
          tridiag_cyclic_solver_complex(d_m(ix,iy),
                                        &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                        &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                        alpha,beta,nzpg);
        }
      }
    }
  } // periodic boundary
  else if (field_boundary_type[2][0]==dirichlet && 
           field_boundary_type[2][1]==dirichlet){
    // Dirichlet-dirichlet solver..
    for (int ix=0; ix<nx_global; ix++){
      for (int iy=0; iy<working_array.ny_transpose; iy++){
        // The boundaries are set to phi=0 (phi(0)=0, phi(nz+zguard_size-1)=0)
        tridiag_solver_complex_dirichlet_dirichlet(field_boundary_values[2][0],
                                                   field_boundary_values[2][1],
                                                   d_m(ix,iy), 
                                                   &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                                   &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                                   nzpg);
      }
    }
  } // dirichlet-dirichlet
  else if (field_boundary_type[2][0]==neumann &&
           field_boundary_type[2][1]==neumann){
    for (int ix=0; ix<nx_global; ix++){
      for (int iy=0; iy<working_array.ny_transpose; iy++){
        
        // Neumann-Neumann solver...
        // The boundaries are set to dphi/dn = field_boundary_values[2][0/1]
        tridiag_solver_complex_neumann_neumann(field_boundary_values[2][0]*dz,
                                               field_boundary_values[2][1]*dz,
                                               d_m(ix,iy),
                                               &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                               &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                               nzpg);                                               
      }
    }
  } 
  else if (field_boundary_type[2][0]==dirichlet &&
           field_boundary_type[2][1]==neumann){
    // Dirichlet-Neumann solver...   
    for (int ix=0; ix<nx_global; ix++){
      for (int iy=0; iy<working_array.ny_transpose; iy++){
        tridiag_solver_complex_dirichlet_neumann(field_boundary_values[2][0],
                                                 field_boundary_values[2][1]*dz,
                                                 d_m(ix,iy),
                                                 &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                                 &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                                 nzpg);
      }
    }
  } 
  else if (field_boundary_type[2][0]==neumann &&
           field_boundary_type[2][1]==dirichlet){
    //Neumann-Dirichlet solver...
    for (int ix=0; ix<nx_global; ix++){
      for (int iy=0; iy<working_array.ny_transpose; iy++){
        tridiag_solver_complex_neumann_dirichlet(field_boundary_values[2][0],
                                                 field_boundary_values[2][1]*dz,
                                                 d_m(ix,iy),
                                                 &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                                 &working_array.cdata.data[(ix+iy*nx_global)*nzpg],
                                                 nzpg);
      }
    }
  }
  else{
    // UNKNOWN FIELD BOUNDARY TYPE!!
    printf("Solver not implemented yet for %d, %d boundaries\n",
           field_boundary_type[2][0], field_boundary_type[2][1]);
  }
  
  // Take inverse FFT
  working_array.invtransform();

  // Now put working_array data into phi
  FTYPE nfft = nx_global*nypg;
  for (int ix=0; ix<nx; ix++){
    for (int iy=0; iy<nypg; iy++){
      for (int iz=0; iz<nzpg; iz++){
        phi(INDICIES(ix,iy,iz)) = 
          working_array.data(RFFTW_INDICIES(working_array.x_start+ix,iy,iz))/nfft;
      }
    }
  }
  
  // pass the guard cells 
  pass_guards(phi, phix_guard_size);
}
