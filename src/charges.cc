/* Calculate the charge on the mesh given the particles and ion behavior 

   The integer, it, is a counter which enables the routine to store 
   densities from distributions which are subcycled */

#include "eppic.h"
#include "eppic-mpi.h"

void charges(FArrayND &rho, particle_dist *pic, fluid *fspecie, int it)
{
  // Let's time the amount of time in this routine:
  clock_t start_time=times(&times_buf);
  
  extern void density(FArrayND &, int, particle_dist *, fluid *, FTYPE);
  //  extern void kfilter(FArrayND &x, FTYPE s, FTYPE n);

  // Using a density workspace array saves cpu time and 
  // is required for mpi passing
  FArrayND den=rho;

  static FArrayND *den_store;

  static int first_entry=TRUE;
  if (first_entry == TRUE) {
    first_entry=FALSE;
    /* Set up density storage for subcycled arrays */
    if (subcycle.max() > 1) den_store=new FArrayND[ndist];
    for (int id=0; id<ndist; ++id) 
      if (subcycle(id) > 1) den_store[id]=FArrayND(den);
  } // END  if (first_entry == TRUE)

  /* Initialize the charge density, rho, to 0 - 
     Note that all DC charges will be eliminated */
  rho = 0.;

  for (int id=0; id<ndist; ++id) {
    // Make sure that the particle has charge to contribute
    if (method[id] >= 0) { // add fluids after mpisum of pic densities
      if (subcycle(id) <= 1) {
        // Density returns the charge density of each species.
        density(den, id, pic, fspecie, qd[id]);
        rho += den;
      } else { // This is a subcycled species:
        
        if (it%subcycle[id] == 0)  // The density needs calculating:
          density(den_store[id], id, pic, fspecie, qd[id]);
        rho += den_store[id];
      }
    }
  } // END for (int id=0; id<ndist; ++id) 

  // Let's time the amount of time in this routine:
  charge_time += times(&times_buf)-start_time;

  start_time=times(&times_buf);
  
#ifdef USE_MPI
  // Sum the arrays around the networked processes 
  {
    void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);

    // Use density as the target array
    den = 0;
    int mpi_err=MPI_Allreduce((void*) &(rho(INDICIES(0,0,0))), 
			      (void*) &(den(INDICIES(0,0,0))), rho.size(),
			      MPI_FTYPE, MPI_SUM,
                              subdomain.internal_comm);
    if ( mpi_err != MPI_SUCCESS) 
      mpi_error(mpi_rank, mpi_err,"Allreduce call in charges");
    // Divide over all subdomains to get average density
    den /= subdomain.np;
    rho=den;

    // now add fluids
    for (int id=0; id<ndist; ++id) {
      if (method[id] < 0) { 
        if (subcycle(id) <= 1) {
          // Density returns the charge density of each species.
          density(den, id, pic, fspecie, qd[id]);
          for (int ix=0;ix<nx+xguard_size;ix++)
            for (int iy=0; iy<ny; iy++)
              for (int iz=0;iz<nz; iz++) {
                rho(INDICIES(ix,iy,iz)) += den(INDICIES(ix,iy,iz));
              }
        } else { // This is a subcycled species:
          
          if (it%subcycle[id] == 0)  // The density needs calculating:
            density(den_store[id], id, pic, fspecie, qd[id]);
          rho += den_store[id];
        }
      }
    } // END for (int id=0; id<ndist; ++id) 
    
#ifdef USE_DOMAINS
    void pass_sum_guard(FArrayND &, int nx);
    pass_sum_guard(rho, nx);
#endif

  }
#endif
 
  /* all processors need to wait for first to finish it's work */
  MPI_Barrier(MPI_COMM_WORLD);
  // Let's time the amount of time in this routine:
  collect_time += times(&times_buf)-start_time;
} /* charges */
