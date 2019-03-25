// Calculate the density of distribution id scaled by "scaler"

#include <math.h>
#include <stdio.h>
#include "eppic.h"
#include "eppic-mpi.h"

void pic_density(FArrayND &den,int id,particle_dist &pic,FTYPE scaler)
{
  
  extern void gather(FArrayND &, INDICIES(PTYPEAVec &, PTYPEAVec &, PTYPEAVec &), 
		     const int, FTYPE);
  extern void gather_domain(FArrayND &, INDICIES(PTYPEAVec &, PTYPEAVec &, PTYPEAVec &),
			    const int, const int, FTYPE);
  extern void gather_domain_inject(FArrayND &, INDICIES(PTYPEAVec &, PTYPEAVec &, PTYPEAVec &),
				   const int, const int, FTYPE);

  if (chargeon[id] == 0) {
    
#ifndef USE_DOMAINS
    gather(den, INDICIES(pic.x, pic.y, pic.z), pic.np, scaler*pic[id].n0avg);
#else
    int nxmin = 0;
    if (boundary_type[0] == INJECT && subdomain.id_number==0) {
      nxmin=-1;
      gather_domain_inject(den, INDICIES(pic.x, pic.y, pic.z), pic.np, nxmin,
			   scaler*pic.n0avg*nx*ny*nz/npd[id]);
    } else {
      gather_domain(den, INDICIES(pic.x, pic.y, pic.z), pic.np, nxmin,
      		    scaler*pic.n0avg*nx*ny*nz/npd[id] );

    }
#endif
  }

  
//   FArrayND den_sum = den;
// #ifdef USE_MPI

//   // Sum the flux arrays accross all processors 

//   int mpi_err=MPI_Allreduce((void*)&den(0),(void*)&den_sum(0), den.length(), 
// 			 MPI_FTYPE, MPI_SUM, subdomain.internal_comm);
//   if ( mpi_err != MPI_SUCCESS) 
//     mpi_error(mpi_rank, mpi_err,"Reduce density call in pic_density failed");
//   den = den_sum;
//   den /= subdomain.np;

#ifdef USE_DOMAINS
  // if (nsubdomains>1) {
  //   void pass_sum_guard(FArrayND &rho, int nx_local);
  //   pass_sum_guard(den, nx);
  // }
  // void pass_sum_guard(FArrayND &rho, int nx_local);
  // pass_sum_guard(den, nx);
#endif

// #endif    

  
}
