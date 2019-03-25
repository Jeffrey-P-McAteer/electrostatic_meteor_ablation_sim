#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif
  

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get the scaler weighted flux on the mesh for a particular pic dist 
   and dimension dim. Sums over MPI copies, returns current in this 
   proc's domain. Modeled after output_flux.cc
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void pic_flux(particle_dist &pic,int dim,FArrayND &fluxDim,FTYPE scaler, 
	      int n_avg)
{
  // Zero fluxDim
  fluxDim=0.;
  
  // Define array limits
  int gn[]={INDICIES(nx/n_avg,ny/n_avg,nz/n_avg)};
  FTYPE gmin[]={INDICIES(0,0,0)};
  FTYPE gmax[]={INDICIES(nx,ny,nz)};
  // Generate a set of pointers to vectors of particles
  int wrap[]={INDICIES(nx/n_avg,ny/n_avg,nz/n_avg)};
  
  // These need to be adjusted if there are subdomains:
  if (nsubdomains > 1) {
    gn[0]+=1;
    wrap[0]+=1; // This eliminates the wrap
    gmax[0]+=1*n_avg; 
  }
  
  
  // Calculate total output gridsize
  FTYPE gridsize=nx/n_avg;
  for (int idim=1; idim<ndim; ++idim) gridsize *= (gmax[idim]-gmin[idim])/n_avg;
      

  // Define the array to sum onto
  ArrayNd<OTYPE,NDIM> flux(gn);
// #ifdef USE_MPI
//   // Define an mpi work array
//   ArrayNd<OTYPE,NDIM> flux_sum(gn);
// #endif
  
  flux=0.;

  FTYPE scale=scaler*gridsize*pic.n0avg/pic.np/dt;

  extern void gather_weight(ArrayNd<OTYPE, NDIM> &den, particle_dist &pic, 
			    PTYPEAVec &b, FTYPE *xmin, 
			    FTYPE *xmax, int *nmesh, int *wrap, FTYPE scale);

  switch(dim)
    {
    case 0: 
      scale*=dx;
      gather_weight(flux, pic, pic.vx, gmin, gmax, gn, wrap, scale);
      break;
    case 1:
      scale*=dy;
      gather_weight(flux, pic, pic.vy, gmin, gmax, gn, wrap, scale);
      break;
    case 2:
      scale*=dz;
      gather_weight(flux, pic, pic.vz, gmin, gmax, gn, wrap, scale);
      break;
    }


// #ifdef USE_MPI
//   // Sum the flux arrays accross all processors 
//   int mpi_err=MPI_Allreduce((void*)&flux(0),(void*)&flux_sum(0), flux.length(), 
// 			 MPI_OTYPE, MPI_SUM, subdomain.internal_comm);
//   if ( mpi_err != MPI_SUCCESS) 
//     mpi_error(mpi_rank, mpi_err,"Reduce flux call in output_fluxes failed");
//   flux = flux_sum;
//   flux /= subdomain.np;


// #ifdef USE_DOMAINS
//   if (nsubdomains>1) {
//     void pass_sum_guard(ArrayNd<OTYPE,NDIM> &, int);
//     pass_sum_guard(flux, nx);
//   }
// #endif
// #endif    


  // Copy flux to fluxDim
  // for(int ix=0;ix<nx;ix++)
  //   for (int iy=0;iy<ny;iy++)
  //     for (int iz=0;iz<nz;iz++)
  // 	fluxDim(INDICIES(ix,iy,iz)) += flux(INDICIES(ix,iy,iz));
  for(int ix=0;ix<nx;ix++)
    for (int iy=0;iy<ny;iy++)
      for (int iz=0;iz<nz;iz++)
  	fluxDim(INDICIES(ix,iy,iz)) = flux(INDICIES(ix,iy,iz));

  
}
