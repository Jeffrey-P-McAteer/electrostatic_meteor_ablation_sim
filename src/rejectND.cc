// Initialize the paticle quantities 
#include <stdio.h>
#include <math.h>
#include "eppic.h"
#include "eppic-mpi.h"

// Basic Formula
// 1) Find Maximum for whole (all domain) grid. 
// 2) Find the total area for particular domain
// 3) Area = Factor*n0peak*nx*ny*ny, solve for factor
//    Factor = 1 if flat distrubution
//    Factor will be used to scale number of particles
//    eg if Factor = 1, we would want np particles
//    Factor < 1, fewer particles, etc.
// 4) Use max from (1) to scale this: new_factor=factor/max
//    Since we don't want more then np particles in domain, and the max will
//    have the max particles, this is emulating an n0peak == max
//    We don't actually reset n0peak, we just use it in the scaling of the np
// 5) Find local max for this domain
//    This is used to scale den_func from 0 to 1, so sampling for this domain
//    is done correctly.
// 6) loop over np*new_factor, and sample using den_func/local_max

// EPPIC no longer uses this density scaling scheme. Background density is set by the input values 
// and not changed. --LKT 06/07/18

void rejectND(particle_dist &dist, int id, int &iran, 
		  PTYPE (den_func)(int id,INDICIES(PTYPE x,PTYPE y,PTYPE z)))
{
  
  FTYPE ran3(int *);
    
  FTYPE local_max_den=0,global_max_den;
  FTYPE total_den=0,max_total_den=0;
  FTYPE curr_denval=0;
  int ix_start=0;
  if ((subdomain.id_number == 0)&&(boundary_type[0] == INJECT)) ix_start = -1;
  for (int ix=ix_start;ix<=nx;++ix) {
    for (int iy=0; iy<=ny; ++iy) {
      for (int iz=0; iz<=nz; ++iz) {
	curr_denval=den_func(id,INDICIES(ix,iy,iz));
	
	total_den +=
	  (
	   // z+
	   // - y+
	   den_func(id,INDICIES(ix-0.5,iy+0.5,iz+0.5))+
	   den_func(id,INDICIES(ix+0.5,iy+0.5,iz+0.5))+
	     // - y-
	   den_func(id,INDICIES(ix-0.5,iy-0.5,iz+0.5))+
	   den_func(id,INDICIES(ix+0.5,iy-0.5,iz+0.5))+
	   // z-
	   // - y+
	   den_func(id,INDICIES(ix-0.5,iy+0.5,iz-0.5))+
	   den_func(id,INDICIES(ix+0.5,iy+0.5,iz-0.5))+
	   // - y-
	   den_func(id,INDICIES(ix-0.5,iy-0.5,iz-0.5))+
	   den_func(id,INDICIES(ix+0.5,iy-0.5,iz-0.5))
	   )/8.;
	
	if (curr_denval> local_max_den) local_max_den = curr_denval;

      }
    }
  }

  if (nsubdomains>1) {
    MPI_Allreduce(&total_den,&max_total_den,1,MPI_FTYPE,MPI_MAX,
		  subdomain.neighbor_comm);
    MPI_Allreduce(&local_max_den,&global_max_den,1,MPI_FTYPE,MPI_MAX,
		  subdomain.neighbor_comm);
  }
  else {
    max_total_den=total_den;
    global_max_den = local_max_den;
  }

  int np_new = int(total_den*dist.np/max_total_den);
  // See if too many particles are needed, and wont fit in allocated arrays
  if (np_new >= dist.x.size()) 
    terminate(-1,"rejectND attempted to place too many particles");
  if (np_new <= 1) {
    if (mpi_rank<nsubdomains) 
      cout << "WARNING: RejectND placing zero particles in subdoimain " <<
	subdomain.id_number << "!!!!\n";
  };
  // find normalization constant, use to compare to random number
  // only over this domain; by adjusting np to np_new, we will get the 
  // correct dist over all the domains
  FTYPE normal_factor = 1/local_max_den;

  dist.np=0;

  // normalization: n0avg is scaled by the fraction less than the global max.
  // in other words, if the domain with the highest point (global_max_den) was
  // also flat, then this normalization would be 1.
  // ---this normalization is not correct, and causes scaling issues for density and E-fields. 
  // ---commenting it out causes a 'particle crosses the entire domain' error at it=0.
  // ---LKT 01/10/17
  // ---max_total_den and global_max_den are wanted values calculated before particles
  // ---are placed. 
  // ---global_max_den = peak density + background density in units of n0
  // --- max_total_den = max wanted total density in a subdomain
  // ---LKT 01/17/17
  //  dist.n0avg*=max_total_den/FTYPE(((abs(ix_start)+nx+1)*(ny+1)*(nz+1)))/(global_max_den);

  if (mpi_rank == 0) {
    //    printf("; WARNING: N0avg for PIC adjusted!\n");
    printf("\tn0avgd%1d = %g\n",id,dist.n0avg);
    printf("\tnp%1d = %d\n",id,np_new);
    printf("\ttotal_den = %g\n",total_den);
    printf("\tn0max_total_den%1d = %g\n",id,max_total_den);
    printf("\tglobal_max_den%1d = %g\n",id,global_max_den);
  }

  int i=0;
  int nx_total = nx+abs(ix_start);
  while(i<np_new) {
    PTYPE x, y, z;
    x = ix_start+ran3(&iran)*nx_total;
    if (ndim >= 2) y = ran3(&iran) * ny;
    if (ndim == 3) z = ran3(&iran) * nz;

    // Keep this particle only when:
    if (den_func(id,INDICIES(x,y,z))*normal_factor > ran3(&iran)) {
      dist.x[dist.np]=x;
      if (ndim >= 2) dist.y[dist.np]=y;
      if (ndim == 3) dist.z[dist.np]=z;
      dist.np++;
      i++;
    }
  }
  
}
