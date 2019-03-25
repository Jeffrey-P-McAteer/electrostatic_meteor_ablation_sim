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

void rejectArrayNd(particle_dist &dist, int id, int &iran, 
		   FArrayND den_arraynd)
{
  
    FTYPE ran3(int *);

  FTYPE n0per=0;//perturbation on n0
  FTYPE n0per_max_domain=0;//perturbation on n0
  FTYPE normal_factor=0;
  FTYPE global_max=0;
  FTYPE local_max=0,tmp_global_max=0;
  FTYPE den_value(FArrayND &den_arraynd, PTYPE x,PTYPE y, PTYPE z );

  cout << den_arraynd;


  for (int ix=0; ix<nx;++ix) {
    for (int iy=0; iy<ny; ++iy) {
      for (int iz=0; iz<nz; ++iz) {
	tmp_global_max=den_value(den_arraynd,ix,iy,iz);//(INDICIES(ix,iy,iz));den_arraynd(INDICIES(ix,iy,iz));
	
	n0per += den_value(den_arraynd,ix,iy,iz);//(INDICIES(ix,iy,iz));
	if (tmp_global_max > local_max) local_max = tmp_global_max;
      }
    }
  }

  if (nsubdomains==1) {
    n0per_max_domain = n0per;
    global_max = local_max;
  } else { /* find which domain has max */
    MPI_Allreduce(&n0per,&n0per_max_domain,
		  1,MPI_FTYPE,MPI_MAX,subdomain.neighbor_comm);
    MPI_Allreduce(&local_max,&global_max,
		  1,MPI_FTYPE,MPI_MAX,subdomain.neighbor_comm);
  }

  printf("Rank = %d; N0per = %g; n0per_max_domain = %g\n",
	 mpi_rank,n0per,n0per_max_domain);
  printf("Rank = %d; local_max = %g; global_max = %g\n",
	 mpi_rank,local_max,global_max);

  /*
    use this to find how many particles are needed 
    particle number is scaled so that the max domain has npd particles
  */
  int np_new = int(n0per*dist.np/n0per_max_domain);
  /* 
     See if too many particles are needed, and wont fit in allocated arrays
     This should never happen 
  */
  
  if (np_new >= dist.x.size()) terminate(-1,"rejectArrayNd attempted to place too many particles");
  if (np_new <= 1) {
    if (mpi_rank<nsubdomains) 
      cout << "WARNING: RejectND placing very few particles for subdoimain " <<
	subdomain.id_number << "!!!!\n";
  };
  // find normalization constant, use to compare to random number
  // only over this domain; by adjusting np to np_new, we will get the 
  // correct dist over all the domains

  normal_factor = 1/local_max;
  //  while (i<npd[id]) {
  dist.np=0;
  // need to work on normalization:
  // dist.n0avg*=n0per_max_domain/FTYPE((nx*ny*nz))/(global_max);
  //  n0peak[id] = dist.n0avg;
  if (mpi_rank == 0) {
    // printf("; WARNING: N0avg for PIC adjusted!\n");
    printf("\tn0avgd%1d = %g\n",id,dist.n0avg);
    printf("\tnp%1d = %d\n",id,np_new);
    printf("\tn0per_max_domain%1d = %g\n",id,n0per_max_domain);
    printf("\tglobal_max%1d = %g\n",id,global_max);
  }
    //cout << "Old N0: " << n0peak[id] <<
    //	       "New N0: " << 
    //		       dist.n0avg << "\n" << "Global max: " << global_max << "\n";
  
  int i=0;
  while(i<np_new) {
    PTYPE x, y, z;
    x = ran3(&iran)*(nx);
    if (ndim >= 2) y = ran3(&iran) * (ny);
    if (ndim == 3) z = ran3(&iran) * (nz);
    
    // Keep this particle only when:
    if (den_value(den_arraynd,x,y,z)*normal_factor > ran3(&iran)) {
      dist.x[dist.np]=x;
      if (ndim >= 2) dist.y[dist.np]=y;
      if (ndim == 3) dist.z[dist.np]=z;
      dist.np++;
      i++;
    }
  }

  
}


FTYPE den_value(FArrayND &den_arraynd, PTYPE x,PTYPE y, PTYPE z ) {
  FTYPE avg_value=0;
  int ix= static_cast<int>(floor(x));
  int iy= static_cast<int>(floor(y));
  int iz= static_cast<int>(floor(z));
  int npts=0;
  int ixp1=ix+1,iyp1=iy+1,izp1=iz+1;
  FTYPE dmx=x-ix,dmy=y-iy,dmz=z-iz;
  FTYPE dpx=ixp1-x,dpy=iyp1-y,dpz=izp1-z;
  if(iyp1>ny-1) iyp1-=ny;
  if(izp1>nz-1) izp1-=nz;

  /* x-1 */
  /* -- y-1 */
  /* -- -- z-1 */
  avg_value += den_arraynd(INDICIES(ix,iy,iz))*dmx*dmy*dmz; npts++;
  
  /* -- -- z+1 */
  avg_value += den_arraynd(INDICIES(ix,iy,izp1))*dmx*dmy*dpz; npts++;

  /* -- y+1 */
  /* -- -- z-1 */
  avg_value += den_arraynd(INDICIES(ix,iyp1,iz))*dmx*dpy*dmz; npts++;
  
  /* -- -- z+1 */
  avg_value += den_arraynd(INDICIES(ix,iyp1,izp1))*dmx*dpy*dpz; npts++;

  /* x+1 */
  /* -- y-1 */
  /* -- -- z-1 */
  avg_value += den_arraynd(INDICIES(ixp1,iy,iz))*dpx*dmy*dmz;   if(iyp1>ny-1) iyp1-=ny;
  if(izp1>nz-1) izp1-=nz;
npts++;
  
  /* -- -- z+1 */
  avg_value += den_arraynd(INDICIES(ixp1,iy,izp1))*dpx*dmy*dpz; npts++;

  /* -- y+1 */
  /* -- -- z-1 */
  
  avg_value += den_arraynd(INDICIES(ixp1,iyp1,iz))*dpx*dpy*dmz; npts++;
  /* -- -- z+1 */
  avg_value += den_arraynd(INDICIES(ixp1,iyp1,iz))*dpx*dpy*dpz; npts++;

  //avg_value /= npts;

  return avg_value;
}
