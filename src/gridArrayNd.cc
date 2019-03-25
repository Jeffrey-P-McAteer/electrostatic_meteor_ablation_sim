// Initialize the paticle quantities 
#include <sys/stat.h>
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
//    have the max particles, this is emulating an n0avg == max
//    We don't actually reset n0peak, we just use it in the scaling of the np
// 5) Find local max for this domain
//    This is used to scale den_func from 0 to 1, so sampling for this domain
//    is done correctly.
// 6) loop over np*new_factor, and sample using den_func/local_max

// EPPIC no longer uses this density scaling scheme. Background density is set by the input values 
// and not changed. --LKT 06/07/18

void gridArrayNd(particle_dist &dist, int id, int &iran, 
		   FArrayND den_arraynd)
{
  
  FTYPE ran3(int *);


  FTYPE local_max_den=0,global_max_den;
  FTYPE total_den=0,max_total_den=0;
  int ix_start=0;
  if ((subdomain.id_number == 0)&&(boundary_type[0] == INJECT)) ix_start = -1;

//   for (int ix=0; ix<nx;++ix) {
//     for (int iy=0; iy<ny; ++iy) {
//       for (int iz=0; iz<nz; ++iz) {
// 	tmp_global_max=den_arraynd(INDICIES(ix,iy,iz));
	
// 	n0per += den_arraynd(INDICIES(ix,iy,iz));
// 	if (tmp_global_max > local_max) local_max = tmp_global_max;
//       }
//     }
//   }


  total_den = den_arraynd.sum();
  local_max_den = den_arraynd.max();
  
  if (nsubdomains==1) {
    max_total_den = total_den;
    global_max_den = local_max_den;
  } else { /* find which domain has max */
    MPI_Allreduce(&total_den,&max_total_den,
		  1,MPI_FTYPE,MPI_MAX,subdomain.neighbor_comm);
    MPI_Allreduce(&local_max_den,&global_max_den,
		  1,MPI_FTYPE,MPI_MAX,subdomain.neighbor_comm);
  }

  int np_new = int(total_den*dist.np/max_total_den);
  // See if too many particles are needed, and wont fit in allocated arrays
  if (np_new >= dist.x.size()) 
    terminate(-1,"gridArrayND attempted to place too many particles");
  if (np_new <= 1) {
    if (mpi_rank<nsubdomains) 
      cout << "WARNING: gridArrayND placing very few particles for subdoimain " <<
	subdomain.id_number << "!!!!\n";
  };
  // find normalization constant, use to compare to random number
  // only over this domain; by adjusting np to np_new, we will get the 
  // correct dist over all the domains
  FTYPE normal_factor = 1./local_max_den;
  


  printf("Rank = %d; N0per = %g; n0per_max_domain = %g\n",
	 mpi_rank,total_den,max_total_den);
  printf("Rank = %d; local_max = %g; global_max = %g\n",
	 mpi_rank,local_max_den,global_max_den);


  dist.np=0;

  // normalization: n0avg is scaled by the fraction less than the global max.
  // in other words, if the domain with the highest point (global_max_den) was
  // also flat, then this normalization would be 1.
  // dist.n0avg*=max_total_den/
    FTYPE(((abs(ix_start)+nx+1)*(ny)*(nz)))/(global_max_den);

  if (mpi_rank == 0) {
    // printf("; WARNING: N0avg for PIC adjusted!\n");

    printf("\tn0avgd%1d = %g\n",id,dist.n0avg);
    printf("\tnp%1d = %d\n",id,np_new);
    printf("\tn0max_total_den%1d = %g\n",id,max_total_den);
    printf("\tglobal_max_den%1d = %g\n",id,global_max_den);


    char name[256];
    sprintf(name,"%seppic.i\0",outdir);
    
    FILE *fpout;
    if ((fpout=fopen(name,"a"))==NULL) {
      terminate(-1,"Cannot open copy of eppic.i for reading");
    };
    char line[130];
    sprintf(line,";; Changes during run:\n\nn0d%1d = %g\n",id,dist.n0avg);
    fputs(line,fpout);
    sprintf(line,"np%1d = %d\n",id,np_new);
    fputs(line,fpout);
    fclose(fpout);
    
  }


  
  int i=0;
  PTYPE nx_total = nx+abs(ix_start);
  while(i<np_new) {
    PTYPE x, y, z;
    x = ix_start+ran3(&iran)*nx;
    if (ndim >= 2) y = (ran3(&iran) * ny);
    if (ndim == 3) z = (ran3(&iran) * nz);

    int ix = static_cast<int>(x)%nx;
    int iy = static_cast<int>(y)%ny;
    int ixp1 = ix+1;
    int iyp1 = (iy+1)%ny;
#if NDIM==3
    int iz = static_cast<int>(z)%nz;
#endif
    
    // Keep this particle only when:
    FTYPE curr_den = 
      (y-iy)   *(
		 (x-ix)*den_arraynd(INDICIES(ix,iy,iz))+
		 (1+ix-x)*den_arraynd(INDICIES(ixp1,iy,iz))
		 )
      +(1+iy-y)*(
		 (x-ix)*den_arraynd(INDICIES(ix,iyp1,iz))+
		 (1+ix-x)*den_arraynd(INDICIES(ixp1,iyp1,iz))
		 );
      
    if (curr_den*normal_factor > ran3(&iran)) {
      dist.x[dist.np]=x;
      if (ndim >= 2) dist.y[dist.np]=y;
      if (ndim == 3) dist.z[dist.np]=z;
      dist.np++;
      i++;
    }
  }
}



