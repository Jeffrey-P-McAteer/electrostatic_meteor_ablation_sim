// DUMP to allow restart after crash or end:

#include <stdio.h> 
#include <string.h> 
FILE* fopensafe(char* filename, char* mode);
FILE* fopensafe(char* filename, char* mode, unsigned long skip);
#include "eppic.h"
#include "eppic-mpi.h"
// #include "eppic-efield.h"

// intermediate array; needs to be external for restart purposes
extern fluid *f_old;

void dump(particle_dist *pic, fluid *fspecie, field &Efield, int it){

  char name[128];
  void dump_fluid(fluid& fspecie, FILE *fdump, int id);
// intermediate array; needs to be external for restart purposes
  //  fluid *f_old;
  if (restart_nonlocal) {
    sprintf(name,"%s%srestart%d.rst",outdir,"restart/",mpi_rank);  
  } else {
    sprintf(name,"%srestart%d.rst","restart/",mpi_rank);
  }
  if (mpi_rank == 0) {
    printf(" Making restart files at t=%6d \n",it);
    //    printf ("%s\n",name); 
  }
  FILE *fdump = fopensafe(name,"wb");
  fwrite(&it,sizeof(it),1,fdump);
  for (int id=0; id<ndist; ++id) {
    if (method[id] >= 0) {
      fwrite(&(pic[id].n0avg),sizeof(pic[id].n0avg),1,fdump);
      fwrite(&(pic[id].np),sizeof(pic[id].np),1,fdump);
      fwrite(&(pic[id].x[0]),sizeof(pic[id].x[0]),pic[id].np,fdump);
      if (ndim >= 2) 
      fwrite(&(pic[id].y[0]),sizeof(pic[id].y[0]),pic[id].np,fdump);
      fwrite(&(pic[id].vx[0]),sizeof(pic[id].vx[0]),pic[id].np,fdump);
      if (ndim >= 2) 
      fwrite(&(pic[id].vy[0]),sizeof(pic[id].vy[0]),pic[id].np,fdump);
      if (vel_dim[id]==3) fwrite(&(pic[id].vz[0]),sizeof(pic[id].vz[0]),
				 pic[id].np,fdump);
      if (ndim==3) 
	fwrite(&(pic[id].z[0]),sizeof(pic[id].z[0]),pic[id].np,fdump);
    } else { // method[id] < 0, ie fluid
      if (subdomain.rank == subdomain.root) {
	if (method[id] != -4) {
	  dump_fluid(fspecie[id],fdump, id);
	  dump_fluid(f_old[id],fdump, id);
	}
      }
    }
  }
  if (efield_algorithm == 2) {
    fwrite(&(Efield.phi_rho(INDICIES(int(-phix_guard_size[0]),0,0))),
    	   sizeof(Efield.phi_rho(INDICIES(int(-phix_guard_size[0]),0,0))),
    	   Efield.phi_rho.length(),fdump);
  }
  fclose(fdump);
}
