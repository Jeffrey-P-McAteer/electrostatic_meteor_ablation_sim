#include <stdio.h> 
#include <string.h> 
#include "eppic.h"
#include "eppic-mpi.h"


void restart_fluid(fluid &fspecie, FILE *fdump, int id) 
{
  
  fread(&(fspecie.den(INDICIES(int(-1*denx_guard_size),0,0))),
	sizeof(fspecie.den(INDICIES(int(-1*denx_guard_size),0,0))),
	fspecie.den.size(),fdump);
  fread(&(fspecie.vx(INDICIES(int(-1*velx_guard_size),0,0))),
	sizeof(fspecie.vx(INDICIES(int(-1*velx_guard_size),0,0))),
	fspecie.vx.size(),fdump);
  
  if (vel_dim[id] >= 2) fread(&(fspecie.vy(INDICIES(int(-1*velx_guard_size),0,0))),
			      sizeof(fspecie.vy(INDICIES(int(-1*velx_guard_size),0,0))),
			      fspecie.vy.size(),fdump);
  
  if (vel_dim[id] == 3) fread(&(fspecie.vz(INDICIES(int(-1*velx_guard_size),0,0))),
			      sizeof(fspecie.vz(INDICIES(int(-1*velx_guard_size),0,0))),
			      fspecie.vz.size(),fdump);
  
}
