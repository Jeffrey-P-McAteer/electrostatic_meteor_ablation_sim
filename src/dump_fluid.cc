
#include <stdio.h> 
#include <string.h> 
#include "eppic.h"
#include "eppic-mpi.h"

void dump_fluid(fluid& fspecie, FILE *fdump, int id){



  fwrite(&(fspecie.den(INDICIES(int(-denx_guard_size),0,0))),
	 sizeof(fspecie.den(INDICIES(int(-denx_guard_size),0,0))),
	 fspecie.den.length(),fdump);
  fwrite(&(fspecie.vx(INDICIES(int(-velx_guard_size),0,0))),
	 sizeof(fspecie.vx(INDICIES(int(-velx_guard_size),0,0))),
	 fspecie.vx.length(),fdump);
  if (vel_dim[id] >= 2) fwrite(&(fspecie.vy(INDICIES(int(-velx_guard_size),0,0))),
			       sizeof(fspecie.vy(INDICIES(int(-velx_guard_size),0,0))),
			       fspecie.vy.length(),fdump);
  
  if (vel_dim[id] == 3) fwrite(&(fspecie.vz(INDICIES(int(-velx_guard_size),0,0))),
			       sizeof(fspecie.vz(INDICIES(int(-velx_guard_size),0,0))),
			       fspecie.vz.length(),fdump);
  

}
