#include <stdio.h> 
#include <string.h> 

#include "eppic.h"
#include "eppic-mpi.h"



void bcast_fluid (fluid &fspecie, int id) 
{

  // one-to-all; form subdomain rank
	
  MPI_Bcast(
	    &(fspecie.den(INDICIES(int(-1*denx_guard_size),0,0))),
	    fspecie.den.size(), 
	    MPI_FTYPE,
	    subdomain.root,
	    subdomain.internal_comm);
  MPI_Bcast(
	    &(fspecie.vx(INDICIES(int(-1*velx_guard_size),0,0))),
	    fspecie.vx.size(), 
	    MPI_FTYPE,
	    subdomain.root,
	    subdomain.internal_comm);
  if (vel_dim[id] >= 2) {
    MPI_Bcast(
	      &(fspecie.vy(INDICIES(int(-1*velx_guard_size),0,0)
			   )),
	      fspecie.vy.size(), 
	      MPI_FTYPE,
	      subdomain.root,
	      subdomain.internal_comm);
  }
  if (vel_dim[id] == 3) {
    MPI_Bcast(
	      &(fspecie.vz(INDICIES(int(-1*velx_guard_size),0,0)
			   )),
	      fspecie.vz.size(), 
	      MPI_FTYPE,
	      subdomain.root,
	      subdomain.internal_comm);
  }
}
