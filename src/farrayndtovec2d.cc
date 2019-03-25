#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif



/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Convert 2D FArrayND object into a Vector
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* The logic of this if...else sequence has been reversed.
#ifndef HAVE_PETSC
*/
/* if you don't have PETSC, here's a dummy routine */
/* and a dummy Vec definition */
/*
typedef int Vec;

void FArrayNDToVec2D(FArrayND &phi, Vec &phiVec, int nx, int ny)
{
  if (mpi_rank==0) {
    cout << "\tThis version of EPPIC was not configured with PETSC\n"
	 << "\tTherefore, FArrayNDToVec2D is not available. Please\n"
	 << "\trecompile with access to PETSC or change your input\n"
	 << "\toptions.";
    terminate(TRUE,"\n");
  }
}
#else
*/

#if HAVE_PETSC

#include "petscvec.h"


void FArrayNDToVec2D(FArrayND &phi, Vec &phiVec, int nx, int ny)
{
  /* 
     For now get Petsc_based full vector.
     In the future, consider only getting the local part of phi, and 
     manage the communication independantly.
     
     Consider adding this communication option to ArrayND class
  */

  PetscScalar *globalPhiVec;
  int ix,iy;
  int rowStart, rowEnd;

  VecGetOwnershipRange( phiVec, &rowStart, &rowEnd );
  //  VecGetArray(phiVec,&globalPhiVec);
  for (int r=rowStart; r<rowEnd; r++) {
    iy = (r % ny);
    ix = (r / ny)-nx*subdomain.id_number;
    /*
    PetscSynchronizedPrintf( subdomain.neighbor_comm,
		"Copy,Domain: (%d,%d)->R: %d, IX: %d, IY: %d Phi: %g\n",
			     subdomain.rank,subdomain.id_number,r,ix,iy,
			     phi(ix,iy));
    PetscSynchronizedFlush( subdomain.neighbor_comm );
    */
    VecSetValue( phiVec, r, phi(ix,iy), INSERT_VALUES );
  }
  VecAssemblyBegin(phiVec);
  VecAssemblyEnd(phiVec);

}

#else
/* If not using PETSc, here is a dummy routine. */
typedef int Vec;

void FArrayNDToVec2D(FArrayND &phi, Vec &phiVec, int nx, int ny)
{ return; }

#endif
