#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Convert vector into a 2D FArrayND object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* The logic of this if...else sequence has been reversed.
#ifndef HAVE_PETSC
*/
/* if you don't have PETSC, here's a dummy routine */
/* and a dummy Vec definition */
/*
typedef int Vec;

void FArrayNDFromVec2D(FArrayND &phi, Vec &phiVec, int nx, int ny)
{
  if (mpi_rank==0) {
    cout << "\tThis version of EPPIC was not configured with PETSC\n"
	 << "\tTherefore, FArrayNDFromVec2D is not available. Please\n"
	 << "\trecompile with access to PETSC or change your input\n"
	 << "\toptions.";
    terminate(TRUE,"\n");
  }
}
#else
*/

#if HAVE_PETSC
#include "petscvec.h"

void FArrayNDFromVec2D(FArrayND &phi, Vec &phiVec, int nx, int ny)
{
  /* 
     For now get Petsc_based full vector.
     In the future, consider only getting the local part of phi, and 
     manage the communication independantly.
     
     Consider adding this communication option to ArrayND class
  */

  PetscScalar *globalPhiVec;
  int ix,iy,r;


  VecGetArray(phiVec,&globalPhiVec);
  r=0;
  for (ix=0;ix<nx;ix++)
    {
      for (iy=0;iy<ny;iy++)
	{
	  //  r = iy+(ix+nx*subdomain.id_number)*ny;
          // changed by glenn to allow complex petsc
#ifdef PETSC_USE_COMPLEX         
          phi(ix,iy)=PetscRealPart(globalPhiVec[r]);
#else
          phi(ix,iy)=globalPhiVec[r];
#endif
	  r++;
	}
    }
      
  VecRestoreArray(phiVec,&globalPhiVec);



}

#else
/* If not using PETSc, here is a dummy routine. */
typedef int Vec;

void FArrayNDFromVec2D(FArrayND &phi, Vec &phiVec, int nx, int ny)
{ return; }

#endif
