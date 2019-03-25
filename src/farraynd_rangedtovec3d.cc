#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Convert 3D FArrayND_ranged object into a Vector
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* The logic of this if...else sequence has been reversed.
#ifndef HAVE_PETSC
*/
/* if you don't have PETSC, here's a dummy routine */
/* and a dummy Vec definition */
/*
typedef int Vec;

void FArrayND_rangedToVec3D(FArrayND_ranged &phi, Vec &phiVec, int nx, int ny, int nz, int n_exclude=-10)
{
  if (mpi_rank==0) {
    cout << "\tThis version of EPPIC was not configured with PETSC\n"
	 << "\tTherefore, FArrayND_rangedToVec3D is not available.\n"
	 << "\tPlease recompile with access to PETSC or change your \n"
	 << "\tinput options.";
    terminate(TRUE,"\n");
  }
}
#else
*/

#if HAVE_PETSC
#include "petscvec.h" 

void FArrayND_rangedToVec3D(FArrayND_ranged &phi, Vec &phiVec, int nx, int ny, int nz, int n_exclude=-10)
{
  /* 
     For now get Petsc_based full vector.
     In the future, consider only getting the local part of phi, and 
     manage the communication independantly.
     
     Consider adding this communication option to ArrayND class
  */

  PetscScalar *globalPhiVec;
  int ix,iy,iz;
  int rowStart, rowEnd;

  VecGetOwnershipRange( phiVec, &rowStart, &rowEnd );
  //  VecGetArray(phiVec,&globalPhiVec);
  for (int r=rowStart; r<rowEnd; r++) {
    if (r==n_exclude) continue;
    iz = r/(nx*ny);
    iy = (r%ny);
    ix = ((r/ny)%nx)-nx*subdomain.id_number;
    //    cout << "R: " << r << " IX: " << ix << " IY: " << iy << "\n";
    /*
    PetscSynchronizedPrintf( subdomain.neighbor_comm,
		"Copy,Domain: (%d,%d)->R: %d, IX: %d, IY: %d PHI:%g\n",
			     subdomain.rank,subdomain.id_number,r,ix,iy
			     ,phi(ix,iy));
    PetscSynchronizedFlush( subdomain.neighbor_comm );
    */
    VecSetValue( phiVec, r, phi(ix,iy,iz), INSERT_VALUES );
  }

  VecAssemblyBegin(phiVec);
  VecAssemblyEnd(phiVec);

}

#else
/* If not using PETSc, here is a dummy routine. */
typedef int Vec;

void FArrayND_rangedToVec3D(FArrayND_ranged &phi, Vec &phiVec, int nx, int ny, int nz, int n_exclude=-10)
{ return; }

#endif
