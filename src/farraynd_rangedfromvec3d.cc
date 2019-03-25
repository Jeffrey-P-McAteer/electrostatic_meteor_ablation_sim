#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif



/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Convert vector into a 3D FArrayND_ranged object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* The logic of this if...else sequence has been reversed.
#ifndef HAVE_PETSC
*/
/* if you don't have PETSC, here's a dummy routine */
/* and a dummy Vec definition */
/*
typedef int Vec;
void FArrayND_rangedFromVec3D(FArrayND_ranged &phi, Vec &phiVec, int nx, int ny, int nz, int n_exclude)
{
  if (mpi_rank==0) {
    cout << "\tThis version of EPPIC was not configured with PETSC\n"
	 << "\tTherefore, FArrayND_rangedFromVec3D is not available.\n"
	 << "\tPlease recompile with access to PETSC or change your \n"
	 << "\tinput options.";
    terminate(TRUE,"\n");
  }
}
#else
*/

#if HAVE_PETSC
#include "petscvec.h" 

void FArrayND_rangedFromVec3D(FArrayND_ranged &phi, Vec &phiVec, int nx, int ny, int nz, int n_exclude=-10)
{
  /* 
     For now get Petsc_based full vector.
     In the future, consider only getting the local part of phi, and 
     manage the communication independantly.
     
     Consider adding this communication option to ArrayND class
  */

  PetscScalar *globalPhiVec;
  int ix,iy,iz,r,localr,xshift=nx*subdomain.id_number;
  int rowStart,rowEnd;

  VecGetArray(phiVec,&globalPhiVec);

  /* if we are using copies, loop through each copy's local value */
  if (subdomain.np !=1) { 
    FTYPEAVec tmpvec=FTYPEAVec(nx*ny*nz/subdomain.np);
    for (int sendcopy=0; sendcopy<subdomain.np;sendcopy++) {
      VecGetOwnershipRange( phiVec, &rowStart, &rowEnd );
      if (sendcopy==subdomain.rank){
	int ivec=0;
	for (int r=rowStart;r<rowEnd;r++){
	  ix = (r/ny)%nx;
	  iy = r%ny;
	  iz = r/(nx*ny);
          // changed by glenn to allow petsc complex numbers
#ifdef PETSC_USE_COMPLEX         
          tmpvec[ivec]=PetscRealPart(globalPhiVec[ivec]);
#else
          tmpvec[ivec]=globalPhiVec[ivec];
#endif
	  ivec++;
	}
      } 
      /* send from sendcopy to others, then add to phi */
      MPI_Bcast(
		&tmpvec[0],
		tmpvec.size(), 
		MPI_FTYPE,
		sendcopy,
		subdomain.internal_comm);
      MPI_Bcast(
		&rowStart,
		1,
		MPI_INT,
		sendcopy,
		subdomain.internal_comm);
      MPI_Bcast(
		&rowEnd,
		1,
		MPI_INT,
		sendcopy,
		subdomain.internal_comm);
    
      localr=-1;
      for (r = rowStart;
	   r < rowEnd;
	   r++) {
	ix = ((r/ny)%nx)-xshift;
	iy = r%ny;
	iz = r/(nx*ny);
	if (r==n_exclude) {
	  localr++;
	  continue;
	} else {
	  localr++;
	}
	phi(ix,iy,iz) = tmpvec[localr];
      }
    }
  } else {
    r=0;
    for (ix=0;ix<nx;ix++){
      for (iy=0;iy<ny;iy++){
	for (iz=0;iz<nz;iz++){
	  if (r==n_exclude) continue;
	  //	  cout << ix << "\t" << iy << "\t" << globalPhiVec[r] << "\n";
#ifdef PETSC_USE_COMPLEX
          phi(ix,iy,iz) = PetscRealPart(globalPhiVec[r]);
#else
          phi(ix,iy,iz) = globalPhiVec[r];
#endif
	  r++;
	}
      }
    }
      
  }

  VecRestoreArray(phiVec,&globalPhiVec);

}

#else
/* If not using PETSc, here is a dummy routine. */
typedef int Vec;

void FArrayND_rangedFromVec3D(FArrayND_ranged &phi, Vec &phiVec, int nx, int ny, int nz, int n_exclude)
{ return; }

#endif
