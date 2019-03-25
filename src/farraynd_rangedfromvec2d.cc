#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif



/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Convert vector into a 2D FArrayND_ranged object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* The logic of this if...else sequence has been reversed.
#ifndef HAVE_PETSC
*/
/* if you don't have PETSC, here's a dummy routine */
/* and a dummy Vec definition */
/*
typedef int Vec;
void FArrayND_rangedFromVec2D(FArrayND_ranged &phi, Vec &phiVec, int nx, int ny, int n_exclude)
{
  if (mpi_rank==0) {
    cout << "\tThis version of EPPIC was not configured with PETSC\n"
	 << "\tTherefore, FArrayND_rangedFromVec2D is not available.\n"
	 << "\tPlease recompile with access to PETSC or change your \n"
	 << "\tinput options.";
    terminate(TRUE,"\n");
  }
}
#else
*/

#if HAVE_PETSC
#include "petscvec.h" 

void FArrayND_rangedFromVec2D(FArrayND_ranged &phi, Vec &phiVec, int nx, int ny, int n_exclude=-10)
{
  /* 
     For now get Petsc_based full vector.
     In the future, consider only getting the local part of phi, and 
     manage the communication independantly.
     
     Consider adding this communication option to ArrayND class
  */

  PetscScalar *globalPhiVec;
  int ix,iy,r,localr,xshift=nx*subdomain.id_number;
  int rowStart,rowEnd;

  VecGetArray(phiVec,&globalPhiVec);

  // PetscInt phiVecSize_g,phiVecSize_l;
  // VecGetSize(phiVec,&phiVecSize_g);
  // VecGetLocalSize(phiVec,&phiVecSize_l);
  // PetscSynchronizedPrintf(PETSC_COMM_WORLD,
  // 			  "[%d] phiVecSize_g = %d  phiVecSize_l = %d\n",
  // 			  mpi_rank,phiVecSize_g,phiVecSize_l);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

  /* if we are using copies, loop through each copy's local value */
  if (subdomain.np !=1) { 

    // FTYPEAVec tmpvec=FTYPEAVec(nx*ny/subdomain.np);
    FTYPEAVec tmpvec=FTYPEAVec(nx*ny);

    // PetscSynchronizedPrintf(PETSC_COMM_WORLD,
    // 			    "[%d] tmpvec.size = %d\n",
    // 			    mpi_rank,tmpvec.size());
    // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    for (int sendcopy=0; sendcopy<subdomain.np;sendcopy++) {
      VecGetOwnershipRange( phiVec, &rowStart, &rowEnd );

      // printf("[%d] sendcopy = %d  subdomain.rank = %d\n",
      // 	     mpi_rank,sendcopy,subdomain.rank);

      if (sendcopy==subdomain.rank){
	int ivec=0;
	for (int r=rowStart;r<rowEnd;r++){
	ix = r/ny;
	iy = r % ny;
        // Changed by Glenn to allow complex petsc numbers
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
	ix = r/ny-xshift;
	iy = r % ny;
	if (r==n_exclude) {
	  localr++;
	  continue;
	} else {
	  localr++;
	}
	phi(ix,iy) = tmpvec[localr];
      }
    }
  } else {

    r=0;
    for (ix=0;ix<nx;ix++){
      for (iy=0;iy<ny;iy++)
	{
	  if (r==n_exclude) continue;
	  //	  cout << ix << "\t" << iy << "\t" << globalPhiVec[r] << "\n";

	  // PetscSynchronizedPrintf(subdomain.neighbor_comm,
	  // 			       "[%d] (r,iy,iz) = (%03d,%03d,%03d)\n", 
	  // 			       mpi_rank,r,ix,iy);
	  // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
          
          // changed by glenn to allow petsc complex numbers
#ifdef PETSC_USE_COMPLEX         
          phi(ix,iy) = PetscRealPart(globalPhiVec[r]);
#else
	  phi(ix,iy) = globalPhiVec[r];
#endif
	  r++;
	}
    }
      
  }

  VecRestoreArray(phiVec,&globalPhiVec);

}

#else
/* If not using PETSc, here is a dummy routine. */
typedef int Vec;

void FArrayND_rangedFromVec2D(FArrayND_ranged &phi, Vec &phiVec, int nx, int ny, int n_exclude)
{ return; }

#endif
