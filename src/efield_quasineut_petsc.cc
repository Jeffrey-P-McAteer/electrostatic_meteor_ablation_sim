#include <stdio.h>
#include <math.h>
#include <sys/times.h> //Clock related functions
#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif

#if USE_QN
#if HAVE_PETSC
#include <petscksp.h> // includes other header files, eg petscmat.h ...
#include <petscpc.h>

#include "efield_quasineut.h"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   This next subroutine deals with all the PETSc stuff; This is 
   separated out for reasons of modularity. 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// int efield_quasineut_petsc(field &Efield,FArrayND_ranged qden,FArrayND_ranged divG)
int efield_quasineut_petsc(field &Efield,FArrayND_ranged qden,
			   INDICIES(FArrayND_ranged Gx,
				    FArrayND_ranged Gy,
				    FArrayND_ranged Gz))
{

  if (SUBDOMAINS_ADJACENT != false)
    terminate(-1,
	      "efield_quasineut_petsc error: " \
	      "Only works for SUBDOMAINS_ADJACENT=false\n");

  /* Local functions */
  void mpi_check(int error_code);
  void write_local_bin(FArrayND_ranged array, const char *name, int *wxguards);
  void write_local_bin(FArrayND array, const char *name);
  PetscErrorCode getPetscLHS(Mat &, FArrayND_ranged &);
  PetscErrorCode getPetscRHS(Vec &, field &, FArrayND_ranged &, 
			     FArrayND_ranged &);
  PetscErrorCode getPetscRHS_flux(Vec &, field &, 
				  FArrayND_ranged &,
				  INDICIES(FArrayND_ranged &,
					   FArrayND_ranged &,
					   FArrayND_ranged));
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);

  int wxguards[] = {0,0}; // Diagnostics
  if (hybrid_diagnostic_subcycle > 0 && 
      it%nout*hybrid_diagnostic_subcycle == 0) {
    char wlb_name[32];
    sprintf(wlb_name,"qden-%06d-",it);
    write_local_bin(qden,wlb_name,wxguards);
    sprintf(wlb_name,"Gx-%06d-",it);
    write_local_bin(Gx,wlb_name,wxguards);
    sprintf(wlb_name,"Gy-%06d-",it);
    write_local_bin(Gy,wlb_name,wxguards);
  }

  /* Grab phi from Efield for convenience */
  FArrayND_ranged &phi = Efield.phi_rho;

  /* MPI and PETSc error codes */
  int merr=0;
  PetscErrorCode perr=0;
  PetscBool optChk=PETSC_FALSE;

  /* Local variables */
  Vec x,b;               // Solution and RHS vectors

  /* Initialize */
  // Create the RHS Vec (b)
  perr=VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(perr);
  perr=VecSetSizes(b,PETSC_DECIDE,global_length);CHKERRQ(perr);
  perr=VecSetFromOptions(b);CHKERRQ(perr);

  // Duplicate the solution Vec (x) from b. Does not copy entries.
  perr=VecDuplicate(b,&x);CHKERRQ(perr);

  /* Build RHS vector (b) */
  // perr=getPetscRHS(b,Efield,qden,divG);CHKERRQ(perr);
  perr=getPetscRHS_flux(b,Efield,qden,INDICIES(Gx,Gy,Gz));CHKERRQ(perr);

  /* Check RHS */
  if (hybrid_diagnostic_subcycle > 0 && 
      it%nout*hybrid_diagnostic_subcycle == 0) {
    char wlb_name[32];
    sprintf(wlb_name,"rhs-%06d-",it);
    FArrayND rhs=FArrayND(INDICIES(nx,ny,nz));
    PetscScalar *bLocal;
    PetscInt ix,iy,iz,ir,rowStart,rowEnd,irl;
    // int nxl = nx/subdomain.np;
    int nxl = nx/petsc_np;
    int xoffset = nxl*subdomain.rank;
    int xshift=nx*subdomain.id_number;
    perr=VecGetArray(b,&bLocal);CHKERRQ(perr);
    perr=VecGetOwnershipRange(b,&rowStart,&rowEnd);CHKERRQ(perr);
    for (ir=rowStart;ir<rowEnd;ir++) {
      irl = ir - rowStart;
      iy = irl%ny;
      iz = irl/(nxl*ny);
      ix = irl/ny - iz*nxl + xoffset;
      rhs(INDICIES(ix,iy,iz)) = bLocal[irl];
    }
    perr=VecRestoreArray(b,&bLocal);CHKERRQ(perr);
    write_local_bin(rhs,wlb_name);
  }
  // perr=VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(perr);

  /* Build LHS operator matrix (A) */
  perr=getPetscLHS(A,qden);CHKERRQ(perr);
  // perr=MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(perr);

  /* Deal with null space of LHS operator */
  if (nullspace_method == 0) {
    // Does the null space contain the constant vector?
    PetscBool hasCnst=PETSC_TRUE;
    // Dimension of null space, excluding constant vector
    PetscInt nullDim=0;
    // PETSc object that will contain the null space
    MatNullSpace nullSpace;
    // Test whether this is actually a valid null space of A
    PetscBool isNull;
    perr=MatNullSpaceCreate(PETSC_COMM_WORLD,hasCnst,nullDim,
			    NULL,&nullSpace);CHKERRQ(perr);
    perr=MatNullSpaceTest(nullSpace,A,&isNull);CHKERRQ(perr);
    // Set the nullspace if it is valid
    if (isNull) {
      perr=MatSetNullSpace(A,nullSpace);CHKERRQ(perr);
      perr=MatSetTransposeNullSpace(A,nullSpace);CHKERRQ(perr);
    } else {
      perr=PetscPrintf(PETSC_COMM_WORLD,
		       "WARNING: Null space is not valid at time step %d\n",
		       it);CHKERRQ(perr);
    }
    // Free memory
    perr=MatNullSpaceDestroy(&nullSpace);CHKERRQ(perr);
  }

  // Update the Krylov-method operators
  perr=KSPSetOperators(ksp,A,A);CHKERRQ(perr);
  // Solve the system
  perr=KSPSolve(ksp,b,x);CHKERRQ(perr);

  /* Check the residual */
  if (hybrid_diagnostic_subcycle > 0 && 
      it%nout*hybrid_diagnostic_subcycle == 0) {
    Vec y;
    PetscReal y_norm;
    PetscViewer viewer;
    perr = VecDuplicate(x,&y);
    perr = MatMult(A,x,y);CHKERRQ(perr);
    perr = VecAXPY(y,-1.0,b);CHKERRQ(perr);
    perr = VecNorm(y,NORM_1,&y_norm);CHKERRQ(perr);
    perr = PetscPrintf(PETSC_COMM_WORLD,
		       "||Ax-b||_1 = %f\n",
		       y_norm);CHKERRQ(perr);
    perr = VecNorm(y,NORM_2,&y_norm);CHKERRQ(perr);
    perr = PetscPrintf(PETSC_COMM_WORLD,
		       "||Ax-b||_2 = %f\n",
		       y_norm);CHKERRQ(perr);
    perr = VecNorm(y,NORM_INFINITY,&y_norm);CHKERRQ(perr);
    perr = PetscPrintf(PETSC_COMM_WORLD,
		       "||Ax-b||_inf = %f\n",
		       y_norm);CHKERRQ(perr);
    perr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Ax-b.txt",
  				&viewer);CHKERRQ(perr);
    perr = VecView(y,viewer);CHKERRQ(perr);
    perr = PetscViewerDestroy(&viewer);CHKERRQ(perr);
    perr = VecDestroy(&y);CHKERRQ(perr);
  }  

  /* Set up arrays for building phi */
  FArrayND_ranged work;
  int phi_size[]={INDICIES(nx+phix_guard_size[0]+phix_guard_size[1],ny,nz)};
  int phi_start[]={INDICIES(0-phix_guard_size[0],0,0)};
  work = FArrayND_ranged(phi_start,phi_size);
  PetscScalar *xLocal;
  PetscInt ix,iy,iz,ir,rowStart,rowEnd,irl;
  // int nxl = nx/subdomain.np;
  int nxl = nx/petsc_np;
  int xoffset = nxl*subdomain.rank;
  int xshift=nx*subdomain.id_number;

  /* Get local portion of PETSc solution */
  perr=VecGetArray(x,&xLocal);CHKERRQ(perr);
  perr=VecGetOwnershipRange(x,&rowStart,&rowEnd);CHKERRQ(perr);

  /* Write to equivalent range in phi.
     ---------------------------------
     This fills in the range of phi equivalent to the range of x that
     this processor owns. In general, that range is NOT nx*ny*nz. 

     For example, with 2 subdomains and 4 processes:
     In subdomain 0, where global ix = [0,nx*subdomains/2),
     process 0 locally fills ix = [0,nx/2), iy = [0,ny), iz = [0,nz) and 
     process 1 locally fills ix = [nx/2,nx), iy = [0,ny), iz = [0,nz).
     In subdomain 1, where global ix = [nx*subdomains/2,nx*nsubdomains),
     process 2 locally fills ix = [0,nx/2), iy = [0,ny), iz = [0,nz) and
     process 3 locally fills ix = [nx/2,nx), iy = [0,ny), iz = [0,nz);

     All other entries in a local phi array are 0.0.
  */
  phi = 0.0;
  for (ir=rowStart;ir<rowEnd;ir++) {
    irl = ir - rowStart;
    iy = irl%ny;
    iz = irl/(nxl*ny);
    ix = irl/ny - iz*nxl + xoffset;
    phi(INDICIES(ix,iy,iz)) = xLocal[irl];
  }
  perr=VecRestoreArray(x,&xLocal);CHKERRQ(perr);

  /* Stitch local copies of phi together. 
     ------------------------------------
     Since partially filled local phi arrays are zero at every point not
     owned by this processor, MPI_Allreduce can build the full local phi
     array for this subdomain simply by adding all subdomain copies. This
     is why this routine only works when SUBDOMAINS_ADJACENT=false.

     Using the example above: 
     Subdomain 0 contains phi on processes 0 and 1. Co-adding them gives
     |------|------|   |------|------|   |------|------|   |------|------|
     | phi0   0.0  | + | 0.0    phi1 | = | phi0   phi1 | = |     phi     |
     |------|------|   |------|------|   |------|------|   |------|------|
     where the final phi represents a local phi array compatible with the
     rest of EPPIC. 

     NB: There is no need to divide by number of processors 
     (e.g., subdomain.np or petsc_np).
  */
  // work = 0.0;
  // merr=MPI_Allreduce((void*)&(phi(INDICIES(0-phix_guard_size[0],0,0))),
  // 		     (void*)&(work(INDICIES(0-phix_guard_size[0],0,0))),
  // 		     phi.size(),
  // 		     MPI_FTYPE,MPI_SUM,
  // 		     subdomain.internal_comm);
  work = 0.0;
  merr=MPI_Allreduce((void*)&(phi(INDICIES(0-phix_guard_size[0],0,0))),
		     (void*)&(work(INDICIES(0-phix_guard_size[0],0,0))),
		     phi.size(),
		     MPI_FTYPE,MPI_SUM,
		     PETSC_COMM_WORLD);
  phi = work;

  // --> From Yann's code
  if((subdomain.id_number+1==nsubdomains) && (boundary_type == INJECT)) {
    phi(INDICIES(nx-1,ny-1,nz-1)) = 0.0;
  }

  // FTYPE phiavg_local = 0;
  // // for (int ix=0;ix<nx;ix++) {
  // //   for (int iy=0;iy<ny; iy++) {
  // //     for (int iz=0;iz<nz; iz++) {
  // // 	phiavg_local+= phi(INDICIES(ix,iy,iz));
  // //     }
  // //   }
  // // }
  // long r_avg = nx*ny*nz;
  // if (nullspace_method == custom) r_avg -= 1;
  // for (ir=0;ir<r_avg;ir++) phiavg_local += phi(ir);
  // phiavg_local /= nx*ny*nz*nsubdomains;
  // FTYPE phiavg=0;
  // if (nsubdomains>1) {
  //   MPI_Allreduce(&phiavg_local,&phiavg,1,MPI_FTYPE,MPI_SUM,
  // 		  subdomain.neighbor_comm);
  // } else phiavg = phiavg_local;
  // // if (boundary_type == PERIODIC) 
  // //   for(int ix=0;ix<nx;ix++)
  // //     for(int iy=0;iy<ny;iy++)
  // // 	for (int iz=0;iz<nz;iz++)
  // // 	  phi(INDICIES(ix,iy,iz)) -= phiavg;
  // if (boundary_type == PERIODIC) 
  //   for (ir=0;ir<r_avg;ir++) phi(ir) -= phiavg;

  /* pass guard cells to neighboring subdomains */
  pass_guards(phi,phix_guard_size);

  // --> From Yann's code
  if (boundary_type == INJECT) {
    if (subdomain.id_number == 0)
      for (int iy = 0;iy<ny;iy++) 
  	for (int iz=0;iz<nz;iz++) 
  	  phi(INDICIES(0,iy,iz)) = 0;
    if (subdomain.id_number == nsubdomains-1)
      for (int iy = 0;iy<ny;iy++) 
  	for (int iz=0;iz<nz;iz++) 
  	  phi(INDICIES(nx-1,iy,iz)) = 0;
  }

  // for (ix=0;ix<nx;ix++)
  //   for (iy=0;iy<ny;iy++)
  //     for (iz=0;iz<nz;iz++) {
  // 	perr=PetscSynchronizedPrintf(PETSC_COMM_WORLD,
  // 				     "[%d] phi(%d,%d,%d) = %f\n",
  // 				     mpi_rank,ix,iy,iz,
  // 				     phi(INDICIES(ix,iy,iz)));
  // 	perr=PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  //     }

  /* Check phi pointwise */
  if (check_solution != 0) {
    int check;
    int check_quasineut_solution(field &Efield,
				 FArrayND_ranged &den,
				 INDICIES(FArrayND_ranged &Gx,
					  FArrayND_ranged &Gy,
					  FArrayND_ranged &Gz),
				 FArrayND_ranged &phi,
				 int check_solution);
    check = check_quasineut_solution(Efield,qden,INDICIES(Gx,Gy,Gz),
				     phi,check_solution);
    if (check != 0) {
      printf("[%d] check_quasineut_solution returned %d at it = %d\n",
	     mpi_rank,check,it);
    }
  }

  if (hybrid_diagnostic_subcycle > 0 && 
      it%nout*hybrid_diagnostic_subcycle == 0) {
    char wlb_name[32];
    sprintf(wlb_name,"phi-%06d-",it);
    write_local_bin(phi,wlb_name,wxguards);
  }

  /* Destroy PETSc objects */
  perr=VecDestroy(&x);CHKERRQ(perr);
  perr=VecDestroy(&b);CHKERRQ(perr);

  /* Electric field = -Grad(phi) */
  if (Efield.Ex.size() == Efield.phi_rho.size()) {
    void gradient(FArrayND &phi, 
  		  INDICIES(FArrayND &Ex, FArrayND &Ey, FArrayND &Ez), 
  		  INDICIES(FTYPE dx, FTYPE dy, FTYPE dz), FTYPE scaler);
    gradient(phi, INDICIES(Efield.Ex, Efield.Ey, Efield.Ez), 
  	     INDICIES(dx, dy, dz), (FTYPE) -1.); 
    // Impose any External or driving E-field 
    //  (changes to phi must be periodic!)
    void externale(field &Efield);
    externale(Efield);
  }

  /* Reset MPI error handler */
  // MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL);

  return(0);

}

#else
/* If not using PETSc, here is a dummy routine.*/
int efield_quasineut_petsc(field &Efield,FArrayND_ranged qden,
			   INDICIES(FArrayND_ranged Gx,
				    FArrayND_ranged Gy,
				    FArrayND_ranged Gz))
{ return(0); }
#endif
#else
/* If not using QN, here is a dummy routine.*/
int efield_quasineut_petsc(field &Efield, FArrayND_ranged qden, FArrayND_ranged divG)
{ return(0); }
#endif
