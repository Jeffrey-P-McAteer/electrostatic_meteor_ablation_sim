 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Solves for phi, the electric potential, using Poisson's equation.
    
    We use the PETSc library, which allows for distributed solutions to
    linear (used in this case, ie Ax=b), non-linear, and matrix free 
    problems. Currently this will use the multigrid petsc method.
    
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif

#if HAVE_PETSC
#include <petscksp.h> // includes other header files, eg petscmat.h ...
#include <petscdmda.h>
#include <petscis.h>
#include <petscdm.h>

extern PetscErrorCode ComputeMatrixMG(KSP,Mat,Mat,void*);
extern PetscErrorCode ComputeRHSMG(KSP,Vec,void*);
extern PetscErrorCode ComputeInitialGuess(KSP,Vec,void*);

void efield_multigrid(field &Efield, FArrayND &rho)
{
  static KSP      ksp;
  static DM       da;
  static bool     first_solve = TRUE;

  FArrayND_ranged &phi = Efield.phi_rho;  // This is done for convenience only
  PetscReal       norm;
  PetscErrorCode  ierr;
  PetscInt        mx,my,mz,xm,ym,zm,xs,ys,zs;
  PetscScalar     Hx,Hy,Hz;
  PetscScalar     ***array;
  Vec             x,b;
  PetscInt        proc_petx, proc_pety, proc_petz, glob_petx, glob_pety, glob_petz;
  int             phi_temp_size;

  proc_petx = 1;
  proc_pety = mpi_np/nsubdomains;
  proc_petz = nsubdomains;
  glob_petx = nz;
  glob_pety = ny;
  glob_petz = nx*nsubdomains;
  phi_temp_size = glob_petx*glob_pety*glob_petz/(proc_petx*proc_pety*proc_petz);
  if (first_solve) {
    PetscPrintf(PETSC_COMM_WORLD,"Building PETSC solver: A dim: %dx%dx%d, procx,y,z: %d, %d, %d\n",
                glob_petx,glob_pety,glob_petz,proc_petx,proc_pety,proc_petz);

    PetscPrintf(PETSC_COMM_WORLD,"Creating KSP.\n");
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); 
    PetscPrintf(PETSC_COMM_WORLD,"Creating DMDA.\n");
    ierr = DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC,
                        DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                        DMDA_STENCIL_STAR,
                        glob_petx,glob_pety,glob_petz,proc_petx,proc_pety,proc_petz,
                        1,1,0,0,0,&da);
    PetscPrintf(PETSC_COMM_WORLD,"Importing options into DM.\n");
    ierr = DMSetFromOptions(da); 
    PetscPrintf(PETSC_COMM_WORLD,"Setting up DM.\n");
    ierr = DMSetUp(da);
    PetscPrintf(PETSC_COMM_WORLD,"Setting DM interoplation type.\n");
    ierr = DMDASetInterpolationType(da, DMDA_Q0); 
    PetscPrintf(PETSC_COMM_WORLD,"Attaching DM to KSP .\n");
    ierr = KSPSetDM(ksp,da); 
    
    //ierr = KSPSetComputeRHS(ksp,ComputeRHSMG,&rho);
    PetscPrintf(PETSC_COMM_WORLD,"Setting KSP RHS.\n");
    //ierr = KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,&phi);
    ierr = KSPSetComputeRHS(ksp,ComputeRHSMG,&rho);
    PetscPrintf(PETSC_COMM_WORLD,"Setting KSP A Matrix.\n");
    ierr = KSPSetComputeOperators(ksp,ComputeMatrixMG,&rho);
    PetscPrintf(PETSC_COMM_WORLD,"Importing options into KSP.\n");
    ierr = KSPSetFromOptions(ksp); 
    first_solve = FALSE;
    PetscPrintf(PETSC_COMM_WORLD,"Finished PETSc matrix setup.  Now doing the first KSP solve.\n");
  }


  //PetscPrintf(PETSC_COMM_WORLD,"Solving KSP.\n");
  ierr = KSPSolve(ksp,NULL,NULL);
  //PetscPrintf(PETSC_COMM_WORLD,"Extracting solution vector.\n");
  ierr = KSPGetSolution(ksp,&x); 
  //ierr = KSPGetRhs(ksp,&b); 
  //PetscPrintf(PETSC_COMM_WORLD, "RHS array:\n");
  //VecView(b, PETSC_VIEWER_STDOUT_WORLD);

  // See if we can store x inside of arrays
  ierr = DMDAGetInfo(da, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0); 
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 
  //PetscPrintf(PETSC_COMM_WORLD,"Building phi_temp of size %d,%d,%d\n", zm,ym,xm);
  //PetscPrintf(PETSC_COMM_WORLD,"Retreiving local array values.\n");
  ierr = DMDAVecGetArray(da, x, &array);    
  // create phi_temp array of dimensions (realx, realy, realz) = (petz, pety, petx)
  FArrayND phi_temp(INDICIES(zm, ym, xm));

  //PetscPrintf(PETSC_COMM_WORLD,"Putting values into phi_temp.\n");
  //int num_print = 5;
  //int counter = 0;
  for (int k=zs; k<zs+zm; k++) {
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        int phi_x = k-zs;
        int phi_y = j-ys;
        int phi_z = i-xs;
#ifdef PETSC_USE_COMPLEX
        phi_temp(INDICIES(phi_x,phi_y,phi_z)) = PetscRealPart(array[k][j][i]);
#else
        phi_temp(INDICIES(phi_x,phi_y,phi_z)) = array[k][j][i];
#endif
        /*
        if (counter++ <= num_print){
          printf("Proc %d: phi_temp(%d,%d,%d)=%f, array = %f rho(%d,%d,%d) = %f)\n", mpi_rank,
                 phi_x, phi_y, phi_z, phi_temp(INDICIES(phi_x,phi_y,phi_z)), array[k][j][i],
                 rho_x, rho_y, rho_z,rho(INDICIES(rho_x, rho_y, rho_z)));
        }
        */
      }
    }
  }
  //PetscPrintf(PETSC_COMM_WORLD,"Restoring solution vector.\n");
  ierr = DMDAVecRestoreArray(da, x, &array); 
  //PetscPrintf(PETSC_COMM_WORLD,"Assembling solution vector.\n");
  //ierr = VecAssemblyBegin(x);
  //ierr = VecAssemblyEnd(x);
  //ierr = VecDestroy(&x); 

  // Send the phi_temp vector to the correct processors

  /*
  int ifirst[] = {INDICIES(0,0,0)};
  int mpi_err = MPI_Allgather(phi_temp.address(ifirst),phi_temp_size,MPI_FTYPE,
                              phi.address(ifirst),phi_temp_size,MPI_FTYPE,
                              subdomain.internal_comm);
  */

  /*
  // Print out certain phi values before and after
  printf("Before allgather Proc %d: phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi_temp(%d,%d,%d)=%f\n",
         mpi_rank, 0,0,0, phi(INDICIES(0,0,0)), 
         0,1,0, phi(INDICIES(0,1,0)),
         0,glob_pety/proc_pety,0, phi(INDICIES(0,glob_pety/proc_pety,0)),
         nx-1,ny-1,nz-1, phi(INDICIES(nx-1,ny-1,nz-1)),
         0,0,0, phi_temp(INDICIES(0,0,0)));
  */
  
  //PetscPrintf(PETSC_COMM_WORLD,"Calling allgather.\n");
  for (int ix = 0; ix<nx; ++ix){
    int ifirst[] = {INDICIES(ix,0,0)};
    int mpi_err = MPI_Allgather(phi_temp.address(ifirst), phi_temp_size/nx, MPI_FTYPE,
                                phi.address(ifirst), phi_temp_size/nx, MPI_FTYPE,
                                subdomain.internal_comm);
    if ( mpi_err != MPI_SUCCESS) 
      mpi_error(mpi_rank, mpi_err,"MPI_Allgather call in efield_glenntest");
  }

  /*
  printf("After allgather Proc %d: phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi_temp(%d,%d,%d)=%f\n",
         mpi_rank, 0,0,0, phi(INDICIES(0,0,0)), 
         0,1,0, phi(INDICIES(0,1,0)),
         0,glob_pety/proc_pety,0, phi(INDICIES(0,glob_pety/proc_pety,0)),
         nx-1,ny-1,nz-1, phi(INDICIES(nx-1,ny-1,nz-1)),
         0,0,0, phi_temp(INDICIES(0,0,0)));
  */
  /*
  printf("Proc %d: before pass_guards phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f\n",
         mpi_rank, -1,0,0, phi(INDICIES(-1,0,0)), 0,0,0, phi(INDICIES(0,0,0)),
         nx-1,0,0, phi(INDICIES(nx-1,0,0)), nx,0,0, phi(INDICIES(nx,0,0)),
         nx+1,0,0, phi(INDICIES(nx+1,0,0)));
  */
 
  //PetscPrintf(PETSC_COMM_WORLD,"Passing guard cells.\n");
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);
  pass_guards(phi, phix_guard_size);

  /*
  printf("Proc %d: after pass_guards phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f, phi(%d,%d,%d)=%f\n",
         mpi_rank, -1,0,0, phi(INDICIES(-1,0,0)), 0,0,0, phi(INDICIES(0,0,0)),1,0,0, phi(INDICIES(1,0,0)),
         nx-1,0,0, phi(INDICIES(nx-1,0,0)), nx,0,0, phi(INDICIES(nx,0,0)),
         nx+1,0,0, phi(INDICIES(nx+1,0,0)));
  */
  /*
  // Error reporting stuff
  ierr = KSPGetOperators(ksp,NULL,&J);
  ierr = VecDuplicate(b,&r); 
  ierr = MatMult(J,x,r); 
  ierr = VecAXPY(r,-1.0,b); 
  ierr = VecNorm(r,NORM_2,&norm); 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Residual norm %g\n",(double)norm); 
  ierr = DMDAGetInfo(da, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0); 
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 
  Hx   = 1.0 / (PetscReal)(mx);
  Hy   = 1.0 / (PetscReal)(my);
  Hz   = 1.0 / (PetscReal)(mz);
  ierr = DMDAVecGetArray(da, x, &array); 
  for (int k=zs; k<zs+zm; k++) {
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        array[k][j][i] -=
          PetscCosScalar(2*PETSC_PI*(((PetscReal)i+0.5)*Hx))*
          PetscCosScalar(2*PETSC_PI*(((PetscReal)j+0.5)*Hy))*
          PetscCosScalar(2*PETSC_PI*(((PetscReal)k+0.5)*Hz));
      }
    }
  }
  ierr = DMDAVecRestoreArray(da, x, &array); 
  ierr = VecAssemblyBegin(x); 
  ierr = VecAssemblyEnd(x); 

  ierr = VecNorm(x,NORM_INFINITY,&norm); 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Error norm %g\n",(double)norm); 
  ierr = VecNorm(x,NORM_1,&norm); 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Error norm %g\n",(double)(norm/((PetscReal)(mx)*(PetscReal)(my)*(PetscReal)(mz)))); 
  ierr = VecNorm(x,NORM_2,&norm); 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Error norm %g\n",(double)(norm/((PetscReal)(mx)*(PetscReal)(my)*(PetscReal)(mz)))); 

  ierr = VecDestroy(&r); 
  ierr = KSPDestroy(&ksp); 
  ierr = DMDestroy(&da); 
  */
  //PetscPrintf(PETSC_COMM_WORLD,"Exiting efield_multigrid.\n");
}

PetscErrorCode ComputeRHSMG(KSP ksp,Vec b,void *ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
  PetscScalar    Hx,Hy,Hz;
  PetscScalar    ***array;
  DM             da;
  MatNullSpace   nullspace;
  FArrayND rho = *(FArrayND*)ctx;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&da);
  ierr = DMDAGetInfo(da, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0); 

  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 
  ierr = DMDAVecGetArray(da, b, &array);
  for (k=zs; k<zs+zm; k++) {
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        // See if we are at a boundary
        if (i==0 && field_boundary_type[2] !=PERIODIC_FIELD){
          if (field_boundary_type[2] == DIRICHLET){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps - field_boundary_values[2][0]/(dz*dz);
          }
          else if (field_boundary_type[2] == NEUMANN){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps + 2.0/(3*dz)*field_boundary_values[2][0];
          }
        }
        else if (i==mx-1 && field_boundary_type[2] !=PERIODIC_FIELD){
          if (field_boundary_type[2] == DIRICHLET){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps - field_boundary_values[2][1]/(dz*dz);
          }
          else if (field_boundary_type[2] == NEUMANN){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps + 2.0/(3*dz)*field_boundary_values[2][1];
          }
        }
        else if (j==0 && field_boundary_type[1] !=PERIODIC_FIELD){
          if (field_boundary_type[1] == DIRICHLET){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps - field_boundary_values[1][0]/(dy*dy);
          }
          else if (field_boundary_type[1] == NEUMANN){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps + 2.0/(3*dy)*field_boundary_values[1][0];
          }
        }
        else if (j==my-1 && field_boundary_type[1] !=PERIODIC_FIELD){
          if (field_boundary_type[1] == DIRICHLET){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps - field_boundary_values[1][1]/(dy*dy);
          }
          else if (field_boundary_type[1] == NEUMANN){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps + 2.0/(3*dy)*field_boundary_values[1][1];
          }
        }
        else if (k==0 && field_boundary_type[0] !=PERIODIC_FIELD){
          if (field_boundary_type[0] == DIRICHLET){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps - field_boundary_values[0][0]/(dx*dx);
          }
          else if (field_boundary_type[0] == NEUMANN){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps + 2.0/(3*dx)*field_boundary_values[0][0];
          }
        }
        else if (k==mz-1 && field_boundary_type[0] !=PERIODIC_FIELD){
          if (field_boundary_type[0] == DIRICHLET){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps - field_boundary_values[0][1]/(dx*dx);
          }
          else if (field_boundary_type[0] == NEUMANN){
            array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps + 2.0/(3*dx)*field_boundary_values[0][1];
          }
        }
        else{
          array[k][j][i] = -rho(INDICIES(k-zs,j,i-xs))/eps;
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(da, b, &array); 
  ierr = VecAssemblyBegin(b); 
  ierr = VecAssemblyEnd(b); 
  
  /* force right hand side to be consistent for singular matrix */
  /* note this is really a hack, normally the model would provide you with a consistent right handside */
  ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace); 
  ierr = MatNullSpaceRemove(nullspace,b); 
  ierr = MatNullSpaceDestroy(&nullspace); 
  PetscFunctionReturn(0);
}


PetscErrorCode ComputeMatrixMG(KSP ksp, Mat J,Mat jac, void *ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,num, numi, numj, numk;
  PetscScalar    v[7],Hx,Hy,Hz,HyHzdHx,HxHzdHy,HxHydHz;
  MatStencil     row, col[7];
  DM             da;
  MatNullSpace   nullspace;

  PetscFunctionBeginUser;
  ierr    = KSPGetDM(ksp,&da); 
  ierr    = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); 
  Hx = 1.0/(dz*dz);
  Hy = 1.0/(dy*dy);
  Hz = 1.0/(dx*dx);
  ierr    = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 

  // Build the matrix one row at a time...
  for (k=zs; k<zs+zm; k++) {
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        row.i = i; row.j = j; row.k = k;
        // See if we are not at a global domain boundary in x direction
        if (i==0 && j==0 && k == 0 && i==mx-1 && j==my-1 && k==mz-1) {
          // i,j,k is not at a global domain boundary
          // do the usual second order finite difference for f''
          v[0] = Hz;                   col[0].i = i;   col[0].j = j;   col[0].k = k-1;
          v[1] = Hy;                   col[1].i = i;   col[1].j = j-1; col[1].k = k;
          v[2] = Hx;                   col[2].i = i-1; col[2].j = j;   col[2].k = k;
          v[3] = -2.0*(Hx + Hy + Hz);  col[3].i = i;   col[3].j = j;   col[3].k = k;
          v[4] = Hx;                   col[4].i = i+1; col[4].j = j;   col[4].k = k;
          v[5] = Hy;                   col[5].i = i;   col[5].j = j+1; col[5].k = k;
          v[6] = Hz;                   col[6].i = i;   col[6].j = j;   col[6].k = k+1;
          ierr = MatSetValuesStencil(jac,1,&row,7,col,v,INSERT_VALUES);
        }
        else {
          // We are at a boundary
          num=1; v[0]=0;
          col[0].i = i; col[0].j=j; col[0].k=k;
          if (k==0 && (field_boundary_type[0] != PERIODIC_FIELD)){
            // We are at the k=0 boundary (really x=0 domain boundary bc petsc transpose)
            // Determine which boundary scheme we need to use
            if (field_boundary_type[0] == DIRICHLET){
              v[0] += -2.0*Hz;
              v[num] = Hz;
              col[num].i = i;
              col[num].j = j;
              col[num].k = k+1;
              num++;
            }
            else if (field_boundary_type[0] == NEUMANN){
              v[0] += -Hz*2.0/3.0;
              v[num] = Hz*2.0/3.0;
              col[num].i = i;
              col[num].j = j;
              col[num].k = k+1;
              num++;
            }
            else{
              printf("Error in efield_multigrid.  Boundary condition %d not supported!\n",
                     field_boundary_type[0]);
              terminate(-1, "Error!");
            }
          } // if k==0
          else if (k==mz-1 && (field_boundary_type[0] != PERIODIC_FIELD)){
            // We are at the k=mz-1 boundary (really x=nx*nsubdomains-1 boundary bc petsc transpose)
            // Determine which boundary scheme we need to use
            if (field_boundary_type[0] == DIRICHLET){
              v[0] += -2.0*Hz;
              v[num] = Hz;
              col[num].i = i;
              col[num].j = j;
              col[num].k = k-1;
              num++;
            }
            else if (field_boundary_type[0] == NEUMANN){
              v[0] += -Hz*2.0/3.0;
              v[num] = Hz*2.0/3.0;
              col[num].i = i;
              col[num].j = j;
              col[num].k = k-1;
              num++;
            }
            else{
              terminate(-1, "Error in efield_multigrid.  Boundary condition not supported!\n");
            }   
          } // if k==mz-1
          else{
            // Not at a k boundary and/or periodic
            v[0] += -2.0*Hz;
            v[num] = Hz;
            col[num].i = i;
            col[num].j = j;
            col[num].k = k-1;
            num++;
            v[num] = Hz;
            col[num].i = i;
            col[num].j = j;
            col[num].k = k+1;
            num++;
          } // else (k boundary)
          if (j==0 && field_boundary_type[1] != PERIODIC_FIELD){
            // We are at the j=0 boundary
            // Determine which boundary scheme we need to use
            if (field_boundary_type[1] == DIRICHLET){
              v[0] += -2.0*Hy;
              v[num] = Hy;
              col[num].i = i;
              col[num].j = j+1;
              col[num].k = k;
              num++;
            }
            else if (field_boundary_type[1] == NEUMANN){
              v[0] += -Hy*2.0/3.0;
              v[num] = Hy*2.0/3.0;
              col[num].i = i;
              col[num].j = j+1;
              col[num].k = k;
              num++;
            }
            else{
              terminate(-1, "Error in efield_multigrid.  Boundary condition not supported!\n");
            }
          } // if j==0
          else if (j==my-1 && field_boundary_type[1] != PERIODIC_FIELD){
            // We are at the j=my-1 boundary
            // Determine which boundary scheme we need to use
            if (field_boundary_type[1] == DIRICHLET){
              v[0] += -2.0*Hy;
              v[num] = Hy;
              col[num].i = i;
              col[num].j = j-1;
              col[num].k = k;
              num++;
              // No other inputs, take care of it with RHS     
            }
            else if (field_boundary_type[1] == NEUMANN){
              v[0] += -Hy*2.0/3.0;
              v[num] = Hy*2.0/3.0;
              col[num].i = i;
              col[num].j = j-1;
              col[num].k = k;
              num++;
            }
            else{
              terminate(-1, "Error in efield_multigrid.  Boundary condition not supported!\n");
            }            
          } // if j==my-1
          else{
            // Not at a j boundary and/or periodic
            v[0] += -2.0*Hy;
            v[num] = Hy;
            col[num].i = i;
            col[num].j = j-1;
            col[num].k = k;
            num++;
            v[num] = Hy;
            col[num].i = i;
            col[num].j = j+1;
            col[num].k = k;
            num++;
          } // else j boundary
          if (i==0 && field_boundary_type[2] != PERIODIC_FIELD){
            // We are at the i=0 boundary
            // Determine which boundary scheme we need to use
            if (field_boundary_type[2] == DIRICHLET){
              v[0] += -2.0*Hx;
              v[num] = Hx;
              col[num].i = i+1;
              col[num].j = j;
              col[num].k = k;
              num++;
            }
            else if (field_boundary_type[2] == NEUMANN){
              v[0] += -Hx*2.0/3.0;
              v[num] = Hx*2.0/3.0;
              col[num].i = i+1;
              col[num].j = j;
              col[num].k = k;
              num++;
            }
            else{
              terminate(-1, "Error in efield_multigrid.  Boundary condition not supported!\n");
            }
          } // if i==0
          else if (i==mx-1 && field_boundary_type[2] != PERIODIC_FIELD){
            // We are at the i=mx-1 boundary
            // Determine which boundary scheme we need to use
            if (field_boundary_type[2] == DIRICHLET){
              v[0] += -2.0*Hx;
              v[num] = Hx;
              col[num].i = i-1;
              col[num].j = j;
              col[num].k = k;
              num++;
              // No other inputs, take care of it with RHS     
            }
            else if (field_boundary_type[2] == NEUMANN){
              v[0] += -Hx*2.0/3.0;
              v[num] = Hx*2.0/3.0;
              col[num].i = i-1;
              col[num].j = j;
              col[num].k = k;
              num++;
            }
            else{
              terminate(-1, "Error in efield_multigrid.  Boundary condition not supported!\n");
            }            
          } // if i==mx
          else{
            // Not at an i boundary and/or periodic
            v[0] += -2.0*Hx;
            v[num] = Hx;
            col[num].i = i-1;
            col[num].j = j;
            col[num].k = k;
            num++;
            v[num] = Hx;
            col[num].i = i+1;
            col[num].j = j;
            col[num].k = k;
            num++;
          } // else i boundary
          ierr = MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);
        } //else (we are at at least one boundary)
      } // loop through x
    } // loop through y
  } // loop through z
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY); 
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY); 

  //MatView(jac, PETSC_VIEWER_STDOUT_WORLD);  
  ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace); 
  ierr = MatSetNullSpace(J,nullspace); 
  ierr = MatNullSpaceDestroy(&nullspace);
  PetscFunctionReturn(0);
}


PetscErrorCode ComputeInitialGuess(KSP ksp, Vec x, void *ctx){
  // Set the guess to the current phi
  //FArrayND_ranged phi = *(FArrayND_ranged*) ctx;
  FArrayND_ranged &phi = *(FArrayND_ranged*) ctx;
  PetscInt *global_coords;
  PetscScalar *values;
  PetscInt num_vals;
  PetscInt xs, ys, zs, xm, ym, zm, mx, my, mz;
  DM da;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscPrintf(PETSC_COMM_WORLD,"Inside initial guess vector.\n");
  
  ierr    = KSPGetDM(ksp,&da); 
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 
  ierr = DMDAGetInfo(da, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0); 
  PetscMalloc1(xm*ym*zm, &values);
  PetscMalloc1(xm*ym*zm, &global_coords);

  
  num_vals = 0;
  for (int k=zs; k<zs+zm; k++) {
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        global_coords[num_vals] = i*(my*mz)+j*(mz)+k;
        values[num_vals++] = phi(INDICIES(k-zs, j, i));
      }
    }
  }
  
  VecSetValues(x, num_vals, global_coords, values, INSERT_VALUES);
  ierr = VecAssemblyBegin(x);
  ierr = VecAssemblyEnd(x);
  
  /*ierr = DMDAVecGetArray(da, x, &array);  
  
  PetscPrintf(PETSC_COMM_WORLD,"array[%d][%d][%d]= %f, phi(%d,%d<5d)= %f\n",
              zs,ys,xs, array[zs][ys][xs], 0,ys,xs, phi(INDICIES(0,ys,xs)));

  for (int k=zs; k<zs+zm; k++) {
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        array[k][j][i] = phi(INDICIES(k-zs, j, i));
      }
    }
  }

  PetscPrintf(PETSC_COMM_WORLD,"Restoring initial guess vector.\n");
  ierr = DMDAVecRestoreArray(da, x, &array); 
  */
  
  //VecSet(x,0);
  PetscPrintf(PETSC_COMM_WORLD,"Leaving initial guess vector.\n");
  return(0);
}
#else
void efield_multigrid(field &Efield, FArrayND &rho){
  terminate(-1, "Error: efield_multigrid can only handle NDIM=3 and HAVE_PETSC=1");
}
#endif
