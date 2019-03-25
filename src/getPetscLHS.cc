/* Build the LHS operator matrix (A) of the quasi-neutral, inertialess
   potential equation */

/* NOTES
   --Started 07May15 (may)
   --Including the line "-test_lhs" (without quotes) in the PETSc options
     file causes a centered-difference Laplacian to be used in place of the
     actual operator. (may)
*/

#if USE_QN
#include "eppic.h"
#include "efield_quasineut.h"

#if HAVE_PETSC
#include <petscksp.h>

PetscErrorCode getPetscLHS(Mat &A, FArrayND_ranged &den)
{

  PetscFunctionBegin;

  /* DEV: Create debug array */
  void write_local_bin(FArrayND_ranged array, const char *name, int *nx_ghost);
  void write_local_bin(FArrayND array, const char *name);
  FArrayND petscRHS=FArrayND(INDICIES(nx,ny,nz));

  int wxguards[] = {0,0}; // Diagnostics
  if (hybrid_diagnostic_subcycle > 0 && 
      it%nout*hybrid_diagnostic_subcycle == 0) {
    char wlb_name[32];
    sprintf(wlb_name,"qden_LHS-%06d-",it);
    write_local_bin(den,wlb_name,wxguards);
  }

  // Error tracking
  PetscErrorCode perr;

  // Grid and stencil variables
  PetscScalar values[stencil_size];
  PetscInt idx,index[stencil_size];
  PetscInt ir,rowStart,rowEnd,irl;
  int ix,iy,iz,ixp1,iyp1,izp1,ixm1,iym1,izm1;
  int ixg,iyg,izg,ixgp1,iygp1,izgp1,ixgm1,iygm1,izgm1;
  // int nxl = nx/subdomain.np;
  int nxl = nx/petsc_np;
  int xoffset = nxl*subdomain.rank;
  int xshift=nx*subdomain.id_number;
  int nxT=nx*nsubdomains;
  int nrows=nxT*ny;
  if (ndim == 3) nrows *= nz;
  int count;
  FTYPE value_diag = 1.0;
  FTYPE value_diagInv = 1.0;
  FTYPE dx_sqr = dx*dx;
  FTYPE dy_sqr = dy*dy;
  FTYPE dz_sqr = dz*dz;
  FTYPE dxInv = 1.0/dx,dyInv = 1.0/dy,dzInv = 1.0/dz; 
  FTYPE dx_sqrInv = 1.0/dx_sqr,dy_sqrInv = 1.0/dy_sqr,dz_sqrInv = 1.0/dz_sqr; 

  // Short-hand for density at grid points
  FTYPE den_xyz,den_xp1,den_yp1,den_zp1,den_xm1,den_ym1,den_zm1;

  // Physical quantities
  FTYPE Omega_e = Bz*fabs(qd[electron_dist])/md[electron_dist];
  FTYPE kappa = Omega_e/coll_rate[electron_dist];
  FTYPE c1 = kappa/(8*dx*dy);

  // Check if user has set "-test_lhs" in the PETSc options file
  PetscBool testLHS=PETSC_FALSE,optChk=PETSC_FALSE;
  //perr=PetscOptionsGetBool(NULL,"-test_lhs",&testLHS,&optChk);CHKERRQ(perr);
  perr=PetscOptionsGetBool(NULL,NULL,"-test_lhs",&testLHS,&optChk);CHKERRQ(perr);

  if (optChk && testLHS) {
    perr=PetscPrintf(PETSC_COMM_WORLD,"WARNING: Using test form of %s\n",
		     __func__);CHKERRQ(perr);
  }
  PetscScalar Hx,Hy,Hz,HyHzdHx,HxHzdHy,HxHydHz;
  Hx = 1.0/(nx*nsubdomains);
  Hy = 1.0/ny;
  Hz = 1.0/nz;
  HyHzdHx = Hy*Hz/Hx;
  HxHzdHy = Hx*Hz/Hy;
  HxHydHz = Hx*Hy/Hz;

  // Get local portion of operator matrix
  perr=MatGetOwnershipRange(A,&rowStart,&rowEnd);CHKERRQ(perr);

  for (ir=rowStart;ir<rowEnd;ir++) {
    /* //Glenn's attempt
    if (ndim == 3) {
      iz = ir/(nx*ny);
      iy = ir%ny;
      ix = (ir/ny)%nx-xshift;
    }
    else {
      iz = 0;
      iy = ir%ny;
      ix = ir/ny-xshift;
    }
    */ 
    izg = ir/(nx*nsubdomains*ny);
    ixg = ir/ny - izg*nx*nsubdomains;
    iyg = ir%ny;
    irl = ir - rowStart;
    iy = irl%ny;
    iz = irl/(nxl*ny);
    ix = irl/ny - iz*nxl + xoffset;

    // Declare adjacent indices. 
    // Guard cells currently take care of BC in X
    // May want to implement for Y and Z.
    ixp1=ix+1;
    ixm1=ix-1;
    if (boundary_type[0] == PERIODIC) {
      if (ixp1 > nx-1 && nsubdomains == 1) ixp1=0;
      if (ixm1 < 0 && nsubdomains == 1) ixm1=nx-1;
    }
    iyp1=iy+1; if (iyp1 > ny-1) iyp1=0;
    iym1=iy-1; if (iym1 < 0) iym1=ny-1;
    izp1=iz+1; if (izp1 > nz-1) izp1=0;
    izm1=iz-1; if (izm1 < 0) izm1=nz-1;

    ixgp1=ixg+1;
    ixgm1=ixg-1;
    if (boundary_type[0] == PERIODIC) {
      if (ixgp1 > nx*nsubdomains-1) ixgp1=0;
      if (ixgm1 < 0) ixgm1=nx*nsubdomains-1;
    }
    iygp1=iyg+1; if (iygp1 > ny-1) iygp1=0;
    iygm1=iyg-1; if (iygm1 < 0) iygm1=ny-1;
    izgp1=izg+1; if (izgp1 > nz-1) izgp1=0;
    izgm1=izg-1; if (izgm1 < 0) izgm1=nz-1;    

    // Define short-hand variables for convenience
    den_xyz = den(INDICIES(ix,iy,iz));
    den_xp1 = den(INDICIES(ixp1,iy,iz));
    den_xm1 = den(INDICIES(ixm1,iy,iz));
    den_yp1 = den(INDICIES(ix,iyp1,iz));
    den_ym1 = den(INDICIES(ix,iym1,iz));
    den_zp1 = den(INDICIES(ix,iy,izp1));
    den_zm1 = den(INDICIES(ix,iy,izm1));

    // Diagonal term
    count = 0;
    index[count] = ir;
    value_diag = -(den_xp1+2*den_xyz+den_xm1)*0.5*dx_sqrInv
      -(den_yp1+2*den_xyz+den_ym1)*0.5*dy_sqrInv;
    if (NDIM==3) value_diag += -(1+kappa*kappa)
    		   *(den_zp1+2*den_xyz+den_zm1)*0.5*dz_sqrInv;
    values[count] = 1.0;
    value_diagInv = 1.0/value_diag;
    if (testLHS) {
      // values[count] = 2.0*(HyHzdHx + HxHzdHy);
      // if (NDIM == 3) values[count] += 2.0*HxHydHz;
      values[count] = 1.0;
    }
    count++;

    // Forward 1 in X
    // index[count] = iyg + ixgp1*ny + izg*ny*nx*nsubdomains;
    idx = iyg + ixgp1*ny + izg*ny*nx*nsubdomains;
    index[count] = (idx < global_length) ? idx : -1;
    values[count] = (den_xyz+den_xp1)*0.5*dx_sqrInv + c1*(den_yp1-den_ym1);
    values[count] *= value_diagInv;
    if (testLHS) {
      // values[count] = -HyHzdHx;
      values[count] = 0.0;
    }
    count++;

    // Backward 1 in X
    // index[count] = iyg + ixgm1*ny + izg*ny*nx*nsubdomains;
    idx = iyg + ixgm1*ny + izg*ny*nx*nsubdomains;
    index[count] = (idx < global_length) ? idx : -1;
    values[count] = (den_xyz+den_xm1)*0.5*dx_sqrInv - c1*(den_yp1-den_ym1);
    values[count] *= value_diagInv;
    if (testLHS) {
      // values[count] = -HyHzdHx;
      values[count] = 0.0;
    }
    count++;

    // Forward 1 in Y
    // index[count] = iygp1 + ixg*ny + izg*ny*nx*nsubdomains;
    idx = iygp1 + ixg*ny + izg*ny*nx*nsubdomains;
    index[count] = (idx < global_length) ? idx : -1;
    values[count] = (den_xyz+den_yp1)*0.5*dy_sqrInv - c1*(den_xp1-den_xm1);
    values[count] *= value_diagInv;
    if (testLHS) {
      // values[count] = -HxHzdHy;
      values[count] = 0.0;
    }
    count++;

    // Backward 1 in Y
    // index[count] = iygm1 + ixg*ny + izg*ny*nx*nsubdomains;
    idx = iygm1 + ixg*ny + izg*ny*nx*nsubdomains;
    index[count] = (idx < global_length) ? idx : -1;
    values[count] = (den_xyz+den_ym1)*0.5*dy_sqrInv + c1*(den_xp1-den_xm1);
    values[count] *= value_diagInv;
    if (testLHS) {
      // values[count] = -HxHzdHy;
      values[count] = 0.0;
    }
    count++;

    if (ndim == 3) {
      // Forward 1 in Z
      // index[count] = iyg + ixg*ny + izgp1*ny*nx*nsubdomains;
      idx = iyg + ixg*ny + izgp1*ny*nx*nsubdomains;
      index[count] = (idx < global_length) ? idx : -1;
      values[count] = (1+kappa*kappa)*(den_xyz+den_zp1)*0.5*dz_sqrInv;
      values[count] *= value_diagInv;
      if (testLHS) {
	// values[count] = -HxHydHz;
      values[count] = 0.0;
      }      
      count++;
      
      // Backward 1 in Z
      // index[count] = iyg + ixg*ny + izgm1*ny*nx*nsubdomains;
      idx = iyg + ixg*ny + izgm1*ny*nx*nsubdomains;
      index[count] = (idx < global_length) ? idx : -1;
      values[count] = (1+kappa*kappa)*(den_xyz+den_zm1)*0.5*dz_sqrInv;
      values[count] *= value_diagInv;
      if (testLHS) {
	// values[count] = -HxHydHz;
      values[count] = 0.0;
      }
      count++;
    }

    // Diagonal: X+1, Y+1
    // index[count] = iygp1 + ixgp1*ny + izg*ny*nx*nsubdomains;
    idx = iygp1 + ixgp1*ny + izg*ny*nx*nsubdomains;
    index[count] = (idx < global_length) ? idx : -1;
    values[count] = c1*(den_yp1-den_xp1);
    values[count] *= value_diagInv;
    if (testLHS) {
      values[count] = 0.0;
    }
    count++;

    // Diagonal: X+1, Y-1
    // index[count] = iygm1 + ixgp1*ny + izg*ny*nx*nsubdomains;
    idx = iygm1 + ixgp1*ny + izg*ny*nx*nsubdomains;
    index[count] = (idx < global_length) ? idx : -1;
    values[count] = c1*(den_xp1-den_ym1);
    values[count] *= value_diagInv;
    if (testLHS) {
      values[count] = 0.0;
    }
    count++;

    // Diagonal: X-1, Y+1
    // index[count] = iygp1 + ixgm1*ny + izg*ny*nx*nsubdomains;
    idx = iygp1 + ixgm1*ny + izg*ny*nx*nsubdomains;
    index[count] = (idx < global_length) ? idx : -1;
    values[count] = c1*(den_xm1-den_yp1);
    values[count] *= value_diagInv;
    if (testLHS) {
      values[count] = 0.0;
    }
    count++;

    // Diagonal: X-1, Y-1
    // index[count] = iygm1 + ixgm1*ny + izg*ny*nx*nsubdomains;
    idx = iygm1 + ixgm1*ny + izg*ny*nx*nsubdomains;
    index[count] = (idx < global_length) ? idx : -1;
    values[count] = c1*(den_ym1-den_xm1);
    values[count] *= value_diagInv;
    if (testLHS) {
      values[count] = 0.0;
    }
    count++;

    // Insert values into current row of the matrix
    perr=MatSetValues(A,1,&ir,count,index,values,INSERT_VALUES);CHKERRQ(perr);

  } // ir

  // Call PETSc assembly routine
  perr=MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(perr);
  perr=MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(perr);

  PetscFunctionReturn(0);
}
#endif // HAVE_PETSC
#endif // USE_QN
