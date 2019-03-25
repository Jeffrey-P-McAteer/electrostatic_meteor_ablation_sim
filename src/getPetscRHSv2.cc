/* Build the RHS forcing vector (b) of the quasi-neutral, inertialess
   potential equation */

/* NOTES
   --Started 07May15 (may)
   --Including the line "-test_rhs" (without quotes) in the PETSc options
     file causes a user-defined test function to be used. Intended to be 
     used with the "-test_lhs" operator to produce analytic solutions. 
     14May15 (may)
*/

#if USE_QN
#include "eppic.h"
#include "efield_quasineut.h"

#if HAVE_PETSC
#include <petscksp.h>

PetscErrorCode getPetscRHSv2(Vec &b, field &Efield, 
			     FArrayND_ranged &den, 
			     INDICIES(FArrayND_ranged &Gx, \
				      FArrayND_ranged &Gy, \
				      FArrayND_ranged &Gz))
{

  PetscFunctionBegin;

  /* Local function declarations */
  void periodic_filter(FArrayND_ranged &array,const int nx_guard[]);
  void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);

  /* DEV: Create debug array */
  void write_local_bin(FArrayND_ranged array, const char *name, int *nx_ghost);
  void write_local_bin(FArrayND array, const char *name);
  FArrayND petscRHS=FArrayND(INDICIES(nx,ny,nz));

  int wxguards[] = {0,0}; // Diagnostics
  // write_local_bin(den,"qden_RHS",wxguards);
  // write_local_bin(Gx,"Gx_RHS",wxguards);
  // write_local_bin(Gy,"Gy_RHS",wxguards);

  // Error tracking
  PetscErrorCode perr;

  // Grid and stencil variables
  PetscInt ir,rowStart,rowEnd,irl;
  int ix,iy,iz;
  int ixp1,iyp1,izp1,ixm1,iym1,izm1;
  // int ixp2,iyp2,izp2,ixm2,iym2,izm2;
  int nxl = nx/subdomain.np;
  int xoffset = nxl*subdomain.rank;
  int xshift=nx*subdomain.id_number;
  PetscScalar value = 0.0;
  PetscScalar Fxp1,Fyp1,Fxm1,Fym1,Fzp1,Fzm1;
  FTYPE value_diag = 1.0;
  FTYPE dx_sqr = dx*dx;
  FTYPE dy_sqr = dy*dy;
  FTYPE dz_sqr = dz*dz;
  FTYPE dxInv = 1.0/dx,dyInv = 1.0/dy,dzInv = 1.0/dz; 
  FTYPE dx_sqrInv = 1.0/dx_sqr,dy_sqrInv = 1.0/dy_sqr,dz_sqrInv = 1.0/dz_sqr; 

  // Density-related quantities
  FTYPE den_xyz;
  FTYPE den_xp1,den_yp1,den_zp1,den_xm1,den_ym1,den_zm1;
  // FTYPE den_xp2,den_yp2,den_zp2,den_xm2,den_ym2,den_zm2;
  FTYPE den_ppz,den_pmz,den_mpz,den_mmz;
  // FTYPE ddenx,ddeny,ddenz,d2denx,d2deny,d2denz;

  // Physical quantities
  FTYPE Omega_e = Bz*fabs(qd[electron_dist])/md[electron_dist];
  FTYPE kappa = Omega_e/coll_rate[electron_dist];
  FTYPE c2 = thermal_gamma[electron_dist]/fabs(qd[electron_dist]);
  FTYPE c3 = (1+kappa*kappa)*md[electron_dist]*coll_rate[electron_dist]/fabs(qd[electron_dist]);
  FTYPE delta_en = 4.e-3;
  FTYPE c4 = 2.*md[electron_dist]/(3.*delta_en);
  FTYPE Te = md[electron_dist]*(Sqr(vxthd[electron_dist])+Sqr(vythd[electron_dist]))/2.;
  FTYPE Tn = Sqr(vxthd_neutral[electron_dist])+Sqr(vythd_neutral[electron_dist]);
  if (NDIM == 3) Tn += Sqr(vzthd_neutral[electron_dist]);
  Tn *= m_neutral/2.;
  FTYPE veSqr = 1./(Bz*Bz);

  // Check if user has set "-test_rhs" in the PETSc options file
  PetscBool testRHS=PETSC_FALSE,optChk=PETSC_FALSE;
  perr=PetscOptionsGetBool(NULL,NULL,"-test_rhs",&testRHS,&optChk);CHKERRQ(perr);
  if (optChk && testRHS) {
    perr=PetscPrintf(PETSC_COMM_WORLD,"WARNING: Using test form of %s\n",
		     __func__);CHKERRQ(perr);
  }

  // Build grad[n]
  // const int nx_guard[] = {1,1};
  int arr_size[] = {INDICIES(nx+qnx_guard_size[0]+qnx_guard_size[1],ny,nz)};
  int arr_start[] = {INDICIES(0-qnx_guard_size[0],0,0)};
  FArrayND_ranged ddenx,ddeny;
  ddenx = FArrayND_ranged(arr_start,arr_size);
  ddeny = FArrayND_ranged(arr_start,arr_size);
#if NDIM == 3
  FArrayND_ranged ddenz;
  ddenz = FArrayND_ranged(arr_start,arr_size);
#endif
  for (ix=0;ix<nx;ix++) {
    for (iy=0;iy<ny;iy++) {
      for (iz=0;iz<nz;iz++) {
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

	ddenx(INDICIES(ix,iy,iz)) = 
	  0.5*dxInv*(den(INDICIES(ixp1,iy,iz))-den(INDICIES(ixm1,iy,iz)));
	ddeny(INDICIES(ix,iy,iz)) = 
	  0.5*dyInv*(den(INDICIES(ix,iyp1,iz))-den(INDICIES(ix,iym1,iz)));
// 	if (boundary_type[1] == PERIODIC) { // Calls pass_guards internally
// 	  periodic_filter(ddenx,qnx_guard_size);
// 	  periodic_filter(ddeny,qnx_guard_size);
// 	} else {
// 	  pass_guards(ddenx,qnx_guard_size);
// 	  pass_guards(ddeny,qnx_guard_size);
// 	}
// #if NDIM == 3
// 	ddenz(INDICIES(ix,iy,iz)) = 
// 	  0.5*dzInv*(den(INDICIES(ix,iy,izp1))-den(INDICIES(ix,iy,izm1)));
// 	if (boundary_type[2] == PERIODIC) // Calls pass_guards internally
// 	  periodic_filter(ddenz,qnx_guard_size);
// 	else
// 	  pass_guards(ddenz,qnx_guard_size);
// #endif
	pass_guards(ddenx,qnx_guard_size);
	pass_guards(ddeny,qnx_guard_size);
#if NDIM == 3
	pass_guards(ddenz,qnx_guard_size);
#endif
      }
    }
  }
  ix = 0; iy = 0; iz = 0;

  // Get local portion of forcing vector
  perr=VecGetOwnershipRange(b,&rowStart,&rowEnd);CHKERRQ(perr);
  for (ir=rowStart;ir<rowEnd;ir++) {
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
    /* This may require redefining # of guard cells.
       That's straightforward but possibly time consuming.
    ixm1 = ix-1; ixm2 = ix-2;
    ixp1 = ix+1; ixp2 = ix+2;
    iym1 = iy-1; iym2 = iy-2;
    iyp1 = iy+1; iyp2 = iy+2;
    izm1 = iz-1; izm2 = iz-2;
    izp1 = iz+1; izp2 = iz+2;
    if (boundary_type[0] == PERIODIC) {
    }
    */

    // Define short-hand variables for convenience
    den_xyz = den(INDICIES(ix,iy,iz));
    den_xp1 = den(INDICIES(ixp1,iy,iz));
    den_xm1 = den(INDICIES(ixm1,iy,iz));
    den_yp1 = den(INDICIES(ix,iyp1,iz));
    den_ym1 = den(INDICIES(ix,iym1,iz));
    den_zp1 = den(INDICIES(ix,iy,izp1));
    den_zm1 = den(INDICIES(ix,iy,izm1));
    den_ppz = den(INDICIES(ixp1,iyp1,iz));
    den_pmz = den(INDICIES(ixp1,iym1,iz));
    den_mpz = den(INDICIES(ixm1,iyp1,iz));
    den_mmz = den(INDICIES(ixm1,iym1,iz));

    // // Define derivatives
    // ddenx = (den_xp1-den_xm1)*0.5*dxInv;
    // ddeny = (den_yp1-den_ym1)*0.5*dyInv;
    // ddenz = (den_zp1-den_zm1)*0.5*dzInv;
    // d2denx = (den_xp1-2*den_xyz+den_xm1)*dx_sqrInv;
    // d2deny = (den_yp1-2*den_xyz+den_ym1)*dy_sqrInv;
    // d2denz = (den_zp1-2*den_xyz+den_zm1)*dz_sqrInv;

    // Need diagonal value of LHS operator
    // value_diag = -(den_xp1+2*den_xyz+den_xm1)*0.5*dx_sqrInv
    //   -(den_yp1+2*den_xyz+den_ym1)*0.5*dy_sqrInv;
    // if (ndim == 3) value_diag += -(1+kappa*kappa)
    // 		   *(den_zp1+2*den_xyz+den_zm1)*0.5*dz_sqrInv;
    // value_diag=1.0;

    // It may work nicely to set d2denz = 0.0 if NDIM == 2 so that 
    // the 3-D expression reduces to the 2-D expression automatically
    // (assuming Ez0_external = 0 in 2-D).
    // if (ndim == 2) {
    //   // veSqr *= Sqr(Ex0_external+Efield.Ex(INDICIES(ix,iy,iz))
    //   // 		   +Ey0_external+Efield.Ey(INDICIES(ix,iy,iz)));
    //   // Te = Tn + c4*veSqr;
    //   value = c2*Te*(d2denx+d2deny) + c3*divGamma(INDICIES(ix,iy,iz))
    //   	+ddenx*(Ex0_external - kappa*Ey0_external)
    //   	+ddeny*(kappa*Ex0_external + Ey0_external);
    //   value /= value_diag;
    // }
    // if (ndim == 3) {
    //   // veSqr *= Sqr(Ex0_external+Efield.Ex(INDICIES(ix,iy,iz))
    //   // 		   +Ey0_external+Efield.Ey(INDICIES(ix,iy,iz))
    //   // 		   +Ez0_external+Efield.Ez(INDICIES(ix,iy,iz)));
    //   // Te = Tn + c4*veSqr;
    //   value = c2*Te*(d2denx+d2deny+(1+kappa*kappa)*d2denz) 
    //   	+c3*divGamma(INDICIES(ix,iy,iz))
    //   	+ddenx*(Ex0_external - kappa*Ey0_external)
    //   	+ddeny*(kappa*Ex0_external + Ey0_external)
    //   	+ddenz*(1+kappa*kappa)*Ez0_external;
    //   value /= value_diag;
    // }

#if NDIM == 2
      // value_diag = -(dx_sqrInv + dy_sqrInv)*den_xyz
      // 	-0.50*(dx_sqrInv + dy_sqrInv)*(den_xp1 + den_yp1 + den_xm1 + den_ym1)
      // 	-0.25*(den_ppz + den_pmz + den_mpz + den_mmz);
    value_diag = -(den_xp1+2*den_xyz+den_xm1)*0.5*dx_sqrInv
      -(den_yp1+2*den_xyz+den_ym1)*0.5*dy_sqrInv;
    Fxm1 = den(INDICIES(ixm1,iy,iz))*(Ex0_external - kappa*Ey0_external)+
      c2*Te*(ddenx(INDICIES(ixm1,iy,iz))-kappa*ddeny(INDICIES(ixm1,iy,iz)))+
      c3*Gx(INDICIES(ixm1,iy,iz));
    Fxp1 = den(INDICIES(ixp1,iy,iz))*(Ex0_external - kappa*Ey0_external)+
      c2*Te*(ddenx(INDICIES(ixp1,iy,iz))-kappa*ddeny(INDICIES(ixp1,iy,iz)))+
      c3*Gx(INDICIES(ixp1,iy,iz));
    Fym1 = den(INDICIES(ix,iym1,iz))*(kappa*Ex0_external + Ey0_external)+
      c2*Te*(kappa*ddenx(INDICIES(ix,iym1,iz))+ddeny(INDICIES(ix,iym1,iz)))+
      c3*Gy(INDICIES(ix,iym1,iz));
    Fyp1 = den(INDICIES(ix,iyp1,iz))*(kappa*Ex0_external + Ey0_external)+
      c2*Te*(kappa*ddenx(INDICIES(ix,iyp1,iz))+ddeny(INDICIES(ix,iyp1,iz)))+
      c3*Gy(INDICIES(ix,iyp1,iz));
    value = 0.5*dxInv*(Fxp1-Fxm1) + 0.5*dyInv*(Fyp1-Fym1);
    value /= value_diag;
#else
    // value_diag = -(dx_sqrInv + dy_sqrInv + (1+kappa*kappa)*dz_sqrInv)*den_xyz
    //   -0.50*((dx_sqrInv + dy_sqrInv)*(den_xp1 + den_yp1 + den_xm1 + den_ym1) +
    // 	     (1+kappa*kappa)*dz_sqrInv*(den_zp1 + den_zm1))
    //   -0.25*(den_ppz + den_pmz + den_mpz + den_mmz);
    value_diag = -(den_xp1+2*den_xyz+den_xm1)*0.5*dx_sqrInv
      -(den_yp1+2*den_xyz+den_ym1)*0.5*dy_sqrInv
      -(1+kappa*kappa)
      *(den_zp1+2*den_xyz+den_zm1)*0.5*dz_sqrInv;
    Fxm1 = den(INDICIES(ixm1,iy,iz))*(Ex0_external - kappa*Ey0_external)+
      c2*Te*(ddenx(INDICIES(ixm1,iy,iz))-kappa*ddeny(INDICIES(ixm1,iy,iz)))+
      c3*Gx(INDICIES(ixm1,iy,iz));
    Fxp1 = den(INDICIES(ixp1,iy,iz))*(Ex0_external - kappa*Ey0_external)+
      c2*Te*(ddenx(INDICIES(ixp1,iy,iz))-kappa*ddeny(INDICIES(ixp1,iy,iz)))+
      c3*Gx(INDICIES(ixp1,iy,iz));
    Fym1 = den(INDICIES(ix,iym1,iz))*(kappa*Ex0_external + Ey0_external)+
      c2*Te*(kappa*ddenx(INDICIES(ix,iym1,iz))+ddeny(INDICIES(ix,iym1,iz)))+
      c3*Gy(INDICIES(ix,iym1,iz));
    Fyp1 = den(INDICIES(ix,iyp1,iz))*(kappa*Ex0_external + Ey0_external)+
      c2*Te*(kappa*ddenx(INDICIES(ix,iyp1,iz))+ddeny(INDICIES(ix,iyp1,iz)))+
      c3*Gy(INDICIES(ix,iyp1,iz));
    Fzm1 = den(INDICIES(ix,iy,izm1))*(1+kappa*kappa)*Ez0_external+
      c2*Te*(1+kappa*kappa)*ddenz(INDICIES(ix,iy,izm1))+
      c3*Gz(INDICIES(ix,iy,izm1));
    Fzp1 = den(INDICIES(ix,iy,izp1))*(1+kappa*kappa)*Ez0_external+
      c2*Te*(1+kappa*kappa)*ddenz(INDICIES(ix,iy,izp1))+
      c3*Gz(INDICIES(ix,iy,izp1));
    value = 0.5*dxInv*(Fxp1-Fxm1) + 0.5*dyInv*(Fyp1-Fym1) + 0.5*dzInv*(Fzp1-Fzm1);
    value /= value_diag;
#endif

    // Use test (e.g. known analytical or simplified physical) form of RHS
    if (testRHS) {
      value = 0.0;
      value_diag = 1.0;
      PetscScalar Hx,Hy,Hz;
      Hx = 1.0/(nx*nsubdomains);
      Hy = 1.0/ny;
      Hz = 1.0/nz;
      if (ndim == 2) {
	value = 8*PETSC_PI*PETSC_PI
	  *PetscCosScalar(2*PETSC_PI*(((PetscReal)(ix+xshift)+0.5)*Hx))
	  *PetscCosScalar(2*PETSC_PI*(((PetscReal)iy+0.5)*Hy))
	  *Hx*Hy;
	// value = iy + (ix-xoffset)*ny + iz*nxl*ny + rowStart;
      } else if (ndim == 3) {
	value = 12*PETSC_PI*PETSC_PI
	  *PetscCosScalar(2*PETSC_PI*(((PetscReal)(ix+xshift)+0.5)*Hx))
	  *PetscCosScalar(2*PETSC_PI*(((PetscReal)iy+0.5)*Hy))
	  *PetscCosScalar(2*PETSC_PI*(((PetscReal)iz+0.5)*Hz))
	  *Hx*Hy*Hz;
	// value = iy + (ix-xoffset)*ny + iz*nxl*ny + rowStart;
      }
    }

    // Write current value to debugging array
    petscRHS(INDICIES(ix,iy,iz)) = value;

    // Insert value into current row of the vector
    perr=VecSetValue(b,ir,value,INSERT_VALUES);CHKERRQ(perr);

  } // ir

  // // Write out debugging array
  // write_local_bin(petscRHS,"petscRHS");

  // Call PETSc assembly routine
  perr=VecAssemblyBegin(b);CHKERRQ(perr);
  perr=VecAssemblyEnd(b);CHKERRQ(perr);  

  PetscFunctionReturn(0);
}
#endif // HAVE_PETSC
#endif // USE_QN
