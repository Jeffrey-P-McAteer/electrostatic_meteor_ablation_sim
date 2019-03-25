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

PetscErrorCode getPetscRHS(Vec &b, field &Efield, 
			   FArrayND_ranged &den, FArrayND_ranged &divGamma)
// PetscErrorCode getPetscRHS(Vec &b, field &Efield, 
// 			   FArrayND_ranged &den, 
// 			   INDICIES(FArrayND_ranged &Gx,FArrayND_ranged &Gy,FArrayND_ranged &Gz))
{

  PetscFunctionBegin;

  /* DEV: Create debug array */
  void write_local_bin(FArrayND_ranged array, const char *name, int *nx_ghost);
  void write_local_bin(FArrayND array, const char *name);
  FArrayND petscRHS=FArrayND(INDICIES(nx,ny,nz));

  // Error tracking
  PetscErrorCode perr;

  // Extract phi for convenience and speed
  FArrayND_ranged phi = Efield.phi_rho;

  // Grid and stencil variables
  PetscInt ir,rowStart,rowEnd,irl;
  int ix,iy,iz,ixp1,iyp1,izp1,ixm1,iym1,izm1;
  int nxl = nx/subdomain.np;
  int xoffset = nxl*subdomain.rank;
  int xshift=nx*subdomain.id_number;
  PetscScalar value = 0.0;
  FTYPE value_diag = 1.0;
  FTYPE dx_sqr = dx*dx;
  FTYPE dy_sqr = dy*dy;
  FTYPE dz_sqr = dz*dz;
  FTYPE dxInv = 1.0/dx,dyInv = 1.0/dy,dzInv = 1.0/dz; 
  FTYPE dx_sqrInv = 1.0/dx_sqr,dy_sqrInv = 1.0/dy_sqr,dz_sqrInv = 1.0/dz_sqr; 

  // Density-related quantities
  FTYPE den_xyz,den_xp1,den_yp1,den_zp1,den_xm1,den_ym1,den_zm1;
  FTYPE ddenx,ddeny,ddenz,d2denx,d2deny,d2denz;

  // Potential-related quantities
  FTYPE dphix,dphiy,dphiz;

  // Other physical quantities
  FTYPE nu_e = coll_rate[electron_dist],m_e = md[electron_dist],q_e = qd[electron_dist];
  FTYPE Omega_e = Bz*fabs(q_e)/m_e;
  FTYPE kappa_e = Omega_e/nu_e;
  FTYPE c2 = thermal_gamma[electron_dist]/fabs(q_e);
  FTYPE c3 = (1+kappa_e*kappa_e)*m_e*nu_e/fabs(q_e);
  FTYPE delta_en = 4.e-3;
  FTYPE c4 = 2.*m_e/(3.*delta_en);
  FTYPE kTe = m_e*(Sqr(vxthd[electron_dist])+Sqr(vythd[electron_dist]))/2.;
  FTYPE kTn = Sqr(vxthd_neutral[electron_dist])+Sqr(vythd_neutral[electron_dist]);
  if (NDIM == 3) kTn += Sqr(vzthd_neutral[electron_dist]);
  kTn *= m_neutral/2.;
  // FTYPE veSqr = 1./(Bz*Bz);
  FTYPE vex = 0.0, vey = 0.0, veSqr = 0.0;
  
  // Check if user has set "-test_rhs" in the PETSc options file
  PetscBool testRHS=PETSC_FALSE,optChk=PETSC_FALSE;

  //perr=PetscOptionsGetBool(NULL,"-test_rhs",&testRHS,&optChk);CHKERRQ(perr);
  perr=PetscOptionsGetBool(NULL,NULL,"-test_rhs",&testRHS,&optChk);CHKERRQ(perr);

  if (optChk && testRHS) {
    perr=PetscPrintf(PETSC_COMM_WORLD,"WARNING: Using test form of %s\n",
		     __func__);CHKERRQ(perr);
  }

  // Get local portion of forcing vector
  perr=VecGetOwnershipRange(b,&rowStart,&rowEnd);CHKERRQ(perr);
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

    // Define short-hand variables for convenience
    den_xyz = den(INDICIES(ix,iy,iz));
    den_xp1 = den(INDICIES(ixp1,iy,iz));
    den_xm1 = den(INDICIES(ixm1,iy,iz));
    den_yp1 = den(INDICIES(ix,iyp1,iz));
    den_ym1 = den(INDICIES(ix,iym1,iz));
    den_zp1 = den(INDICIES(ix,iy,izp1));
    den_zm1 = den(INDICIES(ix,iy,izm1));

    // Define derivatives
    ddenx = (den_xp1-den_xm1)*0.5*dxInv;
    ddeny = (den_yp1-den_ym1)*0.5*dyInv;
    ddenz = (den_zp1-den_zm1)*0.5*dzInv;
    d2denx = (den_xp1-2*den_xyz+den_xm1)*dx_sqrInv;
    d2deny = (den_yp1-2*den_xyz+den_ym1)*dy_sqrInv;
    d2denz = (den_zp1-2*den_xyz+den_zm1)*dz_sqrInv;
    dphix = (phi(INDICIES(ixp1,iy,iz))-phi(INDICIES(ixm1,iy,iz)))*0.5*dxInv;
    dphiy = (phi(INDICIES(ix,iyp1,iz))-phi(INDICIES(ix,iym1,iz)))*0.5*dyInv;
    dphix = (phi(INDICIES(ix,iy,izp1))-phi(INDICIES(ix,iy,izm1)))*0.5*dzInv;
    
    // Need diagonal value of LHS operator
    value_diag = -(den_xp1+2*den_xyz+den_xm1)*0.5*dx_sqrInv
      -(den_yp1+2*den_xyz+den_ym1)*0.5*dy_sqrInv;
    if (ndim == 3) value_diag += -(1+kappa_e*kappa_e)
    		   *(den_zp1+2*den_xyz+den_zm1)*0.5*dz_sqrInv;
    // value_diag=1.0;

    // It may work nicely to set d2denz = 0.0 if NDIM == 2 so that 
    // the 3-D expression reduces to the 2-D expression automatically
    // (assuming Ez0_external = 0 in 2-D).
    if (ndim == 2) {
      // veSqr *= Sqr(Ex0_external+Efield.Ex(INDICIES(ix,iy,iz))
      // 		   +Ey0_external+Efield.Ey(INDICIES(ix,iy,iz)));
      // // kTe = kTn + c4*veSqr;

      // c2 *= kTn + c4*veSqr;
      // value = c2*(d2denx+d2deny) + c3*divGamma(INDICIES(ix,iy,iz))
      // 	+ddenx*(Ex0_external - kappa_e*Ey0_external)
      // 	+ddeny*(kappa_e*Ex0_external + Ey0_external);

      /* DEV
      vex = -(nu_e*nu_e + Omega_e*Omega_e)
	*(fabs(q_e)/m_e/nu_e*(Ex0_external-dphix)+kTe/nu_e/m_e*ddenx/den_xyz +
	  kappa_e*kappa_e*((dphiy-Ey0_external)/Bz+kTe/q_e/Bz*ddeny/den_xyz));
      vey = -(nu_e*nu_e + Omega_e*Omega_e)
	*(fabs(q_e)/m_e/nu_e*(Ey0_external-dphiy)+kTe/nu_e/m_e*ddeny/den_xyz +
	  kappa_e*kappa_e*((Ex0_external-dphix)/Bz-kTe/q_e/Bz*ddenx/den_xyz));
      veSqr = vex*vex + vey*vey;
      c2 *= kTn + c4*veSqr;
      */
      c2 *= kTe;
      value = c2*(d2denx+d2deny) + c3*divGamma(INDICIES(ix,iy,iz))
      	+ddenx*(Ex0_external - kappa_e*Ey0_external)
      	+ddeny*(kappa_e*Ex0_external + Ey0_external);

    }
    if (ndim == 3) {
      // veSqr *= Sqr(Ex0_external+Efield.Ex(INDICIES(ix,iy,iz))
      // 		   +Ey0_external+Efield.Ey(INDICIES(ix,iy,iz))
      // 		   +Ez0_external+Efield.Ez(INDICIES(ix,iy,iz)));
      // kTe = kTn + c4*veSqr;

      value = c2*kTe*(d2denx+d2deny+(1+kappa_e*kappa_e)*d2denz) 
      	+c3*divGamma(INDICIES(ix,iy,iz))
      	+ddenx*(Ex0_external - kappa_e*Ey0_external)
      	+ddeny*(kappa_e*Ex0_external + Ey0_external)
      	+ddenz*(1+kappa_e*kappa_e)*Ez0_external;
    }

    // Use test (e.g. known analytical or simplified physical) form of RHS
    if (testRHS) {
      value = 0.0;
      value_diag = 1.0;
      PetscScalar Hx,Hy,Hz;
      Hx = 1.0/(nx*nsubdomains);
      Hy = 1.0/ny;
      Hz = 1.0/nz;
      if (ndim == 2) {
	// value = 8*PETSC_PI*PETSC_PI
	//   *PetscCosScalar(2*PETSC_PI*(((PetscReal)(ix+xshift)+0.5)*Hx))
	//   *PetscCosScalar(2*PETSC_PI*(((PetscReal)iy+0.5)*Hy))
	//   *Hx*Hy;

	// value = iy + (ix-xoffset)*ny + iz*nxl*ny + rowStart;

	value = ir;
      } else if (ndim == 3) {
	value = 12*PETSC_PI*PETSC_PI
	  *PetscCosScalar(2*PETSC_PI*(((PetscReal)(ix+xshift)+0.5)*Hx))
	  *PetscCosScalar(2*PETSC_PI*(((PetscReal)iy+0.5)*Hy))
	  *PetscCosScalar(2*PETSC_PI*(((PetscReal)iz+0.5)*Hz))
	  *Hx*Hy*Hz;
	// value = iy + (ix-xoffset)*ny + iz*nxl*ny + rowStart;
      }
    }

    value /= value_diag;

    petscRHS(INDICIES(ix,iy,iz)) = value;

    // Insert value into current row of the vector
    perr=VecSetValue(b,ir,value,INSERT_VALUES);CHKERRQ(perr);

  } // ir

  // write_local_bin(petscRHS,"petscRHS");

  // Call PETSc assembly routine
  perr=VecAssemblyBegin(b);CHKERRQ(perr);
  perr=VecAssemblyEnd(b);CHKERRQ(perr);  

  PetscFunctionReturn(0);
}
#endif // HAVE_PETSC
#endif // USE_QN
