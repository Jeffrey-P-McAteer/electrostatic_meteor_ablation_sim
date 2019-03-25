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

PetscErrorCode getPetscRHS_flux(Vec &b, field &Efield, 
				FArrayND_ranged &den, 
				INDICIES(FArrayND_ranged &Gx,	\
					 FArrayND_ranged &Gy,	\
					 FArrayND_ranged &Gz))
{

  PetscFunctionBegin;

  /* Local function declarations */
  // void periodic_filter(FArrayND_ranged &array,const int nx_guard[]);
  // void pass_guards(ArrayNd_ranged<FTYPE,NDIM> &in, const int nx_guard[]);

  /* DEV: Create debug array */
  void write_local_bin(FArrayND_ranged array, const char *name, int *nx_ghost);
  void write_local_bin(FArrayND array, const char *name);
  FArrayND petscRHS=FArrayND(INDICIES(nx,ny,nz));

  int wxguards[] = {0,0}; // Diagnostics
  if (hybrid_diagnostic_subcycle > 0 && 
      it%nout*hybrid_diagnostic_subcycle == 0) {
    char wlb_name[32];
    sprintf(wlb_name,"qden_RHS-%06d-",it);
    write_local_bin(den,wlb_name,wxguards);
    sprintf(wlb_name,"Gx_RHS-%06d-",it);
    write_local_bin(Gx,wlb_name,wxguards);
    sprintf(wlb_name,"Gy_RHS-%06d-",it);
    write_local_bin(Gy,wlb_name,wxguards);
  }

  // Error tracking
  PetscErrorCode perr;

  // Grid and stencil variables
  PetscInt ir,rowStart,rowEnd,irl;
  int ix,iy,iz;
  int ixp1,iyp1,izp1,ixm1,iym1,izm1;
  // int ixp2,iyp2,izp2,ixm2,iym2,izm2;
  // int nxl = nx/subdomain.np;
  int nxl = nx/petsc_np;
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
  FTYPE ddenx,ddeny,ddenz;
  FTYPE d2denx,d2deny,d2denz;
  FTYPE dGxx,dGyy,dGzz;
  FTYPE Gx_xp1,Gx_yp1,Gx_zp1,Gx_xm1,Gx_ym1,Gx_zm1;
  FTYPE Gy_xp1,Gy_yp1,Gy_zp1,Gy_xm1,Gy_ym1,Gy_zm1;
  FTYPE Gz_xp1,Gz_yp1,Gz_zp1,Gz_xm1,Gz_ym1,Gz_zm1;

  // Physical quantities
  FTYPE Omega_e = Bz*fabs(qd[electron_dist])/md[electron_dist];
  FTYPE kappa = Omega_e/coll_rate[electron_dist];
  FTYPE c2 = thermal_gamma[electron_dist]/fabs(qd[electron_dist]);
  FTYPE c3 = (1+kappa*kappa)*md[electron_dist]*
    coll_rate[electron_dist]/fabs(qd[electron_dist]);
  FTYPE delta_en = 4.e-3;
  FTYPE c4 = 2.*md[electron_dist]/(3.*delta_en);
  FTYPE Te = md[electron_dist]*(Sqr(vxthd[electron_dist])+
				Sqr(vythd[electron_dist]))/2.;
  FTYPE Tn = Sqr(vxthd_neutral[electron_dist])+
    Sqr(vythd_neutral[electron_dist]);
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
    if (boundary_type == PERIODIC) {
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
    den_ppz = den(INDICIES(ixp1,iyp1,iz));
    den_pmz = den(INDICIES(ixp1,iym1,iz));
    den_mpz = den(INDICIES(ixm1,iyp1,iz));
    den_mmz = den(INDICIES(ixm1,iym1,iz));
    Gx_xp1 = Gx(INDICIES(ixp1,iy,iz));
    Gx_xm1 = Gx(INDICIES(ixm1,iy,iz));
    Gx_yp1 = Gx(INDICIES(ix,iyp1,iz));
    Gx_ym1 = Gx(INDICIES(ix,iym1,iz));
    Gx_zp1 = Gx(INDICIES(ix,iy,izp1));
    Gx_zm1 = Gx(INDICIES(ix,iy,izm1));
    Gy_xp1 = Gy(INDICIES(ixp1,iy,iz));
    Gy_xm1 = Gy(INDICIES(ixm1,iy,iz));
    Gy_yp1 = Gy(INDICIES(ix,iyp1,iz));
    Gy_ym1 = Gy(INDICIES(ix,iym1,iz));
    Gy_zp1 = Gy(INDICIES(ix,iy,izp1));
    Gy_zm1 = Gy(INDICIES(ix,iy,izm1));
#if NDIM == 3
    Gz_xp1 = Gz(INDICIES(ixp1,iy,iz));
    Gz_xm1 = Gz(INDICIES(ixm1,iy,iz));
    Gz_yp1 = Gz(INDICIES(ix,iyp1,iz));
    Gz_ym1 = Gz(INDICIES(ix,iym1,iz));
    Gz_zp1 = Gz(INDICIES(ix,iy,izp1));
    Gz_zm1 = Gz(INDICIES(ix,iy,izm1));
#endif

    // Define derivatives
    ddenx = (den_xp1-den_xm1)*0.5*dxInv;
    ddeny = (den_yp1-den_ym1)*0.5*dyInv;
    ddenz = (den_zp1-den_zm1)*0.5*dzInv;
    dGxx = (Gx_xp1-Gx_xm1)*0.5*dxInv;
    dGyy = (Gy_yp1-Gy_ym1)*0.5*dyInv;
#if NDIM == 3
    dGzz = (Gz_zp1-Gz_zm1)*0.5*dzInv;
#endif
    d2denx = (den_xp1-2*den_xyz+den_xm1)*dx_sqrInv;
    d2deny = (den_yp1-2*den_xyz+den_ym1)*dy_sqrInv;
    d2denz = (den_zp1-2*den_xyz+den_zm1)*dz_sqrInv;

    // Build RHS term
#if NDIM == 2
    value_diag = -(den_xp1+2*den_xyz+den_xm1)*0.5*dx_sqrInv
      -(den_yp1+2*den_xyz+den_ym1)*0.5*dy_sqrInv;

    value = 
      ddenx*(Ex0_external - kappa*Ey0_external) +
      ddeny*(kappa*Ex0_external + Ey0_external) +
      c2*Te*(d2denx + d2deny) +
      c3*(dGxx + dGyy);

    value /= value_diag;
#else
    value_diag = -(den_xp1+2*den_xyz+den_xm1)*0.5*dx_sqrInv
      -(den_yp1+2*den_xyz+den_ym1)*0.5*dy_sqrInv
      -(1+kappa*kappa)
      *(den_zp1+2*den_xyz+den_zm1)*0.5*dz_sqrInv;

    value = 
      ddenx*(Ex0_external - kappa*Ey0_external) +
      ddeny*(kappa*Ex0_external + Ey0_external) +
      ddenz*(1+kappa*kappa)*Ez0_external +
      c2*Te*(d2denx + d2deny + (1+kappa*kappa)*d2denz) +
      c3*(dGxx + dGyy + dGzz);

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

    // Insert value into current row of the vector
    perr=VecSetValue(b,ir,value,INSERT_VALUES);CHKERRQ(perr);

  } // ir

  // Write out debugging array
  if (hybrid_diagnostic_subcycle > 0 && 
      it%nout*hybrid_diagnostic_subcycle == 0) {
    petscRHS(INDICIES(ix,iy,iz)) = value;
    char wlb_name[32];
    sprintf(wlb_name,"petscRHS-%06d-",it);
    write_local_bin(petscRHS,wlb_name);
  }


  // Call PETSc assembly routine
  perr=VecAssemblyBegin(b);CHKERRQ(perr);
  perr=VecAssemblyEnd(b);CHKERRQ(perr);  

  PetscFunctionReturn(0);
}
#endif // HAVE_PETSC
#endif // USE_QN
