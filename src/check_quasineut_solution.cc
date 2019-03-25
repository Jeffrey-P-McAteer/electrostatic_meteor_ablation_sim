#if USE_QN
#include "eppic.h"
#include "efield_quasineut.h"

int check_quasineut_solution(field &Efield,
			     FArrayND_ranged &den,
			     INDICIES(FArrayND_ranged &Gx,
				      FArrayND_ranged &Gy,
				      FArrayND_ranged &Gz),
			     FArrayND_ranged &phi,
			     int check_solution)
{

  // Grid and stencil variables
  int ix,iy,iz;
  int ixp1,iyp1,izp1,ixm1,iym1,izm1;
  int ix_max=0,iy_max=0,iz_max=0;
  // int nxl = nx/subdomain.np;
  int nxl = nx/petsc_np;
  int xoffset = nxl*subdomain.rank;
  int xshift=nx*subdomain.id_number;
  FTYPE dx_sqr = dx*dx;
  FTYPE dy_sqr = dy*dy;
  FTYPE dz_sqr = dz*dz;
  FTYPE dxInv = 1.0/dx,dyInv = 1.0/dy,dzInv = 1.0/dz; 
  FTYPE dx_sqrInv = 1.0/dx_sqr,dy_sqrInv = 1.0/dy_sqr,dz_sqrInv = 1.0/dz_sqr; 

  // Data-related quantities
  FTYPE den_xyz;
  FTYPE den_xp1,den_yp1,den_zp1,den_xm1,den_ym1,den_zm1;
  FTYPE phi_xyz;
  FTYPE phi_xp1,phi_yp1,phi_zp1,phi_xm1,phi_ym1,phi_zm1;
  FTYPE den_ppz,den_pmz,den_mpz,den_mmz;
  FTYPE phi_ppz,phi_pmz,phi_mpz,phi_mmz;
  FTYPE ddenx,ddeny,ddenz;
  FTYPE d2denx,d2deny,d2denz;
  FTYPE dGxx,dGyy,dGzz;
  FTYPE Gx_xp1,Gx_yp1,Gx_zp1,Gx_xm1,Gx_ym1,Gx_zm1;
  FTYPE Gy_xp1,Gy_yp1,Gy_zp1,Gy_xm1,Gy_ym1,Gy_zm1;
  FTYPE Gz_xp1,Gz_yp1,Gz_zp1,Gz_xm1,Gz_ym1,Gz_zm1;
  FTYPE rhs,lhs;
  FTYPE max_diff=0.0;

  // Physical quantities
  FTYPE Omega_e = Bz*fabs(qd[electron_dist])/md[electron_dist];
  FTYPE kappa = Omega_e/coll_rate[electron_dist];
  FTYPE c1 = kappa/(8*dx*dy);
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
  
  for (long ir=0;ir<nx*ny*nz;ir++) {
    iy = ir%ny;
    iz = ir/(nxl*ny);
    ix = ir/ny - iz*nxl + xoffset;

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
    phi_xyz = phi(INDICIES(ix,iy,iz));
    phi_xp1 = phi(INDICIES(ixp1,iy,iz));
    phi_xm1 = phi(INDICIES(ixm1,iy,iz));
    phi_yp1 = phi(INDICIES(ix,iyp1,iz));
    phi_ym1 = phi(INDICIES(ix,iym1,iz));
    phi_zp1 = phi(INDICIES(ix,iy,izp1));
    phi_zm1 = phi(INDICIES(ix,iy,izm1));
    phi_ppz = phi(INDICIES(ixp1,iyp1,iz));
    phi_pmz = phi(INDICIES(ixp1,iym1,iz));
    phi_mpz = phi(INDICIES(ixm1,iyp1,iz));
    phi_mmz = phi(INDICIES(ixm1,iym1,iz));

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

    // Calculate value of RHS vector (b)
    rhs = 
      ddenx*(Ex0_external - kappa*Ey0_external) +
      ddeny*(kappa*Ex0_external + Ey0_external) +
      c2*Te*(d2denx + d2deny) +
      c3*(dGxx + dGyy);

    // Calculate value of LHS vector (Ax)
    lhs =
      phi_xyz*(-(den_xp1+2*den_xyz+den_xm1)*0.5*dx_sqrInv
	       -(den_yp1+2*den_xyz+den_ym1)*0.5*dy_sqrInv) +
      phi_xp1*((den_xyz+den_xp1)*0.5*dx_sqrInv + c1*(den_yp1-den_ym1)) +
      phi_xm1*((den_xyz+den_xm1)*0.5*dx_sqrInv - c1*(den_yp1-den_ym1)) +
      phi_yp1*((den_xyz+den_yp1)*0.5*dy_sqrInv - c1*(den_xp1-den_xm1)) +
      phi_ym1*((den_xyz+den_ym1)*0.5*dy_sqrInv + c1*(den_xp1-den_xm1)) +
      phi_ppz*(c1*(den_yp1-den_xp1)) +
      phi_pmz*(c1*(den_xp1-den_ym1)) +
      phi_mpz*(c1*(den_xm1-den_yp1)) +
      phi_mmz*(c1*(den_ym1-den_xm1));
    if (ndim == 3) {
      lhs += 
	phi_zp1*((1+kappa*kappa)*(den_xyz+den_zp1)*0.5*dz_sqrInv) +
	phi_zm1*((1+kappa*kappa)*(den_xyz+den_zm1)*0.5*dz_sqrInv);
    }

  //   if (fabs(rhs-lhs)/fabs(lhs) > max_diff) {
  //     ix_max = ix;
  //     iy_max = iy;
  //     iz_max = iz;
  //   }
  //   max_diff = max(max_diff,fabs(rhs-lhs)/fabs(lhs));
  //   if (check_solution == 3)
  //     printf("[%d] %s ir = %d \t max(|Ax-b|/|Ax|) = %e\n",
  // 	     mpi_rank,__func__,ir,max_diff);
  //   if (check_solution == 4)
  //     printf("[%d] %s\n" \
  // 	     "     (ix,iy,iz) = (%d,%d,%d)\n" \
  // 	     "     Ax = %f\n" \
  // 	     "     b = %f\n" \
  // 	     "     |Ax-b|/|Ax| = %e\n",
  // 	     mpi_rank,__func__,ix,iy,iz,lhs,rhs,fabs(lhs-rhs)/fabs(lhs));
  // }
  // if (check_solution == 2)
  //   printf("[%d] %s max(|Ax-b|/|Ax|) = %e at (%d,%d,%d)\n",
  // 	   mpi_rank,__func__,max_diff,ix_max,iy_max,iz_max);
  // FTYPE max_diff_all = 0.0;
  // if (check_solution == 1) {
  //   MPI_Allreduce(&max_diff,&max_diff_all,1,
  // 		  MPI_FTYPE,MPI_MAX,MPI_COMM_WORLD);
  //   if (mpi_rank == 0) 
  //     printf("%s global max(|Ax-b|/|Ax|) = %e\n",
  // 	     __func__,max_diff_all);
  // }
    if (fabs(rhs-lhs)/fabs(rhs) > max_diff) {
      ix_max = ix;
      iy_max = iy;
      iz_max = iz;
    }
    max_diff = max(max_diff,fabs(rhs-lhs)/fabs(rhs));
    if (check_solution == 3)
      printf("[%d] %s ir = %ld \t max(|Ax-b|/|b|) = %e\n",
	     mpi_rank,__func__,ir,max_diff);
    if (check_solution == 4)
      printf("[%d] %s\n" \
	     "     (ix,iy,iz) = (%d,%d,%d)\n" \
	     "     Ax = %f\n" \
	     "     b = %f\n" \
	     "     |Ax-b|/|b| = %e\n",
	     mpi_rank,__func__,ix,iy,iz,lhs,rhs,fabs(lhs-rhs)/fabs(rhs));
  }
  if (check_solution == 2)
    printf("[%d] %s max(|Ax-b|/|b|) = %e at (%d,%d,%d)\n",
	   mpi_rank,__func__,max_diff,ix_max,iy_max,iz_max);
  FTYPE max_diff_all = 0.0;
  if (check_solution == 1) {
    // MPI_Allreduce(&max_diff,&max_diff_all,1,
    // 		  MPI_FTYPE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&max_diff,&max_diff_all,1,
		  MPI_FTYPE,MPI_MAX,PETSC_COMM_WORLD);
    if (mpi_rank == 0) 
      printf("%s global max(|Ax-b|/|b|) = %e\n",
	     __func__,max_diff_all);
  }

  return(0);
}

#endif
