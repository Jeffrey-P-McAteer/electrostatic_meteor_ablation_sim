

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "ArrayNd.h"
#include "eppic.h"
#include "eppic-mpi.h"

int efield_quasineut_solver_test(fluid *fspecie,FArrayND_ranged &phi,FArrayND_ranged &den,FArrayND_ranged &Jx,FArrayND_ranged &Jy, FTYPE err_tol) {


  int failed = 0;
  FTYPE abs_err = 0.0;
  FTYPE rel_err = 0.0;
  FTYPE abs_delta = 0.0;
  FTYPE w_e_cyc = Bz*fabs(qd[0])/md[0];
  FTYPE kappa = w_e_cyc/coll_rate[0];
  // FTYPE kappa = 0.0;

  // Echo some parameters
  if (mpi_rank == 0) {
    // cout << "test w_e_cyc = " << w_e_cyc << "\n";
    // cout << "test coll_rate[0] = " << coll_rate[0] << "\n";
    // cout << "test kappa = " << kappa << "\n";
  }

  int xstart=0,xend=nx;
  if (boundary_type[0] == INJECT) {
    if (subdomain.id_number == 0) {
      xstart=1;
      int ix=0;
      int ixp1 = 1;
      int ixm1 = -1;
      FTYPE value_diag;
      for (int iy=0; iy<ny; iy++) {
	int iyp1 = iy+1;
	if (iyp1==ny) iyp1=0;
	int iym1 = iy-1;
	if (iym1<0) iym1=ny-1;
      
	/* LHS */
	
	FTYPE phi_LHS = 0;
	
	/* Diagonal Term */
	/* changed from periodic case */
	value_diag=-((den(ixp1,iy)+den(ix,iy))/(2*dx*dx)+
	   (den(ix,iyp1)+den(ix,iym1)+2*den(ix,iy))/(2*dy*dy));	
	phi_LHS+=phi(ix,iy);

	/*   Forward 1 in x  */
	/* unchanged from periodic case */
	phi_LHS+=((den(ixp1,iy)+den(ix,iy))/(2*dx*dx)
		  +kappa*( den(ix,iyp1)-den(ix,iym1))/(8*dx*dy))
	  *phi(ixp1,iy)/value_diag;

	/*   Backward 1 in x */
	/* on rhs in inject case */

	/*   Forward 1 in y  */
	/* changed from periodic case */
	phi_LHS+=((den(ix,iyp1)+den(ix,iy))/(2*dy*dy)
		  -kappa*(den(ixp1,iy)+den(ix,iy))/(8*dx*dy))
	  *phi(ix,iyp1)/value_diag;

	/*   Backward 1 in y */
	/* changed from periodic case */
	phi_LHS+=((den(ix,iym1)+den(ix,iy))/(2*dy*dy)
		  +kappa*( den(ixp1,iy)+den(ix,iy))/(8*dx*dy))
	  *phi(ix,iym1)/value_diag;

	/*   Diagonal: (ixp1,iyp1) */
	/* unchanged from periodic case */
	phi_LHS+=kappa*(-den(ixp1,iy)+den(ix,iyp1))/(8*dx*dy)
	*phi(ixp1,iyp1)/value_diag;

	/*   Diagonal: (ixp1,iym1) */
	/* unchanged from periodic case */
	phi_LHS+=kappa*(den(ixp1,iy)-den(ix,iym1))/(8*dx*dy)
	  *phi(ixp1,iym1)/value_diag;

	/*   Diagonal (ixm1,iym1)  */
	/* on rhs in inject case */

	/*   Diagonal: (ixm1,iyp1) */
	/* on rhs in inject case */

	/* RHS */
	FTYPE phi_RHS = 0;
	FTYPE Gamma = (1+kappa*kappa)*
	  fspecie[0].m*fspecie[0].nu/(fspecie[0].q*fspecie[0].q);
	/* FTYPE Gamma = (1+kappa*kappa)*
	   fspecie[0].m*fspecie[0].nu/fabs(fspecie[0].q); */
	phi_RHS=
	  Ex0_external*((den(ixp1,iy)-den(ixm1,iy))/(2*dx))+
	  Ey0_external*((den(ix,iyp1)-den(ix,iym1))/(2*dy))+
	  kappa*Ex0_external*((den(ix,iyp1)-den(ix,iym1))/(2*dy))
	  -kappa*Ey0_external*((den(ixp1,iy)-den(ixm1,iy))/(2*dx))
	  +fspecie[0].gamma*fspecie[0].T0/fabs(fspecie[0].q)*
	  ((den(ixp1,iy)+den(ixm1,iy)-2*den(ix,iy))/(dx*dx)+
	   (den(ix,iyp1)+den(ix,iym1)-2*den(ix,iy))/(dy*dy))+
	  +Gamma*
	  ((Jx(ixp1,iy)-Jx(ixm1,iy))/(2*dx)+
	   (Jy(ix,iyp1)-Jy(ix,iym1))/(2*dy));
	
	phi_RHS/=value_diag;

	abs_delta += (phi_LHS-phi_RHS)*(phi_LHS-phi_RHS);
	/* compare */
	if (fabs(phi_LHS-phi_RHS)/phi_LHS > err_tol) {
	  failed++;
	  /* if (subdomain.rank == 0) printf("Failed on %4d  at %4d, %4d with err %10.5g RHS %10.5g LHS %10.5g\n",
					  subdomain.id_number,ix,iy,
					  abs((phi_LHS-phi_RHS)),phi_RHS,phi_LHS);*/
	} 
      }
    } else if (subdomain.id_number == nsubdomains-1) {
      xend = nx -1;
      int ix=nx-1;
      for (int iy=0; iy<ny; iy++) {
	FTYPE value_diag;
	int ixm1 = ix-1;
	int ixp1 = nx;
	int iyp1 = iy+1;
	if (iyp1==ny) iyp1=0;
	int iym1 = iy-1;
	if (iym1<0) iym1=ny-1;
      
	/* LHS */
	
	FTYPE phi_LHS = 0;
	
	/* Diagonal Term */
	/* changed from periodic case */
	value_diag=-((den(ixm1,iy)+den(ix,iy))/(2*dx*dx)+
		     (den(ix,iyp1)+den(ix,iym1)+2*den(ix,iy))/(2*dy*dy));
	phi_LHS+=phi(ix,iy);

	/*   Forward 1 in x  */
	/* moved to rhs in inject */

	/*   Backward 1 in x */
	/* unchanged from periodic case */
	phi_LHS+=((den(ixm1,iy)+den(ix,iy))/(2*dx*dx)
		  +kappa*(-den(ix,iyp1)+den(ix,iym1))/(8*dx*dy))
	  *phi(ixm1,iy)/value_diag;

	/*   Forward 1 in y  */
	/* changed from periodic case */
	phi_LHS+=((den(ix,iyp1)+den(ix,iy))/(2*dy*dy)
		  +kappa*(den(ix,iy)+den(ixm1,iy))/(8*dx*dy))
	  *phi(ix,iyp1)/value_diag;

	/*   Backward 1 in y */
	/* changed from periodic case */
	phi_LHS+=((den(ix,iym1)+den(ix,iy))/(2*dy*dy)
		  -kappa*( den(ix,iy)+den(ixm1,iy))/(8*dx*dy))
	  *phi(ix,iym1)/value_diag;

	/*   Diagonal: (ixp1,iyp1) */
	/* moved to rhs in inject case */

	/*   Diagonal: (ixp1,iym1) */
	/* moved to rhs in inject case */

	/*   Diagonal (ixm1,iym1)  */
	/* unchanged from periodic case */
	phi_LHS+=kappa*(-den(ixm1,iy)+den(ix,iym1))/(8*dx*dy)
	  *phi(ixm1,iym1)/value_diag;

	/*   Diagonal: (ixm1,iyp1) */
	/* unchanged from periodic case */
	phi_LHS+=kappa*(den(ixm1,iy)-den(ix,iyp1))/(8*dx*dy)
	  *phi(ixm1,iyp1)/value_diag;

	/* RHS */
	FTYPE phi_RHS = 0;
	FTYPE Gamma = (1+kappa*kappa)*
	  fspecie[0].m*fspecie[0].nu/(fspecie[0].q*fspecie[0].q);
	/* FTYPE Gamma = (1+kappa*kappa)*
	   fspecie[0].m*fspecie[0].nu/fabs(fspecie[0].q); */
	phi_RHS=
	  Ex0_external*((den(ixp1,iy)-den(ixm1,iy))/(2*dx))+
	  Ey0_external*((den(ix,iyp1)-den(ix,iym1))/(2*dy))+
	  kappa*Ex0_external*((den(ix,iyp1)-den(ix,iym1))/(2*dy))
	  -kappa*Ey0_external*((den(ixp1,iy)-den(ixm1,iy))/(2*dx))
	  +fspecie[0].gamma*fspecie[0].T0/fabs(fspecie[0].q)*
	  ((den(ixp1,iy)+den(ixm1,iy)-2*den(ix,iy))/(dx*dx)+
	   (den(ix,iyp1)+den(ix,iym1)-2*den(ix,iy))/(dy*dy))+
	  +Gamma*
	  ((Jx(ixp1,iy)-Jx(ixm1,iy))/(2*dx)+
	   (Jy(ix,iyp1)-Jy(ix,iym1))/(2*dy));
	
	phi_RHS/=value_diag;

	abs_delta += (phi_LHS-phi_RHS)*(phi_LHS-phi_RHS);
	/* compare */
	if (abs((phi_LHS-phi_RHS))>err_tol) {
	  failed++;
	  /* if (subdomain.rank == 0) printf("Failed on %4d  at %4d, %4d with err %10.5g RHS %10.5g LHS %10.5g\n",
					  subdomain.id_number,ix,iy,
					  abs((phi_LHS-phi_RHS)),phi_RHS,phi_LHS); */
	} else {
	  /*
	  printf("Passed on %4d at %4d, %4d\n",
		 mpi_rank,ix,iy);
	  */
	}
	if ((ix==nx-1) && (iy==ny-1)) {

	}
      }
    }
  }


  /* setup for algorithm 1 */
  for (int ix=xstart; ix<xend; ix++) {
    for (int iy=0; iy<ny; iy++) {
      FTYPE value_diag;
      
      /*
	
	PetscSynchronizedPrintf(MPI_COMM_WORLD,"PROC(%4d) at (%4d,%4d) den = %10.5g\n",mpi_rank,ix,iy,den(ix,iy));
	PetscSynchronizedPrintf(MPI_COMM_WORLD,"PROC(%4d) at (%4d,%4d) Jx = %10.5g\n",mpi_rank,ix,iy,Jx(ix,iy));
	PetscSynchronizedPrintf(MPI_COMM_WORLD,"PROC(%4d) at (%4d,%4d) Jy = %10.5g\n",mpi_rank,ix,iy,Jy(ix,iy));
	PetscSynchronizedPrintf(MPI_COMM_WORLD,"PROC(%4d) at (%4d,%4d) phi = %10.5g\n",mpi_rank,ix,iy,phitmp(ix,iy));
	PetscSynchronizedFlush(MPI_COMM_WORLD);

	*/
	int ixm1 = ix-1;
	int ixp1 = ix+1;
	if (nsubdomains==1) {
	  if (ixm1<0) ixm1=nx-1;
	  if (ixp1==nx) ixp1=0;
	}

	int iyp1 = iy+1;
	if (iyp1==ny) iyp1=0;
	int iym1 = iy-1;
	if (iym1<0) iym1=ny-1;
      
	/* LHS */
	
	FTYPE phi_LHS = 0;
	
	/* Diagonal Term */
	value_diag=-((den(ixp1,iy)+den(ixm1,iy)+2*den(ix,iy))/(2*dx*dx)+
		     (den(ix,iyp1)+den(ix,iym1)+2*den(ix,iy))/(2*dy*dy));
	phi_LHS+=phi(ix,iy);

	/*   Forward 1 in x  */
	phi_LHS+=((den(ixp1,iy)+den(ix,iy))/(2*dx*dx)
		  +kappa*( den(ix,iyp1)-den(ix,iym1))/(8*dx*dy))
	  *phi(ixp1,iy)/value_diag;

	/*   Backward 1 in x */
	phi_LHS+=((den(ixm1,iy)+den(ix,iy))/(2*dx*dx)
		  +kappa*(-den(ix,iyp1)+den(ix,iym1))/(8*dx*dy))
	  *phi(ixm1,iy)/value_diag;

	/*   Forward 1 in y  */
	phi_LHS+=((den(ix,iyp1)+den(ix,iy))/(2*dy*dy)
		  +kappa*(-den(ixp1,iy)+den(ixm1,iy))/(8*dx*dy))
	  *phi(ix,iyp1)/value_diag;

	/*   Backward 1 in y */
	phi_LHS+=((den(ix,iym1)+den(ix,iy))/(2*dy*dy)
		  +kappa*( den(ixp1,iy)-den(ixm1,iy))/(8*dx*dy))
	  *phi(ix,iym1)/value_diag;

	/*   Diagonal: (ixp1,iyp1) */
	phi_LHS+=kappa*(-den(ixp1,iy)+den(ix,iyp1))/(8*dx*dy)
	*phi(ixp1,iyp1)/value_diag;

	/*   Diagonal: (ixp1,iym1) */
	phi_LHS+=kappa*(den(ixp1,iy)-den(ix,iym1))/(8*dx*dy)
	  *phi(ixp1,iym1)/value_diag;

	/*   Diagonal (ixm1,iym1)  */
	phi_LHS+=kappa*(-den(ixm1,iy)+den(ix,iym1))/(8*dx*dy)
	  *phi(ixm1,iym1)/value_diag;

	/*   Diagonal: (ixm1,iyp1) */
	phi_LHS+=kappa*(den(ixm1,iy)-den(ix,iyp1))/(8*dx*dy)
	  *phi(ixm1,iyp1)/value_diag;

	/* RHS */
	FTYPE phi_RHS = 0;
	FTYPE Gamma = (1+kappa*kappa)*
	   fspecie[0].m*fspecie[0].nu/(fspecie[0].q*fspecie[0].q);
	/* FTYPE Gamma = (1+kappa*kappa)*
	   fspecie[0].m*fspecie[0].nu/fabs(fspecie[0].q); */
	phi_RHS=
	  Ex0_external*((den(ixp1,iy)-den(ixm1,iy))/(2*dx))+
	  Ey0_external*((den(ix,iyp1)-den(ix,iym1))/(2*dy))+
	  kappa*Ex0_external*((den(ix,iyp1)-den(ix,iym1))/(2*dy))
	  -kappa*Ey0_external*((den(ixp1,iy)-den(ixm1,iy))/(2*dx))
	  +fspecie[0].gamma*fspecie[0].T0/fabs(fspecie[0].q)*
	  ((den(ixp1,iy)+den(ixm1,iy)-2*den(ix,iy))/(dx*dx)+
	   (den(ix,iyp1)+den(ix,iym1)-2*den(ix,iy))/(dy*dy))+
	  +Gamma*
	  ((Jx(ixp1,iy)-Jx(ixm1,iy))/(2*dx)+
	   (Jy(ix,iyp1)-Jy(ix,iym1))/(2*dy));
	
	phi_RHS/=value_diag;

	abs_delta += (phi_LHS-phi_RHS)*(phi_LHS-phi_RHS);
	/* compare */
	/*
	if (mpi_rank == 0) {
	  cout << "phi_LHS = " << phi_LHS << "\n";
	  cout << "phi_RHS = " << phi_RHS << "\n";
	  cout << "abs error = " << abs((phi_LHS-phi_RHS)) << "\n";
	}
	*/
	abs_err = abs(phi_LHS-phi_RHS);
	// rel_err = abs_err/phi_RHS;
	// rel_err = sqrt(((phi_LHS-phi_RHS)*(phi_LHS-phi_RHS))/(phi_LHS)*(phi_LHS));
	rel_err = abs(phi_LHS-phi_RHS)/(abs(phi_LHS)*abs(phi_RHS));

	// if (mpi_rank == 0) printf("rel_err = %10.5g\n",rel_err);
	if (rel_err>err_tol) {
	  failed++;
	  if (subdomain.rank == 0) printf("Failed on %4d  at %4d, %4d with rel_err %10.5g RHS %10.5g LHS %10.5g\n",
					  subdomain.id_number,ix,iy,
					  rel_err,phi_RHS,phi_LHS);
	} else {
	  // if (subdomain.rank == 0) printf("Passed on %4d at %4d, %4d\n",mpi_rank,ix,iy);
	  
	}
	if ((ix==nx-1) && (iy==ny-1)) {

	}
      }
    }
  FTYPE abs_delta_sum=0;
  MPI_Allreduce(&abs_delta,&abs_delta_sum,1,MPI_FTYPE,MPI_SUM,
		subdomain.neighbor_comm);
  // if (mpi_rank == 0) printf("Abs delta = %g\n",abs_delta_sum);

  return failed;
}
