#include <stdio.h>
#include <math.h>
#include <complex>
#include <sys/times.h> //Clock related functions
#include "eppic.h"

#include "eppic_fftw.h"

/* -------------------------------------------------------------------------
   helper routines used in the main routine: 
   ------------------------------------------------------------------------- */

void efield_inject_parallel_first_entry(rfftw_plan &plan_fwd,
					rfftw_plan &plan_bckwd,
					ArrayNd<FTYPE,1> &rinv,
					ArrayNd<FTYPE,1> &psol,
					int bc
					);

void efield_inject_parallel_ky0(FArrayND &rho, int bc);

void efield_inject_parallel_kyi_forward(FArrayND &rho, ArrayNd<FTYPE,1> &rinv);

void efield_inject_parallel_kyi_backward(FArrayND &rho,
					 ArrayNd<FTYPE,1> &rinv,
					 ArrayNd<FTYPE,1> &psol);

void efield_inject_parallel_sync_subdomain(FArrayND &rho);

void efield_inject_parallel_kylast(FArrayND &rho,
				   ArrayNd<FTYPE,1> &rinv,
				   ArrayNd<FTYPE,1> &psol);

void efield_inject_parallel_tophi(FArrayND &rho,FArrayND_ranged &phi); 


/* -------------------------------------------------------------------------
   Main routine: efield_inject_parallel
   ------------------------------------------------------------------------- */

void efield_inject_parallel(field &Efield , FArrayND &rho)
//(FArrayND_ranged &phi,FTYPE &rho_dc,FArrayND &rho, int bc)
{

  /* Local variables */

  FArrayND_ranged &phi = Efield.phi_rho;  // This is done for convenience only

  FTYPE rho_dc;

  int bc = boundary_type[0];
  

  static rfftw_plan plan_fwd,plan_bckwd;
  static ArrayNd<FTYPE,1> rinv(ny/2+1), psol(ny/2+1);

  static bool first_entry=TRUE;
  if (first_entry) 
    efield_inject_parallel_first_entry(plan_fwd,plan_bckwd,rinv,psol,bc);

  /* rho dc part calculation */
  rho_dc=0.;
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++) rho_dc+=rho(ix,iy);
  rho_dc/=(nx*ny);
  if (nsubdomains>1)  {
    FTYPE rho_dc_tmp=rho_dc;
    MPI_Allreduce(&rho_dc_tmp,&rho_dc,1,MPI_FTYPE,MPI_SUM,
		  subdomain.neighbor_comm);
    rho_dc/=nsubdomains;
  }
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++) rho(ix,iy)-=rho_dc;


  /* Fourier transform of rho in y direction */
  FArrayND rho_tmp = rho;
  rfftw(plan_fwd,nx,(fftw_real*) &(rho_tmp(0,0)),1,ny,
	(fftw_real*) &(rho(0,0)),1,ny);

  
  /* Calculation of phy(x,ky=0). (Special case). Real values */
  efield_inject_parallel_ky0(rho,bc);

  /* now do kyi>0 except kyi==nNy-1 */
  efield_inject_parallel_kyi_forward(rho,rinv);
  efield_inject_parallel_kyi_backward(rho,rinv,psol);
  /* if subdomain.np > 1, pass ny_subdomain rows to other subdomains */
  if (subdomain.np>1) efield_inject_parallel_sync_subdomain(rho); 

  /* another special case, to avoid double counting */
  efield_inject_parallel_kylast(rho,rinv,psol);

  /* inverse fft */
  rho_tmp = rho;
  rfftw(plan_bckwd,nx,(fftw_real*) &(rho_tmp(0,0)),1,ny,
	(fftw_real*) &(rho(0,0)),1,ny);

  efield_inject_parallel_tophi(rho,phi);



  if(first_entry) first_entry=FALSE;
}

/* ---------------- BEGIN HELPER ROUTINES ---------------------------------*/



void efield_inject_parallel_first_entry(rfftw_plan &plan_fwd,
					rfftw_plan &plan_bckwd,
					ArrayNd<FTYPE,1> &rinv,
					ArrayNd<FTYPE,1> &psol,
					int bc
					)
{
  int nNy=ny/2+1;
  if(bc==0 && mpi_rank == 0) 
    cout<< endl <<";Using 2D Periodic Tridiagonal Solver"<<endl;
  if(bc==1 && mpi_rank == 0) 
    cout<< endl <<";Using 2D Injection Tridiagonal Solver"<<endl;
  
  //In place ffts. Rho will be modified. 
  plan_fwd=rfftw_create_plan(ny,FFTW_FORWARD,FFTW_MEASURE);// | FFTW_IN_PLACE);
  plan_bckwd=rfftw_create_plan(ny,FFTW_BACKWARD,FFTW_MEASURE);// | FFTW_IN_PLACE);

  for(int iy=1;iy<nNy;iy++){
    FTYPE dm=1.+2.*Sqr((dx/dy)*sin(M_PI* (FTYPE)iy /ny));
    FTYPE r=dm+sqrt(Sqr(dm)-1.);
    rinv(iy)=1./r;
    psol(iy)=r*r/(1-r*r);
  }
}



void efield_inject_parallel_ky0(FArrayND &rho,int bc)
{
  //no need of dealing with the imaginary part.
  /* need to pass to last domain from first, if non-periodic */
  FTYPE auxrhoky0,auxrhokynxm1,auxrhokynxm2;
  const static FTYPE fact=-Sqr(dx)/eps;

  if (subdomain.id_number == 0)
    if (bc) 
      auxrhoky0=rho(0,0);
    else 
      auxrhoky0=rho(1,0);

  MPI_Request request1,request2;
  MPI_Status  status1;
  if (nsubdomains>1){
    /* non-blocking send */
    if (subdomain.id_number == 0) 
	MPI_Isend(&auxrhoky0,1, MPI_FTYPE,nsubdomains-1, 2,
		  subdomain.neighbor_comm, &request1);

    if (subdomain.id_number == nsubdomains-1)
      MPI_Irecv(&auxrhoky0, 1, MPI_FTYPE,0, 2,
	       subdomain.neighbor_comm, &request1);

    /* check that send worked */
    if ((subdomain.id_number==0)||(subdomain.id_number==nsubdomains-1))
      MPI_Wait(&request1, &status1);
  }

  /* only last domain needs this */
  if (subdomain.id_number == nsubdomains-1)
    auxrhokynxm1=rho(nx-1,0);
  
  /*
  if (subdomain.rank == 0) 
    printf("Domain: %d auxrhokynxm1: %16.8g\n",subdomain.id_number,auxrhokynxm1);
  */
  FTYPE rho_ky0=0;
  
  /* 0,0 */
  if (subdomain.id_number == 0) rho(0,0)=0.0;

  /* nx-1,0 */
  if (subdomain.id_number == nsubdomains-1) 
    rho(nx-1,0)=auxrhoky0*fact; 

    /* 
       This is the old line: 
      
       rho(nx-1,0)=(rho(1,0)*(FTYPE)(1-bc)+auxrhoky0*(FTYPE)bc)*fact; 

       auxrhoky0 is now defined above based on bc, equivilant outcome, but 
       variable name needs changing
    */

  /* if d.decomposed collect sum in separate var, pass to last domain */
  for(int ix=0;ix<nx-1;ix++) 
    rho_ky0+=rho(ix,0)*(FTYPE)(ix+nx*subdomain.id_number+bc)*fact;  
  if (subdomain.id_number<nsubdomains-1)
    rho_ky0+=rho(nx-1,0)*(FTYPE)(nx-1+nx*subdomain.id_number+bc)*fact;  
  /* if first domain, added too many terms, which were taken care of above */
  if ((!bc) && (subdomain.id_number == 0)) /* only matters in periodic case */
    rho_ky0 -= rho(1,0)*fact;  
  if (nsubdomains > 1) {
    FTYPE rho_ky0_total=0;
    MPI_Reduce(&rho_ky0,&rho_ky0_total,1,MPI_FTYPE,MPI_SUM,
	       nsubdomains-1,subdomain.neighbor_comm);
    rho_ky0 = rho_ky0_total;
  }
  if (subdomain.id_number == nsubdomains-1) {
    rho(nx-1,0)+=rho_ky0;
    rho(nx-1,0)+=auxrhokynxm1*(FTYPE)(nx-1+nx*subdomain.id_number+bc)*fact;  
    rho(nx-1,0)/=(-nx*nsubdomains-bc);
  }

  /* nx-2,0 */
  if (subdomain.id_number == nsubdomains-1) {
    auxrhokynxm2=rho(nx-2,0);
    rho(nx-2,0)=auxrhokynxm1*fact+2.*rho(nx-1,0); 
  }


  /* if d.decomposed loop over domains before ix */
  FTYPE *last = new FTYPE[3];
  if (subdomain.id_number==nsubdomains-1) {
    for(int ix=nx-3;ix>=1;ix--){
      auxrhokynxm1=rho(ix,0);
      rho(ix,0)=auxrhokynxm2*fact+2.*rho(ix+1,0)-rho(ix+2,0);
      auxrhokynxm2=auxrhokynxm1;
    }
    if ((subdomain.id_number > 0) || (bc)) {
      auxrhokynxm1=rho(0,0);
      rho(0,0)=auxrhokynxm2*fact+2.*rho(1,0)-rho(2,0);
      auxrhokynxm2=auxrhokynxm1;
    }
    if (nsubdomains>1) {
      last[0] = auxrhokynxm2;
      last[1] = rho(0,0);
      last[2] = rho(1,0);
      MPI_Isend(last,3, MPI_FTYPE,nsubdomains-2, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);      
    }
  } else {
    MPI_Irecv(last, 3, MPI_FTYPE,subdomain.id_number+1, 2,
	      subdomain.neighbor_comm, &request2);
    MPI_Wait(&request2, &status1);
    auxrhokynxm1=rho(nx-1,0);
    rho(nx-1,0)=last[0]*fact+2.*last[1]-last[2];
    auxrhokynxm2=auxrhokynxm1;
    auxrhokynxm1=rho(nx-2,0);
    rho(nx-2,0)=auxrhokynxm2*fact+2.*rho(nx-1,0)-last[1];
    auxrhokynxm2=auxrhokynxm1;
    for(int ix=nx-3;ix>=1;ix--){
      auxrhokynxm1=rho(ix,0);
      rho(ix,0)=auxrhokynxm2*fact+2.*rho(ix+1,0)-rho(ix+2,0);
      auxrhokynxm2=auxrhokynxm1;
    }
    if ((subdomain.id_number > 0) || (bc)) {
      auxrhokynxm1=rho(0,0);
      rho(0,0)=auxrhokynxm2*fact+2.*rho(1,0)-rho(2,0);
      auxrhokynxm2=auxrhokynxm1;
    }
    if (subdomain.id_number>0) {
      last[0] = auxrhokynxm2;
      last[1] = rho(0,0);
      last[2] = rho(1,0);
      MPI_Isend(last,3, MPI_FTYPE,subdomain.id_number-1, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);      
    }
  }
  delete [] last;
}

void efield_inject_parallel_kyi_forward(FArrayND &rho,
					ArrayNd<FTYPE,1> &rinv) 
{
  MPI_Request request1,request2;
  MPI_Status  status1;
  int nNy=ny/2+1;
  const static FTYPE fact=-Sqr(dx)/eps;
  int ny_subdomain = 2*static_cast<int>(ceil(FTYPE(nNy-1)/FTYPE(subdomain.np)));
  int iy_start = 1+ny_subdomain*subdomain.rank/2;
  int iy_end = iy_start+ny_subdomain/2;
  if (iy_end>nNy-1) iy_end = nNy-1;

  FTYPE *last = new FTYPE[ny_subdomain];
  //Calculation of phi(x,ky!=0) using the Birdsall method (See book pag. 318)
  /* if d.decomposed loop over domains before iy and ix */
  /* forward */
  if ((subdomain.id_number == nsubdomains-1)){
    for(int iy=iy_start;iy<iy_end;iy++){
      rho(nx-1,iy)*=(fact*rinv(iy)); //real part   
      rho(nx-1,ny-iy)*=(fact*rinv(iy)); //imaginary part. See fftw manual.   
      for(int ix=nx-2;ix>=0;ix--){
	rho(ix,iy)=(rho(ix,iy)*fact+rho(ix+1,iy))*rinv(iy); // real part
	rho(ix,ny-iy)=(rho(ix,ny-iy)*fact+rho(ix+1,ny-iy))*rinv(iy); //imaginary part. See fftw manual.
      }   
      last[iy-iy_start] = rho(0,iy);
      last[ny_subdomain-1-(iy-iy_start)] = rho(0,ny-iy);
    }
    if (nsubdomains>1) {
      MPI_Isend(last,ny_subdomain, MPI_FTYPE,nsubdomains-2, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);
    }
  } else {
    MPI_Irecv(last, ny_subdomain, MPI_FTYPE,subdomain.id_number+1, 2,
	     subdomain.neighbor_comm, &request2);
    MPI_Wait(&request2, &status1);
    for(int iy=iy_start;iy<iy_end;iy++){
      rho(nx-1,iy)=(rho(nx-1,iy)*fact+last[iy-iy_start])*rinv(iy); // real part
      rho(nx-1,ny-iy)=(rho(nx-1,ny-iy)*fact+last[(ny_subdomain-1)-(iy-iy_start)])*rinv(iy); //imaginary part. See fftw manual.
	
      for(int ix=nx-2;ix>=0;ix--){
	rho(ix,iy)=(rho(ix,iy)*fact+rho(ix+1,iy))*rinv(iy); // real part
	rho(ix,ny-iy)=(rho(ix,ny-iy)*fact+rho(ix+1,ny-iy))*rinv(iy); //imaginary part. See fftw manual.
      }   
      
      last[iy-iy_start] = rho(0,iy);
      last[ny_subdomain-1-(iy-iy_start)] = rho(0,ny-iy);
    }
    if (subdomain.id_number>0) {
      MPI_Isend(last,ny_subdomain, MPI_FTYPE,subdomain.id_number-1, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);
    }
  }

  delete []last;
}

void efield_inject_parallel_kyi_backward(FArrayND &rho,
					 ArrayNd<FTYPE,1> &rinv,
					 ArrayNd<FTYPE,1> &psol)
{
  MPI_Request request1,request2;
  MPI_Status  status1;
  int nNy=ny/2+1;
  int ny_subdomain = 2*static_cast<int>(ceil(FTYPE(nNy-1)/FTYPE(subdomain.np)));
  int iy_start = 1+ny_subdomain*subdomain.rank/2;
  int iy_end = iy_start+ny_subdomain/2;
  if (iy_end>nNy-1) iy_end = nNy-1;


  FTYPE *last = new FTYPE[ny_subdomain];
  /* back again */
  if (subdomain.id_number == 0) {
    for(int iy=iy_start;iy<iy_end;iy++){
      rho(0,iy)=rho(0,iy)*psol(iy); // real part
      rho(0,ny-iy)=rho(0,ny-iy)*psol(iy); //imaginary part. See fftw manual.
      for(int ix=1;ix<nx;ix++){
	rho(ix,iy)=rho(ix-1,iy)*rinv(iy)-rho(ix,iy); //real part.
	rho(ix,ny-iy)=rho(ix-1,ny-iy)*rinv(iy)-rho(ix,ny-iy); //imaginary part. See fftw manual.
      }
      last[iy-iy_start] = rho(nx-1,iy);
      last[ny_subdomain-1-(iy-iy_start)] = rho(nx-1,ny-iy);
    }
    if (nsubdomains>1) {
      MPI_Isend(last,ny_subdomain, MPI_FTYPE,1, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);
    }
  } else {
    MPI_Irecv(last, ny_subdomain, MPI_FTYPE,subdomain.id_number-1, 2,
	      subdomain.neighbor_comm, &request2);
    MPI_Wait(&request2, &status1);
    for(int iy=iy_start;iy<iy_end;iy++){
      rho(0,iy)=last[iy-iy_start]*rinv(iy)-rho(0,iy); //real part.
      rho(0,ny-iy)=last[ny_subdomain-1-(iy-iy_start)]*rinv(iy)-rho(0,ny-iy); //imaginary part. See fftw manual.
      for(int ix=1;ix<nx;ix++){
	rho(ix,iy)=rho(ix-1,iy)*rinv(iy)-rho(ix,iy); //real part.
	rho(ix,ny-iy)=rho(ix-1,ny-iy)*rinv(iy)-rho(ix,ny-iy); //imaginary part. See fftw manual.
      }
      last[iy-iy_start] = rho(nx-1,iy);
      last[ny_subdomain-1-(iy-iy_start)] = rho(nx-1,ny-iy);
    }
    if (subdomain.id_number<nsubdomains-1) {
      MPI_Isend(last,ny_subdomain, MPI_FTYPE,subdomain.id_number+1, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);
    }
  }

  delete []last;
}

/** Passes nx by ny_subdomain array to other procs in subdomain and assigns to
    correct portion of rho. 
    
    This is used if there are multiple processors per subdomain. In this case, 
    the efield_inject_parallel routines divide the work calculating each iy 
    row amongst the subdomain processors, and this info needs to be 
    interchanged before the calculation finishes.

    \param rho a FArrayND to be modified (input and output)
*/

void efield_inject_parallel_sync_subdomain(FArrayND &rho)
{
  if (subdomain.np==1) return;
  MPI_Barrier(MPI_COMM_WORLD);
  int nNy=ny/2+1;
  int ny_subdomain = 2*static_cast<int>(ceil(FTYPE(nNy-1)/FTYPE(subdomain.np)));
  FArrayND rho_tmp = FArrayND(INDICIES(nx,ny_subdomain,nz));
  for(int isubdomain=0;isubdomain<subdomain.np;isubdomain++){
    int iy_start = 1+ny_subdomain*isubdomain/2;
    int iy_end = iy_start+ny_subdomain/2;
    if (iy_end>nNy-1) iy_end = nNy-1;
    if (subdomain.rank == isubdomain) {
      rho_tmp=0;
      for (int ix=0;ix<nx;ix++) 
	for (int iy=iy_start;iy<iy_end;iy++) 
	  for (int iz=0;iz<nz;iz++) {
	    rho_tmp(INDICIES(ix,iy-iy_start,iz)) = rho(INDICIES(ix,iy,iz));
	    rho_tmp(INDICIES(ix,ny_subdomain-1-(iy-iy_start),iz)) =
	      rho(INDICIES(ix,ny-iy,iz));
	  }
    }
    MPI_Bcast(
	      &(rho_tmp(INDICIES(0,0,0))),
	      rho_tmp.size(), 
	      MPI_FTYPE,
	      isubdomain,
	      subdomain.internal_comm);
    if (subdomain.rank != isubdomain) {
      for (int ix=0;ix<nx;ix++) 
	for (int iy=iy_start;iy<iy_end;iy++) 
	  for (int iz=0;iz<nz;iz++) {
	    rho(INDICIES(ix,iy,iz)) = rho_tmp(INDICIES(ix,iy-iy_start,iz));
	    rho(INDICIES(ix,ny-iy,iz)) =
	      rho_tmp(INDICIES(ix,ny_subdomain-1-(iy-iy_start),iz));
	    
	  }
    }
  }
}

void efield_inject_parallel_kylast(FArrayND &rho,
				   ArrayNd<FTYPE,1> &rinv,
				   ArrayNd<FTYPE,1> &psol)
{
  MPI_Request request1,request2;
  MPI_Status  status1;
  int nNy=ny/2+1;
  const static FTYPE fact=-Sqr(dx)/eps;
  FTYPE *last = new FTYPE;
  if (subdomain.id_number == nsubdomains-1)
    rho(nx-1,nNy-1)*=(fact*rinv(nNy-1));    
  if (nsubdomains>1) {
    if (subdomain.id_number == nsubdomains-1) {
      for(int ix=nx-2;ix>=0;ix--) rho(ix,nNy-1)=(rho(ix,nNy-1)*fact+rho(ix+1,nNy-1))*rinv(nNy-1);
      *last = rho(0,nNy-1);
      MPI_Isend(last,1, MPI_FTYPE,subdomain.id_number-1, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);
    } else if (subdomain.id_number > 0) {
      MPI_Irecv(last, 1, MPI_FTYPE,subdomain.id_number+1, 2,
	       subdomain.neighbor_comm, &request2);
      MPI_Wait(&request2, &status1);
      rho(nx-1,nNy-1)=(rho(nx-1,nNy-1)*fact+*last)*rinv(nNy-1);
      for(int ix=nx-2;ix>=0;ix--) rho(ix,nNy-1)=(rho(ix,nNy-1)*fact+rho(ix+1,nNy-1))*rinv(nNy-1);
      *last  = rho(0,nNy-1);
      MPI_Isend(last,1, MPI_FTYPE,subdomain.id_number-1, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);
    } else {
      MPI_Irecv(last, 1, MPI_FTYPE,subdomain.id_number+1, 2,
		subdomain.neighbor_comm, &request2);
      MPI_Wait(&request2, &status1);
      rho(nx-1,nNy-1)=(rho(nx-1,nNy-1)*fact+*last)*rinv(nNy-1);
      for(int ix=nx-2;ix>=0;ix--) rho(ix,nNy-1)=(rho(ix,nNy-1)*fact+rho(ix+1,nNy-1))*rinv(nNy-1);
    }
  } else for(int ix=nx-2;ix>=0;ix--) rho(ix,nNy-1)=(rho(ix,nNy-1)*fact+rho(ix+1,nNy-1))*rinv(nNy-1);

  if (subdomain.id_number == 0) 
    rho(0,nNy-1)=rho(0,nNy-1)*psol(nNy-1);
  if (nsubdomains>1) {
    if (subdomain.id_number == 0) {
      for(int ix=1;ix<nx;ix++) rho(ix,nNy-1)=rho(ix-1,nNy-1)*rinv(nNy-1)-rho(ix,nNy-1);
      *last = rho(nx-1,nNy-1);
      MPI_Isend(last,1, MPI_FTYPE,subdomain.id_number+1, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);
    } else if (subdomain.id_number<nsubdomains-1) {
      MPI_Irecv(last, 1, MPI_FTYPE,subdomain.id_number-1, 2,
	       subdomain.neighbor_comm, &request2);
      MPI_Wait(&request2, &status1);
      rho(0,nNy-1)=*last*rinv(nNy-1)-rho(0,nNy-1);
      for(int ix=1;ix<nx;ix++) rho(ix,nNy-1)=rho(ix-1,nNy-1)*rinv(nNy-1)-rho(ix,nNy-1);
      *last = rho(nx-1,nNy-1);
      MPI_Isend(last,1, MPI_FTYPE,subdomain.id_number+1, 2,
		subdomain.neighbor_comm, &request1);
      MPI_Wait(&request1, &status1);
    } else {
      MPI_Irecv(last, 1, MPI_FTYPE,subdomain.id_number-1, 2,
	       subdomain.neighbor_comm, &request2);
      MPI_Wait(&request2, &status1);
      rho(0,nNy-1)=*last*rinv(nNy-1)-rho(0,nNy-1);
      for(int ix=1;ix<nx;ix++) rho(ix,nNy-1)=rho(ix-1,nNy-1)*rinv(nNy-1)-rho(ix,nNy-1);

    }
  } else for(int ix=1;ix<nx;ix++) rho(ix,nNy-1)=rho(ix-1,nNy-1)*rinv(nNy-1)-rho(ix,nNy-1);
  if (subdomain.id_number<nsubdomains-1) MPI_Wait(&request1, &status1);
  delete last;
}


void efield_inject_parallel_tophi(FArrayND &rho,FArrayND_ranged &phi)
{
  /* phi dc calculation */
  FTYPE rho_dc=0.;
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++){
    rho_dc+=rho(ix,iy);
    phi(ix,iy)=rho(ix,iy);
  }

  if (nsubdomains>1) {
    FTYPE rho_dc_all=0;
    MPI_Allreduce(&rho_dc,&rho_dc_all,1,MPI_FTYPE,MPI_SUM,
		  subdomain.neighbor_comm);
    rho_dc = rho_dc_all/nsubdomains;
  }

  rho_dc/=(nx*ny);

  /* The division of ny is for normalization of the IFFT. */
  for(int ix=0;ix<nx;ix++) 
    for(int iy=0;iy<ny;iy++) 
      phi(ix,iy)=(phi(ix,iy)-rho_dc)/ny;

  MPI_Request request_west, request_east;
  MPI_Status status_west, status_east;

  /* pass west */
  /* if first domain, set to zero, no pass */
  int size_send=phi.length()/phi.size(0);
  const int index0[]={INDICIES(0,0,0)};
  if (subdomain.id_number > 0) 
    MPI_Isend(phi.address(index0),size_send*phix_guard_size[1],MPI_FTYPE,
	      subdomain.id_number-1,2,subdomain.neighbor_comm,&request_west);
  else 
    for(int ix=-((int)phix_guard_size[0]); ix<0;ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++)
	  phi(INDICIES(ix,iy,iz)) = 0;
  
  /* if last domain, no receive */
  const int index1[]={INDICIES(nx,0,0)};
  if (subdomain.id_number < nsubdomains-1) 
    MPI_Irecv(phi.address(index1),size_send*phix_guard_size[1],MPI_FTYPE,
	      subdomain.id_number+1,2,subdomain.neighbor_comm,&request_west);

  MPI_Wait(&request_west,&status_west);

  /* pass east */
  /* if last domain, set to zero, no pass */
  const int index2[]={INDICIES(nx-phix_guard_size[0],0,0)};
  if (subdomain.id_number < nsubdomains-1) 
    MPI_Isend(phi.address(index2),size_send*phix_guard_size[0],MPI_FTYPE,
	      subdomain.id_number+1,3,subdomain.neighbor_comm,&request_east);
  else
    for (int ix=nx;ix<(nx+int(phix_guard_size[1])); ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++)
	  phi(INDICIES(ix,iy,iz)) = 0;
  
  /* if first domain, no receive */
  const int index3[]={INDICIES(-(int(phix_guard_size[0])),0,0)};
  if (subdomain.id_number > 0) 
    MPI_Irecv(phi.address(index3),size_send*phix_guard_size[0],MPI_FTYPE,
	      subdomain.id_number-1,3,subdomain.neighbor_comm,&request_east);

  
  MPI_Wait(&request_east,&status_east);

}
