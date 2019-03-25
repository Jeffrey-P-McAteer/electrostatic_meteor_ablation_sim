#include <stdio.h>
#include <math.h>
#include <complex>
#include <sys/times.h> //Clock related functions
#include "eppic.h"

#include "eppic_fftw.h"

void efield_inject(field &Efield , FArrayND &rho)
//(FArrayND_ranged &phi,FTYPE &rho_dc,FArrayND &rho, int bc)
{

  if (nsubdomains > 1) {
    terminate(1,"efield_inject algorithm is only setup for 1 domain\nTry efield_inject_parallel\n");
  }
  /* Local variables */

  FArrayND_ranged &phi = Efield.phi_rho;  // This is done for convenience only

  FTYPE rho_dc;

  int bc = boundary_type[0];
  
  const static FTYPE fact=-Sqr(dx)/eps;
  int nNy=ny/2+1;
  
  FTYPE auxrhoky0,auxrhokynxm1,auxrhokynxm2;

  static rfftw_plan plan_fwd,plan_bckwd;
  static ArrayNd<FTYPE,1> rinv(nNy), psol(nNy);

  static bool first_entry=TRUE;
  if (first_entry) {
   cout<<mpi_rank<<endl;
    if(bc==0 && mpi_rank == 1) cout<< endl <<"Using 2D Periodic Tridiagonal Solver"<<endl;
    if(bc==1 && mpi_rank == 1) cout<< endl <<"Using 2D Injection Tridiagonal Solver"<<endl;

    //In place ffts. Rho will be modified. 

    int fftw_plan_type = FFTW_MEASURE;
#ifdef DEBUG 
   fftw_plan_type = FFTW_ESTIMATE;
#endif
    plan_fwd=rfftw_create_plan(ny, FFTW_FORWARD, fftw_plan_type);
    plan_bckwd=rfftw_create_plan(ny, FFTW_BACKWARD, fftw_plan_type);

    for(int iy=1;iy<nNy;iy++){
      FTYPE dm=1.+2.*Sqr((dx/dy)*sin(M_PI* (FTYPE)iy /ny));
      FTYPE r=dm+sqrt(Sqr(dm)-1.);
      rinv(iy)=1./r;
      psol(iy)=r*r/(1-r*r);
    }

    first_entry=FALSE;
  }        

  //rho dc part calculation
  rho_dc=0.;
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++) rho_dc+=rho(ix,iy);
  rho_dc/=(nx*ny);
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++) rho(ix,iy)-=rho_dc;

  //Fourier transform of rho in y direction
  FArrayND rho_tmp = rho;
  rfftw(plan_fwd,nx,(fftw_real*) &(rho_tmp(0,0)),1,ny,(fftw_real*) &(rho(0,0)),1,ny);

  //Calculation of phy(x,ky=0). (Special case). Real values
  //no need of dealing with the imaginary part.
  /* Actually how it works:
     After reading the section on pg. 318 of Birdsall (14-6, Poison's equation 
     solutions for systems bounded in x and Periodic in y), read Appendix D. 
     There he discusses the difference between a periodic solution and a 
     bounded electrodes solution to Poisson's eqn. 

     In problem Da, at the end of appendix, details of how to deal with the 
     case where dm (see equation D.13) equals 1, which is the case for either
     the periodic or non-periodic (in-x) and periodic in-y 2D case, when 
     solving for the iy=0 term. 
     
     Coded below is the correct method for dealing with either the periodic
     or non-periodic cases (where the potential is set to zero at x=-1,Nx.
  */
  auxrhoky0=rho(0,0);
  auxrhokynxm1=rho(nx-1,0);
  
  rho(0,0)=0.0;
  rho(nx-1,0)=(rho(1,0)*(FTYPE)(1-bc)+auxrhoky0*(FTYPE)bc)*fact; 
  
  for(int ix=2-bc;ix<nx-1;ix++) {
    rho(nx-1,0)+=rho(ix,0)*(FTYPE)(ix+bc)*fact;  
  }
  rho(nx-1,0)+=auxrhokynxm1*(FTYPE)(nx-1+bc)*fact;  
  
  rho(nx-1,0)/=(-nx-bc);

  auxrhokynxm2=rho(nx-2,0);
  rho(nx-2,0)=auxrhokynxm1*fact+2.*rho(nx-1,0); 

  for(int ix=nx-3;ix>=1-bc;ix--){
      auxrhokynxm1=rho(ix,0);
      rho(ix,0)=auxrhokynxm2*fact+2.*rho(ix+1,0)-rho(ix+2,0);
      auxrhokynxm2=auxrhokynxm1;
  }
  
  //Calculation of phi(x,ky!=0) using the Birdsall method (See book pag. 318)

  for(int iy=1;iy<nNy-1;iy++){
    rho(nx-1,iy)*=(fact*rinv(iy)); //real part   
    rho(nx-1,ny-iy)*=(fact*rinv(iy)); //imaginary part. See fftw manual.   
    
    for(int ix=nx-2;ix>=0;ix--){
      rho(ix,iy)=(rho(ix,iy)*fact+rho(ix+1,iy))*rinv(iy); // real part
      rho(ix,ny-iy)=(rho(ix,ny-iy)*fact+rho(ix+1,ny-iy))*rinv(iy); //imaginary part. See fftw manual.
    }   
    rho(0,iy)=rho(0,iy)*psol(iy); // real part
    rho(0,ny-iy)=rho(0,ny-iy)*psol(iy); //imaginary part. See fftw manual.
    for(int ix=1;ix<nx;ix++){
      rho(ix,iy)=rho(ix-1,iy)*rinv(iy)-rho(ix,iy); //real part.
      rho(ix,ny-iy)=rho(ix-1,ny-iy)*rinv(iy)-rho(ix,ny-iy); //imaginary part. See fftw manual.
    }
  }
  // The last point of the phiky(ix,ky=nNy-1) has to be calculated outside to avoid 
  // twice calculation inside the loop. Just real values. No need of dealing with imaginary part.
  rho(nx-1,nNy-1)*=(fact*rinv(nNy-1));    
  for(int ix=nx-2;ix>=0;ix--) rho(ix,nNy-1)=(rho(ix,nNy-1)*fact+rho(ix+1,nNy-1))*rinv(nNy-1);
  rho(0,nNy-1)=rho(0,nNy-1)*psol(nNy-1);
  for(int ix=1;ix<nx;ix++) rho(ix,nNy-1)=rho(ix-1,nNy-1)*rinv(nNy-1)-rho(ix,nNy-1);
  
  //Calculation of phi(x,y). This inverse trnsform has not been normalized by ny. This is 
  //taken care when phi_dc is calculated.
  rho_tmp = rho;
  rfftw(plan_bckwd,nx,(fftw_real*) &(rho_tmp(0,0)),1,ny,(fftw_real*) &(rho(0,0)),1,ny);
    
  //phi dc calculation
  rho_dc=0.;
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++){
    rho_dc+=rho(ix,iy);
    phi(ix,iy)=rho(ix,iy);
  }
  rho_dc/=(nx*ny);

  //The division of ny is for normalization of the IFFT.
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++) phi(ix,iy)=(phi(ix,iy)-rho_dc)/ny;

}
