#include <stdio.h>
#include <math.h>
#include <complex>
#include <sys/times.h> //Clock related functions
#include "eppic.h"

#include "eppic_fftw.h"

void efield_inject_orig(field &Efield,FArrayND &rho, int bc)
{

  const static FTYPE fact=-Sqr(dx)/eps;
  int nNy=ny/2+1;

  FArrayND_ranged &phi = Efield.phi_rho;

  FTYPE auxrhoky0,auxrhokynxm1;

  static rfftw_plan plan_fwd,plan_bckwd;
  static ArrayNd<FTYPE,1> rinv(nNy), psol(nNy);

  static bool first_entry=true;
  if (first_entry) {
    if(bc==0 && mpi_rank == 1) cout<< endl <<"Using 2D Periodic Tridiagonal Solver"<<endl;
    if(bc==1 && mpi_rank == 1) cout<< endl <<"Using 2D Injection Tridiagonal Solver"<<endl;

    //In place ffts. Rho will be modified. 
    plan_fwd=rfftw_create_plan(ny,FFTW_FORWARD,FFTW_MEASURE | FFTW_IN_PLACE);
    plan_bckwd=rfftw_create_plan(ny,FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);

    for(int iy=1;iy<nNy;iy++){
      FTYPE dm=1.+2.*Sqr((dx/dy)*sin(M_PI* (FTYPE)iy /ny));
      FTYPE r=dm+sqrt(Sqr(dm)-1.);
      rinv(iy)=1./r;
      psol(iy)=r*r/(1-r*r);
    }

    first_entry=false;
  }        

  //rho dc part calculation
  FTYPE rho_dc=0.;
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++){
     rho_dc+=rho(ix,iy);
     phi(ix,iy)=rho(ix,iy);
  }
  rho_dc/=(nx*ny);
  for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++) phi(ix,iy)-=rho_dc;
  
 //Fourier transform of rho in y direction
 rfftw(plan_fwd,nx,(fftw_real*) &(phi(0,0)),1,ny,(fftw_real*) &(phi(0,0)),1,ny);
  //for(int ix=-1;ix<nx+2;ix++) cout<<"rho "<<ix<<" "<<phi(ix,0)<<endl;
  
  //Calculation of phi(x,ky=0). (Special case). With phi(0)=0 and E(0)=0.
  auxrhoky0=phi(0,0);
  auxrhokynxm1=phi(1,0);
  phi(0,0)=0.;
  phi(1,0)=auxrhoky0*fact/2.;
  phi(-1,0)=phi(1,0);
  //phi(0,0)=0.;
  auxrhoky0=phi(2,0);
  phi(2,0)=auxrhokynxm1*fact+2.*phi(1,0);
  
  for(int ix=3;ix<nx+1;ix++){
    auxrhokynxm1=auxrhoky0;
    auxrhoky0=phi(ix,0);
    phi(ix,0)=auxrhokynxm1*fact+2.*phi(ix-1,0)-phi(ix-2,0);
  }
    phi(nx+1,0)=2.*phi(nx,0)-phi(nx-1,0);
  //for(int ix=-1;ix<nx+2;ix++) cout<<ix<<" "<<phi(ix,0)<<endl;

  //Calculation of phi(x,ky=0). (Special case). Real values
  //no need of dealing with the imaginary part.
  /*auxrhoky0=phi(0,0);
  auxrhokynxm1=phi(nx-1,0);
  
  phi(0,0)=0.0;
  phi(nx-1,0)=(phi(1,0)*(FTYPE)(1-bc)+auxrhoky0*(FTYPE)bc)*fact; 
  
  for(int ix=2-bc;ix<nx-1;ix++) phi(nx-1,0)+=phi(ix,0)*(FTYPE)(ix+bc)*fact;  
  phi(nx-1,0)+=auxrhokynxm1*(FTYPE)(nx-1+bc)*fact;  
  
  phi(nx-1,0)/=(-nx-bc);

  auxrhokynxm2=phi(nx-2,0);
  phi(nx-2,0)=auxrhokynxm1*fact+2.*phi(nx-1,0); 

  for(int ix=nx-3;ix>=1-bc;ix--){
      auxrhokynxm1=phi(ix,0);
      phi(ix,0)=auxrhokynxm2*fact+2.*phi(ix+1,0)-phi(ix+2,0);
      auxrhokynxm2=auxrhokynxm1;
  }*/
  
  //Calculation of phi(x,ky!=0) using the Birdsall method (See book pag. 318)
  for(int iy=1;iy<nNy-1;iy++){
    phi(nx-1,iy)*=(fact*rinv(iy)); //real part   
    phi(nx-1,ny-iy)*=(fact*rinv(iy)); //imaginary part. See fftw manual.   
    
    for(int ix=nx-2;ix>=0;ix--){
      phi(ix,iy)=(phi(ix,iy)*fact+phi(ix+1,iy))*rinv(iy); // real part
      phi(ix,ny-iy)=(phi(ix,ny-iy)*fact+phi(ix+1,ny-iy))*rinv(iy); //imaginary part. See fftw manual.
    }   
    phi(0,iy)=phi(0,iy)*psol(iy); // real part
    phi(0,ny-iy)=phi(0,ny-iy)*psol(iy); //imaginary part. See fftw manual.
    for(int ix=1;ix<nx;ix++){
      phi(ix,iy)=phi(ix-1,iy)*rinv(iy)-phi(ix,iy); //real part.
      phi(ix,ny-iy)=phi(ix-1,ny-iy)*rinv(iy)-phi(ix,ny-iy); //imaginary part. See fftw manual.
    }
  phi(-1,iy)=phi(0,iy)*rinv(iy);
  phi(-1,ny-iy)=phi(0,ny-iy)*rinv(iy);  
  phi(nx,iy)=phi(nx-1,iy)*rinv(iy);
  phi(nx,ny-iy)=phi(nx-1,ny-iy)*rinv(iy);  
  phi(nx+1,iy)=phi(nx,iy)*rinv(iy);
  phi(nx+1,ny-iy)=rinv(iy)*phi(nx,ny-iy);  
  }
  
  // The last point of the phiky(ix,ky=nNy-1) has to be calculated outside to avoid 
  // twice calculation inside the loop. Just real values. No need of dealing with imaginary part.
  phi(nx-1,nNy-1)*=(fact*rinv(nNy-1));    
  for(int ix=nx-2;ix>=0;ix--) phi(ix,nNy-1)=(phi(ix,nNy-1)*fact+phi(ix+1,nNy-1))*rinv(nNy-1);
  phi(0,nNy-1)=phi(0,nNy-1)*psol(nNy-1);
  for(int ix=1;ix<nx;ix++) phi(ix,nNy-1)=phi(ix-1,nNy-1)*rinv(nNy-1)-phi(ix,nNy-1);
  
  phi(-1,nNy-1)=phi(0,nNy-1)*rinv(nNy-1);
  phi(nx,nNy-1)=phi(nx-1,nNy-1)*rinv(nNy-1);
  phi(nx+1,nNy-1)=phi(nx,nNy-1)*rinv(nNy-1);
  
  //Calculation of phi(x,y). This inverse trnsform has not been normalized by ny. This is 
  //taken care when phi_dc is calculated.
  rfftw(plan_bckwd,nx+3,(fftw_real*) &(phi(-1,0)),1,ny,(fftw_real*) &(phi(-1,0)),1,ny);
    
  //phi dc calculation
  rho_dc=0.;
  for(int ix=-1;ix<nx+2;ix++) for(int iy=0;iy<ny;iy++) rho_dc+=phi(ix,iy);
  rho_dc/=((nx+3)*ny);

  //The division of ny is for normalization of the IFFT.
  for(int ix=-1;ix<nx+2;ix++) for(int iy=0;iy<ny;iy++) phi(ix,iy)/=ny;//phi(ix,iy)=(phi(ix,iy)-rho_dc)/ny;
 
 /*cout<<"phi at the end of efield_injection:"<<endl; 
 for(int ix=-1;ix<nx+2;ix++){
   cout<<endl;
   for(int iy=0;iy<ny;iy++) cout<<phi(ix,iy)<<" ";
 }*/

}
