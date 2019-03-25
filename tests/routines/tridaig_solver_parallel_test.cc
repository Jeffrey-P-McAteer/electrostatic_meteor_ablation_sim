
#include "eppic-types.h"
#include "tridiag.h"



int main(int argc, char* argv[])
{

  int nx=100;
  FTYPE xmax=1.0,xmin=0;
  FTYPE delx = (xmax-xmin)/FTYPE(nx+1);
  FTYPEAVec a=FTYPEAVec(nx)=1;
  FTYPEAVec b=FTYPEAVec(nx)=-2;
  FTYPEAVec c=FTYPEAVec(nx)=1;
  FTYPEAVec rhs=FTYPEAVec(nx)=0;
  FTYPEAVec lhs=FTYPEAVec(nx)=0;
  FTYPEAVec solution=FTYPEAVec(nx+1)=0;

  /* 
     model problem taking from Comp. Phys. course found online:
     http://farside.ph.utexas.edu/teaching/329/lectures/node66.html
     It's a mixed boundary condition example
  */

  /* parameters defining boundary conditions */
  /*
    // just another set of boundary conditions
  FTYPE alpha_l = 1;
  FTYPE beta_l =  0;
  FTYPE gamma_l = 0;
  FTYPE alpha_h = 1;
  FTYPE beta_h =  0;
  FTYPE gamma_h = 0;
  */

  FTYPE alpha_l = 1;
  FTYPE beta_l = -1;
  FTYPE gamma_l = 1;
  FTYPE alpha_h = 1;
  FTYPE beta_h =  1;
  FTYPE gamma_h = 1;

  /* parameters used in analytic solution */

  FTYPE g = (gamma_l*(alpha_h+beta_h)-beta_l*(gamma_h-(alpha_h+beta_h)/3.0))/
    (alpha_l*alpha_h+alpha_l*beta_h-beta_l*alpha_h);
  FTYPE h = (alpha_l*(gamma_h-(alpha_h+beta_h)/3.0)-gamma_l*alpha_h)/
    (alpha_l*alpha_h+alpha_l*beta_h-beta_l*alpha_h);

  for (int i=1;i<nx+1;i++) {
    FTYPE x=(xmin+i*delx);
    solution[i] = g+h*x+(x*x)/2.0-((x*x*x*x)/6.0);
    rhs[i-1] = (1-2*x*x)*delx*delx;
  }
  solution[0] = g+h*xmin+(xmin*xmin)/2.0-((xmin*xmin*xmin*xmin)/6.0);
  //  solution[nx+1] = g+h*xmax+(xmax*xmax)/2.0-((xmax*xmax*xmax*xmax)/6.0);

  /* Adjust matrix values for boundary conditions */

  c[nx-1] = 0;
  a[0] = 0;

  b[0] -= beta_l/
    (alpha_l*delx-beta_l);
  b[nx-1] += beta_h/
    (alpha_h*delx+beta_h);

  rhs[0] -= gamma_l*delx/
    (alpha_l*delx-beta_l);
  rhs[nx-1] -= gamma_h*delx/
    (alpha_h*delx+beta_h);
  
  tridiag_solver(a,b,c,rhs,lhs);


  char buffer[80];
  std::ofstream gnuplot_out("tridiag_solver_test.dat");
  for (int i=1;i<nx+1;i++) {
    sprintf(buffer,"%8f\t%16.8e\t%16.8e\n",i*delx,lhs(i-1),solution(i));
    gnuplot_out << buffer;
  }
  gnuplot_out.close();


}
