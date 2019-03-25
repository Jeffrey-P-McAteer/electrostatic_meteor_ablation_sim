
#include "eppic-types.h"
#include "eppic-mpi.h"
#include "tridiag.h"



int main(int argc, char* argv[])
{

  
  int mpi_np,mpi_rank;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_np);

  int nx=40/mpi_np;
  FTYPE xmax=1.0,xmin=0;
  FTYPEAVec a=FTYPEAVec(nx)=1;
  FTYPEAVec b=FTYPEAVec(nx)=-2;
  FTYPEAVec c=FTYPEAVec(nx)=1;
  FTYPEAVec rhs=FTYPEAVec(nx)=0;
  FTYPEAVec lhs=FTYPEAVec(nx)=0;


  FTYPE delx = (xmax-xmin)/FTYPE(nx*mpi_np+1);
  FTYPEAVec solution=FTYPEAVec(nx*mpi_np+1)=0;



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

  if (mpi_rank ==0){
    FTYPE g = (gamma_l*(alpha_h+beta_h)-beta_l*(gamma_h-(alpha_h+beta_h)/3.0))/
      (alpha_l*alpha_h+alpha_l*beta_h-beta_l*alpha_h);
    FTYPE h = (alpha_l*(gamma_h-(alpha_h+beta_h)/3.0)-gamma_l*alpha_h)/
      (alpha_l*alpha_h+alpha_l*beta_h-beta_l*alpha_h);
    
    for (int i=1;i<nx*mpi_np+1;i++) {
      FTYPE x=(xmin+i*delx);
      solution[i] = g+h*x+(x*x)/2.0-((x*x*x*x)/6.0);
    }
    solution[0] = g+h*xmin+(xmin*xmin)/2.0-((xmin*xmin*xmin*xmin)/6.0);
  }

  for (int i=1;i<nx+1;i++) {
    FTYPE x=(xmin+i*delx)+mpi_rank*nx*delx;
    rhs[i-1] = (1-2*x*x)*delx*delx;
  }

  //  solution[nx+1] = g+h*xmax+(xmax*xmax)/2.0-((xmax*xmax*xmax*xmax)/6.0);

  /* Adjust matrix values for boundary conditions */

  if (mpi_rank==mpi_np-1) {
    c[nx-1] = 0;
    b[nx-1] += beta_h/
      (alpha_h*delx+beta_h);
    rhs[nx-1] -= gamma_h*delx/
      (alpha_h*delx+beta_h);
  }

  if (mpi_rank == 0) {
    a[0] = 0;
    b[0] -= beta_l/
      (alpha_l*delx-beta_l);
    rhs[0] -= gamma_l*delx/
      (alpha_l*delx-beta_l);
  }
  
  FTYPEAVec lhs_all=FTYPEAVec(nx*mpi_np)=0;

  if (mpi_np>1) {
    tridiag_solver_parallel(a,b,c,rhs,lhs,MPI_COMM_WORLD);
    MPI_Gather(lhs.address(0),nx,MPI_FTYPE,
	       lhs_all.address(0),nx,MPI_FTYPE,0,MPI_COMM_WORLD);
  } else 
    tridiag_solver(a,b,c,rhs,lhs_all);

  if (mpi_rank==0) {
    char buffer[80];
    sprintf(buffer,"tridiag_solver_parallel_test_%d.dat\0",mpi_np);
    std::ofstream gnuplot_out(buffer);
    for (int i=1;i<nx*mpi_np+1;i++) {
      sprintf(buffer,"%8f\t%16.8e\t%16.8e\n",i*delx,lhs_all(i-1),solution(i));
      gnuplot_out << buffer;
    }
    gnuplot_out.close();
  }

  MPI_Finalize();
}
