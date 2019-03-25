
#include "eppic-mpi.h"
#include "tridiag.h"
#include "math.h"
#include "eppic-efield.h"
#define PI  3.1415926535897932385

int main(int argc, char* argv[])
{

  
  int mpi_np,mpi_rank;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_np);

#ifdef DEBUG
  int nx=128/mpi_np;
  int nit = 1;
#else
  if (argc<3) {
    MPI_Finalize();
    return 1;
  }
  int nx=atoi(argv[1])/mpi_np;
  int nit = atoi(argv[2]);
#endif  
  FTYPE xmax=1.0,xmin=0;
  FTYPEAVec a=FTYPEAVec(nx)=1;
  FTYPE new_b_const = 2.58578643762691;
  FTYPEAVec b=FTYPEAVec(nx)=-new_b_const;
  FTYPEAVec c=FTYPEAVec(nx)=1;
  FTYPEAVec rho=FTYPEAVec(nx)=0;
  FTYPEAVec rhs=FTYPEAVec(nx)=0;
  FTYPEAVec lhs=FTYPEAVec(nx)=0;


  FTYPE delx = (xmax-xmin)/FTYPE(nx*mpi_np/2);
  FTYPE eps = 1;
  FTYPE rho_dc = 0;

  for (int i=0;i<nx;i++) {
    rho[i] = cos(2*PI*3*(i+mpi_rank*nx)/(nx*mpi_np))+10;//-9.72*eps/delx/delx ;
    //    FTYPE x=(xmin+i*delx)+mpi_rank*nx*delx;
    //sin(3.14*x*2);//(x-(xmin+xmax)/2)*(x-(xmin+xmax)/2);//cos(2*PI*3*(i+mpi_rank*nx)/(nx*mpi_np))*
    //(x-10)*(x-10)*(x-10);//(x-(xmax+xmax)/4.0);
    //    rho[i] *=rho[i];
    rho_dc += rho[i];
  }

  FTYPE rho_dc_all=0;
  MPI_Allreduce( &rho_dc, &rho_dc_all, 1,
		 MPI_FTYPE, MPI_SUM, MPI_COMM_WORLD);

  rho_dc =rho_dc_all/(nx*mpi_np);


  for (int i=0;i<nx;i++) {
    rho[i] -= rho_dc;
    rhs[i] = -rho[i]*delx*delx/eps;
  }

  /* boundary conditions */
  // bc_type LHS_boundary = neumann;
//   bc_type RHS_boundary = open;
//   FTYPE LHS_dirichlet = -.41;
//   FTYPE RHS_dirichlet = -.4;

//  FTYPE rhs_open_all = 0;




  if (mpi_rank==0) {
    //b[0] = 1;
    //rhs [0] = 0;
    //c[0] = 0;
    //a[0] = 0;
  }
  if (mpi_rank==mpi_np-1) {
    //c[nx-1] = 0;
  }
  

  FTYPEAVec lhs_all=FTYPEAVec(nx*mpi_np)=0;
  FTYPEAVec rho_all=FTYPEAVec(nx*mpi_np)=0;
  FTYPEAVec rho_calc_all=FTYPEAVec(nx*mpi_np)=0;

  FTYPE time=0;
  for (int it=0; it<nit; it++) {
    if (mpi_np>1) {
      
      FTYPE time_before=MPI_Wtime();
      tridiag_cyclic_solver_parallel(a,b,c,rhs,lhs,1,1,MPI_COMM_WORLD);
      //tridiag_solver_parallel(a,b,c,rhs,lhs,MPI_COMM_WORLD);
      time += MPI_Wtime()-time_before;

#ifdef DEBUG
      MPI_Gather(lhs.address(0),nx,MPI_FTYPE,
		 lhs_all.address(0),nx,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Gather(rho.address(0),nx,MPI_FTYPE,
		 rho_all.address(0),nx,MPI_FTYPE,0,MPI_COMM_WORLD);
#endif
    } else {
      tridiag_cyclic_solver(a,b,c,rhs,lhs_all,1,1);
      //tridiag_solver(a,b,c,rhs,lhs_all);
      rho_all = rho;
    }
    
#ifdef DEBUG
    if (mpi_rank==0) {
      /* average phi to zero */
      /*
	FTYPE phi_sum = lhs_all.sum()/(mpi_np*nx);
	for (int ix=0;ix<mpi_np*nx;ix++)
	lhs_all[ix]-=phi_sum;
      */
      /* calculate rho */
      /** first point **/
      rho_calc_all[0] = -(lhs_all[1]-new_b_const*lhs_all[0]+lhs_all[nx*mpi_np-1])
	*eps/(delx*delx);
      /** other points **/
      for (int ix=1;ix<nx*mpi_np-1;ix++) {
	//      rho_calc_all
	rho_calc_all[ix] = -(lhs_all[ix-1]-new_b_const*lhs_all[ix]+lhs_all[ix+1])*
	  eps/(delx*delx);
      }
      /** last point **/
      int ix=nx*mpi_np-1;
      rho_calc_all[ix] = -(lhs_all[ix-1]-new_b_const*lhs_all[ix]+lhs_all[0])
	*eps/(delx*delx);

      /* report */
      char buffer[80];
      sprintf(buffer,"tridiag_cyclic_solver_parallel_test_%d.dat\0",mpi_np);
      std::ofstream gnuplot_out(buffer);
      for (int i=0;i<nx*mpi_np;i++) {
      
	sprintf(buffer,"%8f\t%8d\t%16.8e\t%16.8e\t%16.8e\n",
		i*delx,
		i,
		lhs_all(i),
		rho_calc_all(i),
		rho_all(i));
	gnuplot_out << buffer;
      }
      gnuplot_out.close();
    }
#endif

  }
#ifndef DEBUG
  if (mpi_rank==0) {
    char buffer[80];
    sprintf(buffer,"%4d %16d %4d %16.8e %16.8e\n",
	    mpi_np,nx,nit,time,time/nit);
    cout << buffer;
  }
#endif
  MPI_Finalize();
}
