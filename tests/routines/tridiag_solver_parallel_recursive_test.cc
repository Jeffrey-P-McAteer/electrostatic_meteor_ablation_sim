
#include "eppic-mpi.h"
#include "tridiag.h"
#include "math.h"
#include "eppic-efield.h"

int main(int argc, char* argv[])
{

  
  int mpi_np,mpi_rank;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_np);

  int max_np=4;

#ifdef DEBUG
  int nx=100/mpi_np;
  int nit = 1;
#else
  if (argc<3) {
    MPI_Finalize();
    return 1;
  }
  int nx=atoi(argv[1])/mpi_np;
  int nit = atoi(argv[2]);
  if (argc==4) 
    max_np = atoi(argv[3]);
#endif  
  FTYPE xmax=1.0,xmin=0;
  FTYPEAVec a=FTYPEAVec(nx)=1;
  FTYPEAVec b=FTYPEAVec(nx)=-2;
  FTYPEAVec c=FTYPEAVec(nx)=1;
  FTYPEAVec rho=FTYPEAVec(nx)=0;
  FTYPEAVec rhs=FTYPEAVec(nx)=0;
  FTYPEAVec lhs=FTYPEAVec(nx)=0;


  FTYPE delx = (xmax-xmin)/FTYPE(nx*mpi_np-1);
  FTYPE eps = 1;
  FTYPE rho_dc = 0;

  for (int i=0;i<nx;i++) {
    FTYPE x=(xmin+i*delx)+mpi_rank*nx*delx;
    rho[i] = (x-(xmin+xmax)/2)*(x-(xmin+xmax)/2);
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
  bc_type LHS_boundary = neumann;
  bc_type RHS_boundary = open;
  FTYPE LHS_dirichlet = -.41;
  FTYPE RHS_dirichlet = -.4;

  FTYPE rhs_open_all = 0;
  if (RHS_boundary == open) {
    /* calculate lhs boundary */
    FTYPE rho_open = 0;
    for (int ix=0; ix<nx; ix++)
      rho_open += (nx*mpi_np-ix-nx*mpi_rank)*rho[ix];

    MPI_Allreduce(&rho_open,&rhs_open_all,1,MPI_FTYPE,MPI_SUM,MPI_COMM_WORLD);
  }
  

  if (mpi_rank==mpi_np-1) {
    switch(RHS_boundary) {
    case neumann:
      /*
      cout << "Neuman RHS " << endl;
      */
      /* neumann */
      c[nx-1] = 0;
      a[nx-1] = 1;
      b[nx-1] = -1;
      /*
      if (LHS_boundary == RHS_boundary) {
	rhs[nx-2]=0;
	b[nx-1]=1;
	a[nx-1]=0;
	c[nx-1]=0;
	  }

      */
      break;
    case dirichlet:
      /*
      cout << "Dirichlet RHS " << endl;
      */
      /* dirichlet */
      rhs[nx-1] -=RHS_dirichlet;
      break;
    case open:
      /*
	cout << "Open RHS " << endl;
	cout << "open boundary = " << rhs_open_all << endl;
      */
      /* dirichlet */
      rhs[nx-1] -=rhs_open_all;
      break;
    }
  }

  FTYPE lhs_open_all=0;
  if (LHS_boundary == open) {
    /* calculate lhs boundary */
    FTYPE rho_open = 0;
    for (int ix=0; ix<nx; ix++)
      rho_open += (ix+nx*mpi_rank+1)*rho[ix];
    MPI_Reduce(&rho_open,&lhs_open_all,1,MPI_FTYPE,MPI_SUM,0,MPI_COMM_WORLD);
  }
  


  if (mpi_rank == 0) {
    switch(LHS_boundary) {
    case neumann:
      /*
      cout << "Neuman LHS " << endl;
      */
      /* neumann */
      a[0] = 0;
      c[0] = 1;
      b[0] = -1;
      break;
    case dirichlet:
      /*
      cout << "Dirichlet LHS " << endl;
      */
      /* dirichlet */
      rhs[0] -=LHS_dirichlet;
      break;
    case open:
      /*
      cout << "Open LHS " << endl;
      cout << "open boundary = " << lhs_open_all << endl;
      */
      /* dirichlet */
      rhs[0] -=lhs_open_all;
      break;
    }
  }
  
  FTYPEAVec lhs_all=FTYPEAVec(nx*mpi_np)=0;
  FTYPEAVec rho_all=FTYPEAVec(nx*mpi_np)=0;
  FTYPEAVec rho_calc_all=FTYPEAVec(nx*mpi_np)=0;

  FTYPE time=0;
  for (int it=0; it<nit; it++) {
    if (mpi_np>1) {
      
      FTYPE time_before=MPI_Wtime();
      tridiag_solver_parallel_recursive(a,b,c,rhs,lhs,MPI_COMM_WORLD,max_np);
      time += MPI_Wtime()-time_before;

#ifdef DEBUG
      MPI_Gather(lhs.address(0),nx,MPI_FTYPE,
		 lhs_all.address(0),nx,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Gather(rho.address(0),nx,MPI_FTYPE,
		 rho_all.address(0),nx,MPI_FTYPE,0,MPI_COMM_WORLD);
#endif
    } else {
      tridiag_solver(a,b,c,rhs,lhs_all);
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
      switch(LHS_boundary) {
      case neumann:
	rho_calc_all[0] = -(lhs_all[1]-1*lhs_all[0])*eps/(delx*delx);
	break;
      case dirichlet:
	rho_calc_all[0] = -(lhs_all[1]-2*lhs_all[0]+LHS_dirichlet)
	  *eps/(delx*delx);
	break;
      case open:
	rho_calc_all[0] = -(lhs_open_all-2*lhs_all[0]+lhs_all[1])
	  *eps/(delx*delx);
	break;
      }
      /** other points **/
      for (int ix=1;ix<nx*mpi_np-1;ix++) {
	//      rho_calc_all
	rho_calc_all[ix] = -(lhs_all[ix-1]-2*lhs_all[ix]+lhs_all[ix+1])*
	  eps/(delx*delx);
      }
      /** last point **/
      int ix=nx*mpi_np-1;
      switch(RHS_boundary) {
      case neumann:
	rho_calc_all[ix] = -(lhs_all[ix-1]-1*lhs_all[ix])*eps/(delx*delx);
	break;
      case dirichlet:
	rho_calc_all[ix] = -(lhs_all[ix-1]+RHS_dirichlet-2*lhs_all[ix])
	  *eps/(delx*delx);
	break;
      case open:
	rho_calc_all[ix] = -(lhs_all[ix-1]-2*lhs_all[ix]+rhs_open_all)
	  *eps/(delx*delx);
	break;
      }
      cout << "rhs_boundary = " << rho_calc_all[nx*mpi_np-1] << endl;

      /* report */
      char buffer[80];
      sprintf(buffer,"tridiag_solver_parallel_recursive_test_%d.dat\0",mpi_np);
      std::ofstream gnuplot_out(buffer);
      for (int i=0;i<nx*mpi_np;i++) {
      
	sprintf(buffer,"%8f\t%8d\t%16.8e\t%16.8e\t%16.8e\n",i*delx,i,
		lhs_all(i),
		rho_calc_all(i),
		rho_all(i));
	gnuplot_out << buffer;
      }
      gnuplot_out.close();

      /* print out what boundaries should be for this rho */
      /* rhs */
      cout << "RHS should be ";
      cout << -rho_all[nx*mpi_np-1]/eps*delx*delx+2*lhs_all[nx*mpi_np-1]-
	lhs_all[nx*mpi_np-2] << endl;
      cout << "LHS should be ";
      cout << -rho_all[0]/eps*delx*delx+2*lhs_all[0]-
	lhs_all[1] << endl;
    }
#endif

  }
#ifndef DEBUG
  if (mpi_rank==0) {
    char buffer[80];
    sprintf(buffer,"%4d %4d %16d %4d %16.8e %16.8e\n",
	    mpi_np,max_np,nx,nit,time,time/nit);
    cout << buffer;
  }
#endif
  MPI_Finalize();
}
