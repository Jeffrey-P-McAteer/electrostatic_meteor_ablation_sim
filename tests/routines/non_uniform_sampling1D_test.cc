#include "eppic-pic.h"
#include "eppic-mpi.h"
#include "eppic-math.h"
#include "eppic-misc.h"

ArrayNd<MPI_Comm,1> comm_across_subdomains;
ArrayNd<MPI_Comm,1> comm_within_subdomain;
int proc_west;
int proc_east;
int boundary_type;
int mpi_rank=0,mpi_np=1;
int nsubdomains=1;
subdomain_traits subdomain;

#if NDIM == 2

int main(int argc, char* argv[])
{
  nsubdomains=1;

  // setup MPI
  MPI_Errhandler mpi_err_hand;
  void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
  if (mpi_rank == 0) {
    MPI_Errhandler_get(MPI_COMM_WORLD,&mpi_err_hand);
  }
  if (mpi_np>1) nsubdomains=2;
  domain_decomp();

  particle_dist pic;
  eppic_system sys_params;
  sys_params.ndim = 2;
  sys_params.nx = 100;
  sys_params.ny = 100;
  sys_params.nz = 1;
  sys_params.dx = 1;
  sys_params.dy = 1;
  sys_params.dz = 1;
  sys_params.eps = 1;
  sys_params.boundary_type = PERIODIC;
  int extra_x = 10;
  int extra_y = 10;

  init_particle_dist(pic, 1.0, 1.0, 1.0, 
		     sys_params.nx*extra_x*sys_params.ny*extra_y
		     , 1.01, 3, sys_params);

  FTYPEAVec prob_dist;
  if (subdomain.id_number==0) {
    prob_dist= FTYPEAVec(sys_params.nx+2);
  } else
    prob_dist= FTYPEAVec(sys_params.nx+1);

  int xstart=0;
  FTYPE slope = 1e-4;
  int xshift = subdomain.id_number*sys_params.nx+xstart;
  if (subdomain.id_number==0)
    xstart=-1;
  for (int ix=xstart; ix<sys_params.nx+1; ix++) {
    prob_dist(ix-xstart) =
      1.0+exp(-slope*Sqr(FTYPE(1+ix+xshift)));
  }

  /*
  for (int ix=0;ix<sys_params.nx;ix++)
    prob_dist(ix) = 1.0;
  */
  //periodic 
  //prob_dist(sys_params.nx) = prob_dist(0);
  FTYPE global_max=0,global_avg_max=0, rightarea=0,leftarea=0;
  prob_dist_1D_ddecomposed(prob_dist,global_max,global_avg_max,
			   rightarea,leftarea);
  
  pic.np = int(pic.np/prob_dist.size()*
		   prob_dist.sum()/global_avg_max);
  pic.n0avg*=global_avg_max/global_max;


  non_uniform_sampling1D(pic, prob_dist, sys_params.nx*extra_x,
			 sys_params, xstart, 0);





  ArrayNd<FTYPE,2> bin = ArrayNd<FTYPE,2>(sys_params.nx+1,sys_params.ny);
  void gather_domain_inject(FArrayND &den, PTYPEAVec &x, PTYPEAVec &y, const int np, 
			    int nxmin, FTYPE scale);

  gather_domain_inject(bin,pic.x,pic.y,pic.np,-1,1);
		/*
  for (int i=0;i<pic.np;i++) {
    int ix = int(floor(pic.x(i)));
    int iy = int(floor(pic.y(i)));
    if (pic.x(i)>=sys_params.nx) 
      cout << "Particle "
	   << i
	   <<" is out of range: "
	   <<pic.x(i)
	   <<endl;
    else if (pic.x(i)>sys_params.nx-1) 
      cout << "Particle "
	   << i
	   <<" is just in range: "
	   <<pic.x(i)
	   <<endl;
    bin(ix,iy)=bin(ix,iy)+1;
  }
		*/
  
  double std_dev=0;

  char buffer[80];
  std::ofstream gnuplot_out;
  sprintf(buffer,"non_uniform_sampling1D_test_%1d.dat\0",mpi_rank);
  gnuplot_out.open(buffer);
  for (int ix=0; ix<sys_params.nx;ix++){
  for (int iy=0; iy<sys_params.ny;iy++){
    sprintf(buffer,"%8d\t%8d\t%16.8e\n",ix+subdomain.id_number*sys_params.nx,iy,bin(ix,iy));
    gnuplot_out << buffer;
    std_dev+=Sqr(FTYPE(bin(ix,iy))/(extra_x*extra_y)-prob_dist(ix));
  }
  gnuplot_out << endl;
  }
  gnuplot_out.close();

  cout << "Standard Dev = " << sqrt(std_dev/(bin.size())) << endl;
  sprintf(buffer,"non_uniform_sampling1D_test_pic_%1d.dat\0",mpi_rank);
  gnuplot_out.open(buffer);
  for (int i=0;i<pic.np;i++) {
    sprintf(buffer,"%8d\t%16.8e\t%16.8e\n",i,pic.x(i)+subdomain.id_number*sys_params.nx,
	    pic.y(i));
    gnuplot_out << buffer;
  }
  gnuplot_out.close();


  MPI_Finalize();
  return 0;

}

#else

int main(int argc, char* argv[]){
  cout << "3D test not implemented!" << endl;
  return 1;

}

#endif
