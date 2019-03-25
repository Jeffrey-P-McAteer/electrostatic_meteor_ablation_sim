#include "eppic-pic.h"
#include "eppic-mpi.h"
#include "eppic-math.h"

ArrayNd<MPI_Comm,1> comm_across_subdomains;
ArrayNd<MPI_Comm,1> comm_within_subdomain;
int proc_west;
int proc_east;
int boundary_type;
int mpi_rank=0,mpi_np=1;
int nsubdomains=1;
subdomain_traits subdomain;


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
  int extra_x = 3;
  int extra_y = 3;

  init_particle_dist(pic, 1.0, 1.0, 1.0, 
		     sys_params.nx*extra_x*sys_params.ny*extra_y
		     , 1.01, 3, sys_params);

  FTYPEAVec prob_dist = FTYPEAVec(sys_params.nx+1);
  for (int ix=0; ix<sys_params.nx/2; ix++) {
    prob_dist(ix) = ix;
  }
  for (int ix=sys_params.nx/2; ix<sys_params.nx; ix++) {
    prob_dist(ix) = sys_params.nx/2 - (ix-sys_params.nx/2);
  }

  /*
  for (int ix=0;ix<sys_params.nx;ix++)
    prob_dist(ix) = 1.0;
  */
  //periodic 
  prob_dist(sys_params.nx) = prob_dist(0);

  non_uniform_sampling1D(pic, prob_dist, sys_params.nx*extra_x,
			 sys_params, 0, 0);

  maxwel_vel_dist_loading(pic,3,100.0,100.0,100.0,0);

  ArrayNd<int,2> bin = ArrayNd<int,2>(sys_params.nx,sys_params.ny);
  ArrayNd<int,1> binvx = ArrayNd<int,1>(100);
  ArrayNd<int,1> binvy = ArrayNd<int,1>(100);
  ArrayNd<int,1> binvz = ArrayNd<int,1>(100);
  FTYPE dv = 0.08;
  int pseed=prime(5);
  for (long long int i=0;i<pic.np;i++) {
    int ix = int(floor(pic.x(i)));
    int iy = int(floor(pic.y(i)));
    int vx = int(pic.vx(i)/100/dv+50);if (vx<0) vx=0; if (vx>99) vx=99;
    int vy = int(pic.vy(i)/100/dv+50);if (vy<0) vy=0; if (vy>99) vy=99;
    int vz = int(pic.vz(i)/100/dv+50);if (vz<0) vz=0; if (vz>99) vz=99;
    bin(ix,iy)=bin(ix,iy)+1;
    binvx(vx)=binvx(vx)+1;
    binvy(vy)=binvy(vy)+1;
    binvz(vz)=binvz(vz)+1;
  }

  
  double std_dev=0;

  char buffer[80];
  std::ofstream gnuplot_out;
  sprintf(buffer,"maxwel_vel_dist_loading_test_%1d.dat\0",mpi_rank);
  gnuplot_out.open(buffer);
  for (int ix=0; ix<sys_params.nx;ix++){
  for (int iy=0; iy<sys_params.ny;iy++){
    sprintf(buffer,"%8d\t%8d\t%8d\n",ix,iy,bin(ix,iy));
    gnuplot_out << buffer;
    std_dev+=Sqr(FTYPE(bin(ix,iy))/(extra_x*extra_y)-prob_dist(ix));
  }
  gnuplot_out << endl;
  }
  gnuplot_out.close();

  sprintf(buffer,"maxwel_vel_dist_loading_test_vel_%1d.dat\0",mpi_rank);
  gnuplot_out.open(buffer);
  for (int vx=0; vx<100;vx++) {
    sprintf(buffer,"%8d\t%16.8e\t%16.8e\t%16.8e\n",vx,
	    FTYPE(binvx(vx))/FTYPE(binvx(50)),
	    FTYPE(binvy(vx))/FTYPE(binvy(50)),
	    FTYPE(binvz(vx))/FTYPE(binvz(50)));

    gnuplot_out << buffer;
  }
  gnuplot_out.close();

  //  cout << "Standard Dev = " << sqrt(std_dev/(bin.size())) << endl;
  sprintf(buffer,"maxwel_vel_dist_loading_test_pic_%1d.dat\0",mpi_rank);
  gnuplot_out.open(buffer);
  for (int i=0;i<pic.np;i++) {
    sprintf(buffer,"%8d\t%16.8e\t%16.8e\t%16.8e\t%16.8e\t%16.8e\n",
	    i,
	    pic.x(i),pic.y(i),
	    pic.vx(i),pic.vy(i),pic.vz(i)
	    );
    gnuplot_out << buffer;
  }
  gnuplot_out.close();


  MPI_Finalize();
  return 0;

}
