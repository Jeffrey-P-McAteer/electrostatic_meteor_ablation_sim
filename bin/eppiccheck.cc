// a small program to check the input for an eppic run

#include <cstdio>
#include <iostream>
#include <ctime>
//#include <locale.h>
//#include <signal.h>

#include "eppic.h"
#include "eppic-mpi.h"

// Parameters which must be known globally: 
// (see eppic.h for definitions) 
char title[132]="\0"; 
FTYPE dt=1;
FTYPE dx=1;
FTYPE dy=dx;
FTYPE dz=dx;
FTYPE eps = 1;
FTYPE Bx = 0, Bz = 0;
FTYPE damp_nu = 0;
int damp_nx = -1;
FTYPE Ermag=0;
FTYPE Ey0_external=0.; 
int nEext=0;
int ndim_space=2;
int nx=64;
int ny=1;
int nz=1;
int nt=1000;
int it=0;
#ifdef DEBUG_LOGGING
int debug_log_level=10,debug_log_time=0,debug_tab=0;
time_t debug_time;
FILE * DEBUG_LOG_FILE;
#endif
int nout=16;
int nout_avg=1;
int npout=1;
int iwrite=25;
int iread=0;
int divj_out_subcycle=-1; 
int charge_out_subcycle=-1; 
int phi_out_subcycle=1; 
int E_out_subcycle=-1; 
char outdir[256];      // Directory for output (entered through argv[1])

int iontype=0;
int ndist=0;
FTYPEAVec n0peak, vxthd, vythd, vzthd, vx0d, vy0d, vz0d;
FTYPEAVec qd, md;
FTYPEAVec species_Bx, species_Bz;
FTYPEAVec vxrd, vyrd, vzrd;
FTYPEAVec n0b, vxthb, vythb, vzthb, vx0b, vy0b, vz0b; 
FTYPE fwidth=0,fsteep=0;
int local=TRUE;
int no_parallel_efields=FALSE;
  int efield_algorithm=0;
  int efield_petsc_algorithm=1;
int kill_modes_ang=0;
int density_addition=TRUE;
intAVec npd;
intAVec chargeon;
intAVec init_dist;
FTYPEAVec part_pad; //This is the amount of extra space to be allocated 
FTYPEAVec coll_rate;
FTYPE v0_neutral[3]={0.,0.,0.}; 
FTYPEAVec vx0d_neutral,vy0d_neutral,vz0d_neutral;
FTYPE vth_neutral=1.0, m_neutral=1.0;
FTYPEAVec vxthd_neutral, vythd_neutral, vzthd_neutral, massd_neutral;
FTYPEAVec vrelmax;

intAVec method;
extern int boundary_type=PERIODIC; 

FTYPE w0=0.;
FTYPEAVec gvxmin, gvxmax;
FTYPEAVec gvymin, gvymax;
FTYPEAVec gvzmin, gvzmax;
intAVec pnvx, pnvy, pnvz;    
FTYPEAVec pvxmin, pvxmax, pvymin, pvymax, pvzmin, pvzmax;
FTYPEAVec pdamp_nu;
FTYPEAVec param1, param2, param3, param4, param5, param6,param7,param8,param9,
  param10;
intAVec subcycle;
intAVec supercycle;
intAVec den_out_subcycle, part_out_subcycle, vdist_out_subcycle, 
  flux_out_subcycle, nvsqr_out_subcycle, phistore_out_subcycle;
intAVec species_dim;
intAVec vel_dim;
FTYPEAVec thermal_gamma, diffc;

int fndist=0;
ArrayNd<FTYPE,2> fnmax, fnmin; // Range to be evaluated 
ArrayNd<int,2> fn;             // Mesh size to be used
intAVec fnout; 
intAVec fnof;  
intAVec fof[MAXDIST]; 

//Parallel processing variables will be known globally:
int mpi_rank=0, mpi_np=1;

//dls/////////////////////////////////////
// parallel
int nsubdomains=1;
#ifdef USE_MPI
ArrayNd<MPI_Comm,1> comm_within_subdomain;
ArrayNd<MPI_Comm,1> comm_across_subdomains;

subdomain_traits subdomain;

int iabsent = -1;   // ending index of absent array
int proc_east = -1;
int proc_west = -1;
#endif
//dls end////////////////////////////////

FTYPE *workarray; // A workspace array will be declared once the array sizes are known.
int nworkarray;

tms times_buf;
clock_t 
  vadvance_time=0,
  xadvance_time=0,
  charge_time=0,    
  collect_time=0,   
  efield_time=0,
  output_time=0,
  fluid_time=0;

//This will be set to true if the proper signal (currently, SIGUSR2) is sent
int end_asap=FALSE;

// Set the math exception handeler to stop on floating point errors 
//typedef void (*sighandler_t)(int);
//sighandler_t signal(int signum, sighandler_t handler);
//void fpe_exception_handler(int){raise(SIGTRAP);terminate(-1,"fpe_error");}

main(int argc, char* argv[])
{
  particle_dist *pic;
  fluid *fspecie;
  field Efield;
  particle_misc *misc;
  FArrayND rho;
  int it0=1;

  extern void infile(char *,int);
  extern void check(void);

  // Detect an end of job signal passed by batch processors:
  //  void finish_asap(int sig);
  //  signal(SIGUSR2, finish_asap);

  // Start the parallel processor 
#ifdef USE_MPI
  int mpi_err;
  MPI_Errhandler mpi_err_hand;
  void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg);
  
  MPI_Init(&argc,&argv);
  mpi_err=MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  if ( mpi_err != MPI_SUCCESS) 
    mpi_error(mpi_rank, mpi_err,"rank call in main");
  mpi_err=MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
  if ( mpi_err != MPI_SUCCESS) 
    mpi_error(mpi_rank, mpi_err,"size call in main");
  if (mpi_rank == 0) {
    mpi_err=MPI_Errhandler_get(MPI_COMM_WORLD,&mpi_err_hand);
    printf("\nMPI STARTING %d PROCESSORS\n MPI Error Handler: %d\n\n", 
	   mpi_np, mpi_err_hand);
  }
#endif


  // Syncronize the c++ and c print buffers:
  ios::sync_with_stdio();
      
  //  Read parameters from the input file: 
  infile(argv[1],0);
  
  if ((fsteep != 0 ) && (fwidth != 0)) 
    printf(";WARNING: Electric Field artificial damping in effect\n");

  

  mpi_np = nsubdomains;

  // Check the input parameters to make sure no obvious errors were made 
  check();
  mpi_np = 1;
  terminate(end_asap,"\n");

}
