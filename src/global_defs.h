#ifndef GLOBALDEFS
#define GLOBALDEFS

// Parameters which must be known globally: 
// (see eppic.h for definitions) 
char title[132]="\0"; 
eppic_system sys;
FTYPE dt=1;
FTYPE dx=1;
FTYPE dy=dx;
FTYPE dz=dx;
FTYPE eps = 1;
FTYPE Bx = 0, By = 0, Bz = 0;
FTYPE Bmag = 0, Binc = 0, Bdec = 0;
FTYPE damp_nu = 0;
int damp_nx = -1;
FTYPE Ermag=0;
FTYPE Ex0_external=0.; 
FTYPE Ex0_rate=0.; 
FTYPE Ey0_external=0.; 
FTYPE Ey0_rate=0.; 
FTYPE Ez0_external=0.;
FTYPE Ez0_rate=0.;
int nEext=0;
int ndim_space=2;
int nx=64;
int ny=1;
int nz=1;
int nt=1000;
int it=0;
int nout=16;
int full_array_nout = 128;
int nout_avg=1;
int npout=1;
int iwrite=25;
int iread=0;
int divj_out_subcycle=-1; 
int charge_out_subcycle=-1; 
int hybrid_diagnostic_subcycle=-1; 
int phi_out_subcycle=1; 
int E_out_subcycle=1; 
char outdir[256];      // Directory for output (entered through argv[1])
bool restart_nonlocal = false;
char domain_outdir[256];      // Directory for output (entered through argv[1])
int hdf_output_arrays=0; // 0 = binary output, 1 = domain-decomposed HDF, 2 = parallel HDF
int ft_output_arrays=0; //0 = regular output, 1 = Fourier transformed output

// Extra guard cells if non-periodic y or z boundaries
int yguard_size=0;
int zguard_size=0;

int iontype=0;
int ndist=0;
FTYPEAVec n0d, vxthd, vythd, vzthd, vx0d, vy0d, vz0d;
FTYPEAVec qd, md;
FTYPEAVec species_Bx, species_By, species_Bz;
FTYPEAVec vxrd, vyrd, vzrd;
FTYPEAVec n0b, vxthb, vythb, vzthb, vx0b, vy0b, vz0b; 
FTYPE fwidth=0,fsteep=0;
int local=TRUE;
int no_parallel_efields=FALSE;
int efield_algorithm=0;    // See efield_wrapper.cc for corresponding routines
int efield_zero_dc=0;
int electron_dist=0;
int dust_shape=0;
FTYPE dust_charge=0.0;
FTYPE dust_den=1.0;
FTYPE dust_sigma=0.25;
FTYPE dust_sigma_x=0.25;
FTYPE dust_sigma_y=0.25;
FTYPE dust_mid=0.50;
int stencil_size=5;
long int global_length=1;
int nullspace_method=0; // -1: None (For debugging/development)
                        //  0: Use PETSc's routines
                        //  1: Use custom routine
int check_solution=0; // 0: No output
                      // 1: Print global value of max(|Ax-b|)
                      // 2: Print running value of max(|Ax-b|)
                      // 3: Print Ax, b, & |Ax-b| for each global row
#if HAVE_PETSC
#include <petscksp.h>
Mat A;
PC pc;
KSP ksp;
MPI_Comm petsc_comm;
int petsc_np=-1; // Use all available MPI processors
#endif // HAVE_PETSC
int MAX_TRIDIAG_SOLVE=-1; // default to don't use
int kill_modes_ang=0;
int density_addition=TRUE;
intAVec npd;
intAVec chargeon;
intAVec init_dist;
FTYPEAVec part_pad; //This is the amount of extra space to be allocated 
FTYPEAVec coll_rate, coll_start_time;
intAVec coll_type;
FTYPEAVec crosssec;
intAVec crosssec_m_model;
intAVec beta_model;
intAVec e_collision_model;
FTYPE vsim_to_kmsec;
FTYPE lsqsim_to_msq;
char coulomb_type[10][30];
intAVec coulomb_subcycle;
FTYPEAVec cc_den;
int cc_rand=0;
FTYPE l_debye;
FTYPEAVec vth_gb; // Global thermal speed 
intAVec coll_create_id;
intAVec coll_create_a_id;
intAVec coll_create_b_id;
FTYPEAVec coll_create_vthb;
FTYPEAVec creation_rate, annihilation_rate, v0_shell, vth_shell, n0_shell; 
FTYPEAVec background_neutral_dens, coll_cross_section;
FTYPEAVec boost_rate, v0_boost, vth_boost;
FTYPEAVec create_v0, create_vth, create_radius, create_posx, create_posy, create_posz;
intAVec create_pair_dist; 
FTYPE v0_neutral[3]={0.,0.,0.}; 
FTYPEAVec vx0d_neutral,vy0d_neutral,vz0d_neutral;
FTYPE vth_neutral=1.0, m_neutral=1.0;
FTYPEAVec vxthd_neutral, vythd_neutral, vzthd_neutral, massd_neutral;
FTYPEAVec vrelmax;

intAVec method;
bool unscale_density=false;
bool inject_prop_dc=false; // If particle boundary type is inject, inject particles proportional to rho_dc to zero it out
int boundary_type[3]={PERIODIC,PERIODIC,PERIODIC}; 
//Added by Glenn on 2/6/2018 to allow different boundaries for different dimensions
int field_boundary_type[3][2]={{PERIODIC_FIELD,PERIODIC_FIELD},
                               {PERIODIC_FIELD,PERIODIC_FIELD},
                               {PERIODIC_FIELD,PERIODIC_FIELD}};
FTYPE field_boundary_values[3][2]={{0,0},{0,0},{0,0}};
// Commented out by Glenn on 20180413.  field_boundary_type[0][0,1] deal with this
//bc_type rhs_boundary=open;
//bc_type lhs_boundary=open;
FTYPEAVec n0lhsd,n0rhsd;
FTYPEAVec vx0lhsd,vx0rhsd,vy0lhsd,vy0rhsd,vz0lhsd,vz0rhsd;
FTYPE w0=0.;
FTYPEAVec gvxmin, gvxmax;
FTYPEAVec gvymin, gvymax;
FTYPEAVec gvzmin, gvzmax;
intAVec pnvx, pnvy, pnvz;    
FTYPEAVec pvxmin, pvxmax, pvymin, pvymax, pvzmin, pvzmax;
FTYPEAVec pdamp_nu;
FTYPEAVec param1, param2, param3, param4, param5, param6,param7,param8,param9,
  param10;
FTYPEAVec start_col, stop_col;
stringVec den_file;
intAVec subcycle;
intAVec supercycle;
intAVec den_out_subcycle, part_out_subcycle, vdist_out_subcycle, 
  flux_out_subcycle, nvsqr_out_subcycle, phistore_out_subcycle;
intAVec species_dim;
intAVec vel_dim;
FTYPEAVec denft_out_min, denft_out_kmax;
FTYPE phift_out_min,phift_out_kmax;
FTYPE Eft_out_min,Eft_out_kmax;
FTYPEAVec thermal_gamma, diffc;

int fndist=0;
ArrayNd<FTYPE,2> fnmax, fnmin; // Range to be evaluated 
ArrayNd<int,2> fn;             // Mesh size to be used
intAVec fnout; 
intAVec fnof;  
intAVec fof[MAXDIST]; 

//Parallel processing variables will be known globally:
int mpi_rank=0, mpi_np=1;

// parallel
int nsubdomains=1;
/* Developing unequal numbers of processor per domain
const int ndomain_index_max=10;
// List of domains which will have extra processors working on them:
ArrayNd<int,1> ndomain_index(ndomain_index_max,-1);
// List of multiplacation factors more processors working on them:
ArrayNd<int,1> ndomain_mult(ndomain_index_max,1);
*/

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
time_t sim_start_time;
FTYPE time_limit=0;
//This will be set to true if the proper signal (currently, SIGUSR2) is sent
int end_asap=FALSE;

// Set the math exception handeler to stop on floating point errors 
//typedef void (*sighandler_t)(int);
//sighandler_t signal(int signum, sighandler_t handler);
//void fpe_exception_handler(int){raise(SIGTRAP);terminate(-1,"fpe_error");}

hid_t H5_FID;
hid_t h5_gid;

#endif
