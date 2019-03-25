/*
Copyright (C) 1989, 1991 Meers Oppenheim Boston University, Boston, MA
02215, USA 

As a condition of using this software, the user must agree
to share any additions or modifications of the software with the
copyright holder and must obtain permission of the copyright holder
before distributing or transferring the code.  Any publications
derived from using this software must include the copyright holder as
a coauthor for a period of 3 years after the submission date of the
first paper which derived from using this software.  Beyond that time,
an acknowledgment is required.
			   
NO WARRANTY

THERE IS NO WARRANTY FOR THIS PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT
WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME
THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN
WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY
AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU
FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR
CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.

*/


// Global includes for program eppic
#ifndef EPPIC_H
#define EPPIC_H

// header automatically created and set by autotools
#ifndef SRCCONFIG_H
#define SRCCONFIG_H
#include "srcconfig.h"
#endif
#include "eppic-types.h"
#include "eppic-mpi.h"
#include "eppic-system.h"
#include <limits>

//// COMMON PRE-STUFF
#ifdef CHECK
#ifndef DEBUG
#define DEBUG
#endif
#else
#ifndef NOCHECK
#define NOCHECK
#endif
#endif
#ifndef VERBOSE
#define VERBOSE 0
#endif

#ifndef HAVE_PETSC
#define HAVE_PETSC 0
#endif

#include "eppic-fluid.h"
#include "eppic-io.h"

//// WITH DOMAIN ONLY STUFF
#ifdef USE_DOMAINS

// domain stuff for fluid
static const int denx_guard_size=2;
static const int velx_guard_size=1;



const int nfiles_parallel=8; // Write at most nfiles_parallel restart files at once.
                                 // for migrating particles.
static const int xguard_size=1;


//// WITHOUT DOMAIN STUFF
#else

// domain stuff for fluid; set to non-domain values
static const int denx_guard_size=0;
static const int velx_guard_size=0;

void output_array(FILE*, FArrayND &);


// domain stuff for field values; set to non-domain values
FTYPEAVec part_pad;
static const int xguard_size=0;

#endif//def USE_DOMAINS


//// COMMON POST-STUFF

#include "eppic-pic.h"
#include "eppic-efield.h"
#include "eppic-math.h"

// Universal parameters which do not vary 
//     qe   -real- Ion charge (coulombs) 
//     qi   -real- Ion charge (coulombs) (qe=-qi) 
//     me   -real- Electron mass (kg) 
//     mi   -real- Ion mass (kg) 
//     eps  -real- dielectric constant 
//     kb   -real- Boltzmann constant 

#if NDIM==1
static const int ndim = 1;
#elif NDIM==2
static const int ndim = 2;
#elif NDIM==3
static const int ndim = 3;
#endif 

#ifndef PI
#define PI  3.1415926535897932385
#endif
#define QI  1.6022e-19
#define QE -1.6022e-19
#define ME  9.1094e-31
#define MI  1.6726e-27
#define EPS0  8.8542e-12
#define KB  1.380658e-23
#define SQRT2 1.41421356237309504880
#define SQRT2PI 2.5066282746310007
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

const PTYPE fabsent=( std::numeric_limits<PTYPE>::max()/(-2) );
const PTYPE fabsent2=( std::numeric_limits<PTYPE>::max()/(-4) );

#define MAXDIST 64
// Parameters which must be known globalally: 
// (all parameters are in MKS units) 
// Add z and y guard sizes in case they are not periodic particle boundaries
extern int yguard_size;
extern int zguard_size;
extern char title[132]; //Run title
extern FTYPE dt;     // time step in seconds 
extern int nt;        // Number of timesteps to evaluate 
extern FTYPE dx, dy, dz; // grid step size in meters 
extern int ndim_space; // Number of spactial dimensions (default=2)
extern int nx, ny, nz;    // number of grid cells 
extern eppic_system sys;  // system variable grid
extern int it;  // current time step
extern FTYPE eps;    // dielectric constant (default=8.8542e-12) 
extern FTYPE Bx;     // x directed Magnetic field in Tesla (default=0) 
extern FTYPE By;     // y directed Magnetic field in Tesla (default=0)
extern FTYPE Bz;     // z directed Magnetic field in Tesla (default=0) 
extern FTYPE Bmag;   // Magnitude of magnetic field in Tesla (default=0)
extern FTYPE Binc;   // Inclination of magnetic field in degrees (default=0 - equatorial)
extern FTYPE Bdec;   // Declination of magnetic field in degrees (default=0 - meridional)
extern FTYPE damp_nu;// Velocity Damping coefficient (default nu=0) 
extern int damp_nx;   // Velocity Damping starting position (default dx=0) 
extern FTYPE Ermag;     // Mag. of random Efield applied to system, defined 
                         // in externalE.cc (default=0)
extern FTYPE Ex0_external; // Mag. of external Efield applied to system in X dir(default=0)
extern FTYPE Ey0_external; // Mag. of external Efield applied to system in Y dir(default=0)
extern FTYPE Ez0_external; // Mag. of external Efield applied to system in Z dir(default=0)
extern int nEext;     // number of time steps to impose the Efield defined 
                      // in externalE.cc (default=0)
extern int ndist;     // number of particle distributions 
extern FTYPEAVec n0d, vxthd, vythd, vzthd, vx0d, vy0d, vz0d, qd, md;
                      // density, thermal velocity, 
                      // drift velocity (default=0), charge (default=QE), and
                      // mass (default=ME) in distribution, all in MKS units 
extern FTYPEAVec n0b, vxthb, vythb, vzthb, vx0b, vy0b, vz0b; 
                      // density, thermal and drift velocity of second
                      // gaussian distribution 
extern FTYPEAVec species_Bx, species_By, species_Bz;
extern FTYPE fwidth,fsteep;
                      // Filtering parameters - damping goes as
                      // exp(-1*pow((i*fwidth)/(x.len/2.),2.*fsteep))
                      // where i=1,x.len/2. Default: fwidth=0, fsteep=0 
extern int local;     // If local==TRUE, calculate the fields locally
                      // otherwise, calculate the fields spectrally 
extern int kill_modes_ang; // if 90 <= kill_modes_ang < 0, then all E modes
                           // aligned closer than this angle to z-axis = 0
extern int density_addition; //density_addition==TRUE saves one array at a
                             //cost of (slightly) slower preformance.
extern bool unscale_density; // how to deal with outputting charge density
                             // if false, use original method (rho-n0avg)/n0avg
                             // if true, use rho
extern bool inject_prop_dc; // If particle boundary type is inject, inject 
                            // particles proportional to rho_dc to zero net charge
extern intAVec npd;   // number of particles in distribution 
extern intAVec chargeon; // if 1, turn the charge of a particle on only after
                         // passing through a vertical wall
extern intAVec init_dist; // Type of initial particle distribution, 
                          // 0=gaussian, 1=bit reversed 
extern FTYPEAVec part_pad; //This is the amount of extra space to be allocated 
extern FTYPEAVec coll_rate; // Collision rate for elastic scattering off a 
                            // population of neutral atoms (0 = collisionless).
extern intAVec coll_type; // Collision type (0=elastic, 1=maxwell, 2=fast_elastic, 3=nu_prop_T, 4=ionizing)
extern FTYPEAVec crosssec; // collisional cross section of a particle with the background atmosphere
extern intAVec crosssec_m_model; // model type used to calculate the cross section
extern intAVec beta_model; // model type used to calculate the ionization probability of a collusion
extern intAVec e_collision_model; // model type used to calculate the post collision electron velocity in momentum_beta_2species
extern FTYPEAVec vbth; // thermal velocity magnitude for momentum_beta_2species
extern FTYPE vsim_to_kmsec; // scaling used to convert vel simulation units to km/sec
extern FTYPE lsqsim_to_msq; // scaling used to convert length simulation units to m^2
extern FTYPEAVec coll_start_time; // time at which collisons begin
extern intAVec coll_create_id; // For ionizing, dist id of collision products
extern intAVec coll_create_a_id; // For ionizing, dist id of collision products
extern intAVec coll_create_b_id; // For ionizing, dist id of collision products
extern FTYPEAVec coll_create_vthb; // For ionizing, dist id of collision products
extern char coulomb_type[10][30]; // Coulomb collision type and pairings (0=off, 1=electron off ion)
                                  // Need a better way than hard coding the max of 10 species
extern intAVec coulomb_subcycle; // Subcyle the coulomb collision routine
extern FTYPEAVec cc_den;  // Allows scaling of the Coulomb collision frequency through field particle density

extern int cc_rand;  // Random seed for Coulomb collisions
extern FTYPE l_debye; // Global debye length calculated using the 0th distribution (needed in coulomb_scatter)

extern FTYPEAVec vth_gb; // Global thermal speed for each distribution

extern FTYPEAVec creation_rate, annihilation_rate, v0_shell, vth_shell, n0_shell; 
                           // particle creation and annihilation rates
                           // v0_shell is the speed of created particles in a shell
                           // vth_shell is the thermal width of the shell
                           // n0_shell is the density of the shell
extern FTYPEAVec background_neutral_dens, coll_cross_section;
                           // These variables are used to calculate the probability
                           // of a collision.  background_neutral_dens is the 
                           // background atmospheric density.  coll_cross_section
                           // is the collisional cross section for the species.
                           // This is used for meteor ablation simulations.

extern FTYPEAVec boost_rate, v0_boost, vth_boost;
                           // particle velocity boost rate
                           // v0_boost is the speed of boosted particles in a boost
                           // vth_boost is the thermal width of the boost
extern intAVec create_pair_dist; //Creation typically done in pairs, so this is the id of 
                                 // distribution where the opposite momentum 
                                 // and charged particles will be placed.
extern FTYPEAVec create_v0, create_vth, create_radius, create_posx, create_posy, create_posz;
                 //Creation of particles from a sphere of radius create_radius located at
                 // (create_posx, create_posy, create_posz) with average radial speed of create_v0 
                 // and a thermal speed of create_vth at a rate of creation_rate
extern FTYPE v0_neutral[3]; // Neutral drift velocity; default for all species
extern FTYPEAVec vx0d_neutral, vy0d_neutral, vz0d_neutral; // neutral drift 
                                   // velocity; distinct for all species
extern FTYPE vth_neutral, m_neutral; // temp and mass of neutrals
                            //default for all species
extern FTYPEAVec vxthd_neutral, vythd_neutral, vzthd_neutral; //temp 
                             // of neutrals; distinct for each species.
extern FTYPEAVec  massd_neutral; //mass of neutrals - distinct for each species
extern FTYPEAVec vrelmax ; //  Relative vel where particle collides
extern int dust_shape; // Shape of simple background dust (0 to suppress)
extern FTYPE dust_charge; // Static charge on simple background dust
extern FTYPE dust_den; // Relative density of simple background dust
extern FTYPE dust_sigma; // Symmetric dust e-folding length 
                         // (fraction of num. of cells)
extern FTYPE dust_sigma_x; // X-axis dust e-folding length 
                           // (fraction of num. of cells)
extern FTYPE dust_sigma_y; // Y-axis dust e-folding length 
                           // (fraction of num. of cells)
extern FTYPE dust_mid; // Dust mid point (fraction of num. of cells)
extern int stencil_size; // Number of elements in discretization 
                         // stencil when using PETSc
extern long global_length; // Length of global MPI vectors
extern int nullspace_method; // Method for handling the null space of the
                             // LHS operator
extern int check_solution; // If 1, check the PETSc solution pointwise and
                           // print the maximum absolute difference
#if HAVE_PETSC
#include <petscksp.h>
/* These PETSc objects will have global scope so that they can
   reuse information across time steps (e.g., reuse operator
   factoring for direct solvers). */
extern Mat A;   // Operator matrix
extern PC pc;   // Preconditioner
extern KSP ksp; // Krylov subspace object.
extern MPI_Comm petsc_comm; // PETSc subcommunicator
extern int petsc_np; // Number of ranks (MPI processes) on which to run PETSc
#endif // HAVE_PETSC


extern intAVec method;/* Type of Algorithm: 
			 0  = pic 
			 3  = df-calc (calculated delta-f - no time-stepping of weights, 
                         23 = df-calc & choosen marker particles - see gdist.cc
		         -1 = fluid (full dynamics - not implemented)
			 -2 = fluid (inertialess dynamics - continuity eqn) 
			 -3 = fluid (inertial dynamics - momentum eqn) 
			 -4 = fluid (parameters only - e.g. quasineutral electrons)*/

extern FTYPEAVec n0lhsd,n0rhsd; // x-boundary density conditions 
extern FTYPEAVec vx0lhsd,vx0rhsd; // x-boundary velocity conditions 
extern FTYPEAVec vy0lhsd,vy0rhsd; // x-boundary velocity conditions 
extern FTYPEAVec vz0lhsd,vz0rhsd; // x-boundary velocity conditions 
//extern bc_type rhs_boundary, lhs_boundary; // x-boundary field conditions Commented out by Glenn 4/13/2018.  field_boundary_type[0][0,1] deals with this now.


extern FTYPE w0;     // Initial magnitude of w - w will be randomly chosen 
                      // between -w0 and w0 
//IO Variables:
extern int nout;      // number of timesteps between output 
extern int full_array_nout; // number of timesteps between full array output when using FT output
extern int nout_avg;  // number of grid points to avg.
extern int npout;     // fraction of particles to print 
extern int iwrite;    // timesteps between restart dumps 
extern int iread;     // restart from dump? (0=no, 1=yes, 2=yes w/ new files)
extern int divj_out_subcycle; // Output divj only every nth output cycle
extern int charge_out_subcycle; // Output rho only every nth output cycle
extern int phi_out_subcycle; // Output phi only every nth output cycle
extern int E_out_subcycle; // Output divj only every nth output cycle
extern int hybrid_diagnostic_subcycle; // Output hybrid diagnostic arrays only every nth output cycle
extern char outdir[256];      // Directory for output (entered through argv[1])
extern bool restart_nonlocal; // true if restart dir is to be made in outdir
extern char domain_outdir[256]; // Directory for output (entered through argv[1])
extern int hdf_output_arrays; // 0 = binary, 1 = domain-decomposed HDF, 2 = collective HDF
extern int ft_output_arrays; //0 = regular output, 1 = Fourier transformed output

extern int iontype;   // ion type (0=frozen, 1=user defined PIC) 
extern intAVec pnvx, pnvy, pnvz;    
                      // velocity resolution (default=nx) for each dist.
extern FTYPEAVec pvxmin, pvxmax, pvymin, pvymax, pvzmin, pvzmax;
extern FTYPEAVec pdamp_nu; // Dist. specific vel. damping coeff. 
extern FTYPEAVec param1, param2, param3, param4, param5,param6,param7,param8,
  param9,param10;
extern FTYPEAVec start_col, stop_col;
extern stringVec den_file;
                      // Parameters for each distribution 
extern FTYPEAVec thermal_gamma, diffc;
                      // Fluid gamma (1=isothermal, 5/3=adiabatic) 
                      // Fluid Hyperdiffusion coefficient 
                      // (Note can make simulator unstable if too large - try 0.01);
extern intAVec species_dim; // Dimension of the particular species 
                            // (allows a species to have 1D behavior in a 2D or 3D run)
extern intAVec vel_dim; // Dimension of the velocity array 
extern intAVec subcycle;
                      // update this particle distribution every subcycle timesteps. 
extern intAVec supercycle;
                      // update this fluid distribution supercycle times per timestep.

// Distribution function IO:
extern intAVec den_out_subcycle, part_out_subcycle, vdist_out_subcycle,
  flux_out_subcycle, nvsqr_out_subcycle, phistore_out_subcycle;
                     //Output only every nth IO cycle (default=1)

// Relative minimum amplitude of the fourier transformed densities to output
extern FTYPEAVec denft_out_min;
extern FTYPE phift_out_min;
extern FTYPE Eft_out_min;

// Maximum abs(k) of the fourier transformed densities to output (in physical units)
extern  FTYPEAVec denft_out_kmax;
extern  FTYPE phift_out_kmax;
extern  FTYPE Eft_out_kmax;

extern int fndist;    // Declares the number of f(x,y,vx,vy) distributions 
                      // to be printed 
extern ArrayNd<FTYPE,2> fnmax, fnmin;   // Range to be evaluated 
extern ArrayNd<int,2> fn;             // Mesh size to be used

                      // number of cells in each direction 
extern intAVec fnout; // Fraction of outputs to ouput distribution 
extern intAVec fnof;  // number of particle distributions to sum 
extern intAVec fof[MAXDIST]; // array containing id numbers of distributions to sum 


// Routines which we'll declare globally 
extern FTYPE reverse(int, int);
extern FTYPE ran3(int *);
extern intAVec intset(int);
extern void charges(FArrayND &, particle_dist *, fluid *, int);
/* extern void quasineutral_den_flux(FArrayND_ranged &,FArrayND_ranged &, */
/* 				  particle_dist *, fluid *, int); */
extern void quasineutral_den_flux(FArrayND_ranged &,
				  INDICIES(FArrayND_ranged &,
					   FArrayND_ranged &,
					   FArrayND_ranged &),
				  particle_dist *, fluid *, int);
extern void quasineutral_moments(FArrayND_ranged &,
				 INDICIES(FArrayND_ranged &,FArrayND_ranged &,FArrayND_ranged &),
				 particle_dist *, fluid *, int);
FTYPE wteval_df3(FTYPE, FTYPE, int);
void terminate(int n, const char *message);



extern FTYPE *workarray; // A workspace array will be declared once the array sizes are known.
extern int nworkarray;

#include "eppic-times.h"

#endif
