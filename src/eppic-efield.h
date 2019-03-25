
/* This header holds declarations for all routines involved in solving 
   for the electric potential, and thus the electric field. There are 
   several different algorithms, and some may be solved for non-periodic 
   conditions. 
*/

#ifndef EPPIC_EFIELD_H
#define EPPIC_EFIELD_H

#include "eppic-types.h"
#include "eppic-pic.h"
#include "eppic-fluid.h"


extern int no_parallel_efields; // If no_parallel_efields==TRUE, damp out 
                                // modes where ikx=0
extern int efield_algorithm; /* which method to use for calc. of phi */
extern int efield_zero_dc; // force total neutrality when calculating phi
extern int electron_dist; /* distribution number from which parameters are obtained */
extern FTYPE Ex0_external; 
extern FTYPE Ex0_rate; 
extern FTYPE Ey0_external; 
extern FTYPE Ey0_rate; 
extern FTYPE Ez0_external;
extern FTYPE Ez0_rate;

// for inject_tridiag routines, to cap number of tridiags solved at a time
extern int MAX_TRIDIAG_SOLVE;


/** bc_type defines the boundary conditions in non-periodic simulations.
 */
enum bc_type{
  // Changed by Glenn 3/1/2018 to be consistent with boundary_type for efield
  open = -1, /**< Laplace's eqn outside simulation region */
  periodic = 0, /**< phi(-1,k) or phi(Nx,k) = 0 */
  dirichlet = 1, /**< Grad_x phi(0,k) or Grad_x(Nx,k) = 0 */
  neumann = 2   /**< phi(0,k) = phi(Nx,k) */
  //open = -1, /**< Laplace's eqn outside simulation region */
  //dirichlet = 0, /**< phi(-1,k) or phi(Nx,k) = 0 */
  //neumann = 1,   /**< Grad_x phi(0,k) or Grad_x(Nx,k) = 0 */
  //periodic = 2  /**< phi(0,k) = phi(Nx,k) */
};

/** This structure holds the definition of the electric potential and
    electric field. */
typedef struct { 
  FArrayND_ranged phi_rho; /**< Electric potential with ghost cells.
			      See phix_guard_size for more info about ghost 
			      cells. */
  FArrayND Ex; ///< Electric field in x 
  FArrayND Ey; ///< Electric field in y 
  FArrayND Ez; ///< Electric field in z
  int nx;  ///< number of grid cells in x 
  int ny;  ///< number of grid cells in y 
  int nz;  ///< number of grid cells in z 
  FTYPE dx;///< grid step size in meters in x
  FTYPE dy;///< grid step size in meters in y 
  FTYPE dz;///< grid step size in meters in z
  FTYPE eps;///< dielectric constant (default=8.8542e-12) 
  bc_type boundary[2]; ///< boundary conditions on the left side (0) and right side (1) 
  FTYPE  boundary_values[4]; ///< values of the boundary conditions, left(0,1),
  //                                                         right(2,3)
  FTYPE Ex0_amplitude;
  FTYPE Ey0_amplitude;
  FTYPE Ez0_amplitude;
  FTYPE Ex0_rate;
  FTYPE Ey0_rate;
  FTYPE Ez0_rate;
  int algorithm;
  bool zero_dc;
} field;

/** Defines the number of guard cells for the potential. 
    This is only applicable when Domain decomposition is used. */
#ifdef USE_DOMAINS
static const int phix_guard_size[]={1,2};
#else
static const int phix_guard_size[]={0,0};
#endif
/** Defines the number of guard cells for quasineutral density
    and flux divergence. This is applicable regardless of whether
    domain decomposition is used. */
/* #if USE_QN // <-- Necessary? 25Jan2017 (may). */
static const int qnx_guard_size[]={1,1};
/* #endif */

/// Initializes the efield field structure and the corresponding rho array
void init_field(field &Efield, /**< [in,out] electric field */ 
		FArrayND &rho, /**< [in,out] charge density array */
		int nx, /**< [in] number of grid cells in x */
		int ny, /**< [in] number of grid cells in y */
		int nz,/**< [in] number of grid cells in z */
		FTYPE dx, /**< [in] grid spacing in meters, x direction */
		FTYPE dy, /**< [in] grid spacing in meters, y direction */
		FTYPE dz, /**< [in] grid spacing in meters, z direction */
		FTYPE eps, ///< [in] dielectric constant (default=8.8542e-12) 
		bc_type boundary[2], /**< [in] boundary conditions for 
				       left (0) and right (1) sides. 
				       See bc_type.
				    */
		FTYPE boundary_values[4], /**< [in] see field def. */
		FTYPE Ex0_amplitude,
		FTYPE Ey0_amplitude,
		FTYPE Ez0_amplitude,
		FTYPE Ex0_rate,
		FTYPE Ey0_rate,
		FTYPE Ez0_rate,
		int efield_algorithm,
		int efield_zero_dc
		);

/// Initializes the efield field structure, 
/// and the corresponding quasineutral density and flux-divergence arrays
void init_field(field &Efield, /**< [in,out] electric field */ 
		FArrayND_ranged &qden, /**< [in,out] quasineutral density */
		FArrayND_ranged &divG, /**< [in,out] flux divergence */
		int nx, /**< [in] number of grid cells in x */
		int ny, /**< [in] number of grid cells in y */
		int nz,/**< [in] number of grid cells in z */
		FTYPE dx, /**< [in] grid spacing in meters, x direction */
		FTYPE dy, /**< [in] grid spacing in meters, y direction */
		FTYPE dz, /**< [in] grid spacing in meters, z direction */
		FTYPE eps, ///< [in] dielectric constant (default=8.8542e-12) 
		bc_type boundary[2], /**< [in] boundary conditions for 
				       left (0) and right (1) sides. 
				       See bc_type.
				    */
		FTYPE boundary_values[4], /**< [in] see field def. */
		FTYPE Ex0_amplitude,
		FTYPE Ey0_amplitude,
		FTYPE Ez0_amplitude,
		FTYPE Ex0_rate,
		FTYPE Ey0_rate,
		FTYPE Ez0_rate,
		int efield_algorithm,
		int efield_zero_dc
		);

/// Initializes the efield field structure, 
/// and the corresponding quasineutral density and flux arrays
void init_field(field &Efield, /**< [in,out] electric field */ 
		FArrayND_ranged &qden, /**< [in,out] quasineutral density */
		INDICIES(FArrayND_ranged &Gx,  /**< [in,out] x flux */
			 FArrayND_ranged &Gy,  /**< [in,out] y flux */
			 FArrayND_ranged &Gz), /**< [in,out] z flux */
		int nx, /**< [in] number of grid cells in x */
		int ny, /**< [in] number of grid cells in y */
		int nz,/**< [in] number of grid cells in z */
		FTYPE dx, /**< [in] grid spacing in meters, x direction */
		FTYPE dy, /**< [in] grid spacing in meters, y direction */
		FTYPE dz, /**< [in] grid spacing in meters, z direction */
		FTYPE eps, ///< [in] dielectric constant (default=8.8542e-12) 
		bc_type boundary[2], /**< [in] boundary conditions for 
				       left (0) and right (1) sides. 
				       See bc_type.
				    */
		FTYPE boundary_values[4], /**< [in] see field def. */
		FTYPE Ex0_amplitude,
		FTYPE Ey0_amplitude,
		FTYPE Ez0_amplitude,
		FTYPE Ex0_rate,
		FTYPE Ey0_rate,
		FTYPE Ez0_rate,
		int efield_algorithm,
		int efield_zero_dc
		);


/* void efield_wrapper(field &Efield, */
/* 		    FArrayND &rho,FArrayND_ranged &qden, */
/* 		    FArrayND_ranged &divG, */
/* 		    FTYPE itime); */
void efield_wrapper(field &Efield,
		    FArrayND &rho,FArrayND_ranged &qden,
		    INDICIES(FArrayND_ranged &Gx,
			     FArrayND_ranged &Gy,
			     FArrayND_ranged &Gz),
		    FTYPE itime);

/** A modified version of Marcos's efield_inject algorithm. 
    The original is probably very similar, and can be found in 
    efield_inject_orig 
*/
void efield_inject(
		   field &Efield, /**< [in,out] electric field */ 
		   FArrayND &rho /**< [in,out] charge density */
		   );


/** The original solver written by Marcos. Modifications were made to allow
    the routine to be called independently of efield. This involved 
    changing the input parameter from phi to Efield, to mimic the other 
    Poisson equation solvers. Also, calls to the timing routines were added */
void efield_inject_orig(
		   field &Efield, /**< [in,out] electric field */ 
		   FArrayND &rho, /**< [in,out] charge density */
		   int bc);



/** This is a first attempt at parallelizing efield_inject. */ 
void efield_inject_parallel(
		   field &Efield, /**< [in,out] electric field */ 
		   FArrayND &rho /**< [in,out] charge density */
		   );


/** This is a better parallelization and improvement on the efield_inject
    routine. It adds more possibilities for boundary conditions, see bc_type. 
*/
void efield_inject_tridiag(
		   field &Efield, /**< [in,out] electric field */ 
		   FArrayND &rho /**< [in,out] charge density */
		   );

void efield_inject_tridiag_new(
		   field &Efield, /**< [in,out] electric field */ 
		   FArrayND &rho /**< [in,out] charge density */
		   );


/** This is the original spectral, periodic efield solver. */
void efield(
		   field &Efield, /**< [in,out] electric field */ 
		   FArrayND &rho /**< [in,out] charge density */
		   );


#include "eppic-pic.h"
#include "eppic-fluid.h"

/** This routine solves for the potential using the electron momentum
    equation with the quasi-neutral approximation */
// -->Consider deprecating. 25Jan2017 (may)
void efield_quasineut(// [in,out] electric field
		      field &Efield,
		      // [in] vector of particle distributions
		      particle_dist *pic,
		      // [in] vector of fluid distributions
		      fluid *fspecie);

/** This is a helper routine for the efield_quasineut algorithm. It does
    the actual solving once all of the temporary arrays are defined. */
// -->Consider deprecating. 25Jan2017 (may)
int efield_quasineut_solver(// [in,out] electric field
			    field &Efield,
			    // [in] total charge density
			    FArrayND_ranged den,
			    // [in] total flux divergence
			    //      (excl. quasineutral, inertialess fluids)
			    FArrayND_ranged divG);

/** This routine solves the quasineutral potential equation with static,
    inertialess, fluid electrons and PIC ions, using PETSc.

    It should be able to accommodate other (inertial) fluid and 
    PIC distributions but has not been tested for such yet. 25Jan2017 (may) */
/* int efield_quasineut_petsc(// [in,out] electric field */
/* 			   field &Efield, */
/* 			   // [in] total charge density */
/* 			   FArrayND_ranged den, */
/* 			   // [in] total flux divergence */
/* 			   //      (excl. quasineutral, inertialess fluids) */
/* 			   FArrayND_ranged divG); */
int efield_quasineut_petsc(// [in,out] electric field
			   field &Efield,
			   // [in] total charge density
			   FArrayND_ranged den,
			   // [in] x, y, & z fluxes
			   INDICIES(FArrayND_ranged Gx,
				    FArrayND_ranged Gy,
				    FArrayND_ranged Gz));

void efield_test_phi_rho(field &Efield, FArrayND rho, eppic_system sys, int it);

void efield_p3dfft(field &Efield, FArrayND &rho);

void efield_multigrid(field &Efield, FArrayND &rho);

void efield_glenntest( field &Efield, FArrayND &rho);

#endif
