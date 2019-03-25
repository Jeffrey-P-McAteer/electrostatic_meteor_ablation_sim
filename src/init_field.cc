/* Initialize E and phi, and either rho or qden & divG: */

#include "eppic.h"
#include "eppic-mpi.h"

// For quasineutral runs (v1), which need qden & divG
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
		)
{
  
#ifdef USE_DOMAINS
  // Phi needs 2 guard cells on the right side but only one on the left side:
  const int nsize_phi[]={INDICIES(nx+phix_guard_size[0]+phix_guard_size[1],ny,nz)};
  int nstart_phi[]={INDICIES(0-phix_guard_size[0],0,0)};
  Efield.phi_rho =  FArrayND_ranged(nstart_phi, nsize_phi);
#else
  // It is more efficient to use electric fields only in 1 and 2-D
  const int nsize_phi[]={INDICIES(nx*nsubdomains,ny,nz)};
  Efield.phi_rho = FArrayND(nsize_phi);
#endif
  const int nsize_qn[]={INDICIES(nx+qnx_guard_size[0]+qnx_guard_size[1],ny,nz)};
  int nstart_qn[]={INDICIES(0-qnx_guard_size[0],0,0)};
  qden = FArrayND_ranged(nstart_qn,nsize_qn);
  divG = FArrayND_ranged(nstart_qn,nsize_qn);

  // Allocate workarray space:
  //  nworkarray=rho.size();
  //  workarray = new FTYPE[nworkarray];

  if (ndim == 1) Efield.Ex =  FArrayND(INDICIES(nx, ny, nz));
  Efield.nx = nx;
  Efield.ny = ny;
  Efield.nz = nz;
  Efield.dx = dx;
  Efield.dy = dy;
  Efield.dz = dz;
  Efield.eps = eps;
  Efield.boundary[0] = boundary[0];
  Efield.boundary[1] = boundary[1];
  /* if either boundary is periodic, make both periodic */
  if ((Efield.boundary[0]==periodic) ||
      (Efield.boundary[1]==periodic)) {
    Efield.boundary[0] = periodic;
    Efield.boundary[1] = periodic;
  }

  /* deal with Ex0_external */
  switch(Efield.boundary[0]) {
  case(open):
    Efield.boundary_values[0] = 0.0;
    Efield.boundary_values[0] = 0.0;
    Efield.boundary_values[0] = 0.0;
    Efield.boundary_values[0] = 0.0;
    break;
  case(dirichlet):

    /* the ny factor comes from the fact that the fftw routines do not 
       normalize the forward and reverse ffts, so we do the normalization
       in the y direction;
    */
    Efield.boundary_values[0] = Ex0_external*(nx*nsubdomains)*dx*ny;
    Ex0_external=0;
    break;
  case(neumann):
    break;
  }

  Efield.Ex0_amplitude = Ex0_external;
  Efield.Ey0_amplitude = Ey0_external;
  Efield.Ez0_amplitude = Ez0_external;
  Efield.Ex0_rate = Ex0_rate;
  Efield.Ey0_rate = Ey0_rate;
  Efield.Ez0_rate = Ez0_rate;
  Efield.algorithm = efield_algorithm;
  Efield.zero_dc = static_cast<bool>(efield_zero_dc);
  if (Efield.algorithm == 5) 
    if (Efield.boundary[0] == periodic)
      Efield.zero_dc = true;
}

// For quasineutral runs (v2), which need qden & fluxes (Gx, Gy, & Gz)
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
		)
{
  
#ifdef USE_DOMAINS
  // Phi needs 2 guard cells on the right side but only one on the left side:
  const int nsize_phi[]={INDICIES(nx+phix_guard_size[0]+phix_guard_size[1],ny,nz)};
  int nstart_phi[]={INDICIES(0-phix_guard_size[0],0,0)};
  Efield.phi_rho =  FArrayND_ranged(nstart_phi, nsize_phi);
#else
  // It is more efficient to use electric fields only in 1 and 2-D
  const int nsize_phi[]={INDICIES(nx*nsubdomains,ny,nz)};
  Efield.phi_rho = FArrayND(nsize_phi);
#endif
  const int nsize_qn[]={INDICIES(nx+qnx_guard_size[0]+qnx_guard_size[1],ny,nz)};
  int nstart_qn[]={INDICIES(0-qnx_guard_size[0],0,0)};
  qden = FArrayND_ranged(nstart_qn,nsize_qn);
  Gx = FArrayND_ranged(nstart_qn,nsize_qn);
#if NDIM > 1
  Gy = FArrayND_ranged(nstart_qn,nsize_qn);
#if NDIM > 2
  Gz = FArrayND_ranged(nstart_qn,nsize_qn);
#endif
#endif

  // Allocate workarray space:
  //  nworkarray=rho.size();
  //  workarray = new FTYPE[nworkarray];

  if (ndim == 1) Efield.Ex =  FArrayND(INDICIES(nx, ny, nz));
  Efield.nx = nx;
  Efield.ny = ny;
  Efield.nz = nz;
  Efield.dx = dx;
  Efield.dy = dy;
  Efield.dz = dz;
  Efield.eps = eps;
  Efield.boundary[0] = boundary[0];
  Efield.boundary[1] = boundary[1];
  /* if either boundary is periodic, make both periodic */
  if ((Efield.boundary[0]==periodic) ||
      (Efield.boundary[1]==periodic)) {
    Efield.boundary[0] = periodic;
    Efield.boundary[1] = periodic;
  }

  /* deal with Ex0_external */
  switch(Efield.boundary[0]) {
  case(open):
    Efield.boundary_values[0] = 0.0;
    Efield.boundary_values[0] = 0.0;
    Efield.boundary_values[0] = 0.0;
    Efield.boundary_values[0] = 0.0;
    break;
  case(dirichlet):

    /* the ny factor comes from the fact that the fftw routines do not 
       normalize the forward and reverse ffts, so we do the normalization
       in the y direction;
    */
    Efield.boundary_values[0] = Ex0_external*(nx*nsubdomains)*dx*ny;
    Ex0_external=0;
    break;
  case(neumann):
    break;
  }

  Efield.Ex0_amplitude = Ex0_external;
  Efield.Ey0_amplitude = Ey0_external;
  Efield.Ez0_amplitude = Ez0_external;
  Efield.Ex0_rate = Ex0_rate;
  Efield.Ey0_rate = Ey0_rate;
  Efield.Ez0_rate = Ez0_rate;
  Efield.algorithm = efield_algorithm;
  Efield.zero_dc = static_cast<bool>(efield_zero_dc);
  if (Efield.algorithm == 5) 
    if (Efield.boundary[0] == periodic)
      Efield.zero_dc = true;
}

// For non-quasineutral runs, which need rho
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
		)
{
#ifdef USE_DOMAINS
  // Phi needs 2 guard cells on the right side but only one on the left side:
  const int nsize_phi[]={INDICIES(nx+phix_guard_size[0]+phix_guard_size[1],
                                  ny+yguard_size,nz+zguard_size)};
  int nstart[]={INDICIES(-phix_guard_size[0],0,0)};
  Efield.phi_rho =  FArrayND_ranged(nstart, nsize_phi);
  const int nsize_rho[]={INDICIES(nx+xguard_size,ny+yguard_size,nz+zguard_size)};
  rho =  FArrayND(nsize_rho);
#else
  // It is more efficient to use electric fields only in 1 and 2-D
  const int nsize[]={INDICIES(nx*nsubdomains,ny,nz)};
  Efield.phi_rho = FArrayND(nsize);
  rho =  FArrayND(nsize);
#endif

  // Allocate workarray space:
  //  nworkarray=rho.size();
  //  workarray = new FTYPE[nworkarray];

  if (ndim == 1) Efield.Ex =  FArrayND(INDICIES(nx, ny, nz));
  Efield.nx = nx;
  Efield.ny = ny;
  Efield.nz = nz;
  Efield.dx = dx;
  Efield.dy = dy;
  Efield.dz = dz;
  Efield.eps = eps;
  Efield.boundary[0] = boundary[0];
  Efield.boundary[1] = boundary[1];
  /* if either boundary is periodic, make both periodic */
  if ((Efield.boundary[0]==periodic) ||
      (Efield.boundary[1]==periodic)) {
    Efield.boundary[0] = periodic;
    Efield.boundary[1] = periodic;
  }

  /* deal with Ex0_external */
  switch(Efield.boundary[0]) {
  case(open):
    Efield.boundary_values[0] = 0.0;
    Efield.boundary_values[0] = 0.0;
    Efield.boundary_values[0] = 0.0;
    Efield.boundary_values[0] = 0.0;
    break;
  case(dirichlet):

    /* the ny factor comes from the fact that the fftw routines do not 
       normalize the forward and reverse ffts, so we do the normalization
       in the y direction;
    */
    Efield.boundary_values[0] = Ex0_external*(nx*nsubdomains)*dx*ny;
    Ex0_external=0;
    break;
  case(neumann):
    break;
  }

  Efield.Ex0_amplitude = Ex0_external;
  Efield.Ey0_amplitude = Ey0_external;
  Efield.Ez0_amplitude = Ez0_external;
  Efield.Ex0_rate = Ex0_rate;
  Efield.Ey0_rate = Ey0_rate;
  Efield.Ez0_rate = Ez0_rate;
  Efield.algorithm = efield_algorithm;
  Efield.zero_dc = static_cast<bool>(efield_zero_dc);
  if (Efield.algorithm == 5) 
    if (Efield.boundary[0] == periodic && !Efield.zero_dc){
      Efield.zero_dc = true;
      if (mpi_rank==0) printf("Warning! Periodic boundary in x but the dc field was not set to 0.  We are fixing that...\n");
    }
}
