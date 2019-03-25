
#include "eppic-efield.h"
#include "eppic-math.h"
#include "eppic-io.h"
#include "eppic.h"

// void efield_wrapper(field &Efield,
// 		    FArrayND &rho,FArrayND_ranged &qden,FArrayND_ranged &divG,
// 		    FTYPE itime) 
void efield_wrapper(field &Efield,
		    FArrayND &rho,FArrayND_ranged &qden,
		    INDICIES(FArrayND_ranged &Gx, \
			     FArrayND_ranged &Gy, \
			     FArrayND_ranged &Gz),
		    FTYPE itime) 
{

  /* Log the amount of time spent on this routine */
  clock_t start_time=times(&times_buf);

  // update External fields
  if (Efield.Ex0_rate > 0)
    Ex0_external = Efield.Ex0_amplitude*cos(2*PI*Efield.Ex0_rate*itime);
  if (Efield.Ey0_rate > 0)
    Ey0_external = Efield.Ey0_amplitude*cos(2*PI*Efield.Ey0_rate*itime);
  if (Efield.Ez0_rate > 0)
    Ez0_external = Efield.Ez0_amplitude*cos(2*PI*Efield.Ez0_rate*itime);

  // calculated Field on grid
  switch(Efield.algorithm)
    {
    case 0:
      efield(Efield,rho);
      break;
    case 1:
      if (mpi_rank==0) {
	cout << "\tefield_algorithm 1 is currently unused.\n"
	     << "\tPlease choose a different value for efield_algorithm.\n"
	     << "\tSee efield_wrapper for available routines.\n";
	terminate(TRUE,"\n");
      }
      break;
    case 2:

      // Make sure PETSc is compiled in
      if (HAVE_PETSC != 1) {
	if (mpi_rank==0) {
	  cout << "\tefield_quasineut_petsc requires PETSc. \n"
	       << "\tPlease reconfigure EPPIC with PETSc and rerun\n";
	  terminate(TRUE,"\n");
	}
      }
#if HAVE_PETSC
      // If this rank is on the PETSc subcomm, proceed to field solver
      if (mpi_rank < petsc_np) {
	int status=0;
	// status = efield_quasineut_petsc(Efield,qden,divG);
	status = efield_quasineut_petsc(Efield,qden,
					INDICIES(Gx,Gy,Gz));
	if (status != 0) {
	  if (mpi_rank == 0) 
	    printf("%s:%d efield_quasineut_petsc exited with status %d\n",
				    __func__,__LINE__,status);
	}
      }
      // Let all ranks catch up to each other
      MPI_Barrier(MPI_COMM_WORLD);

      // Broadcast phi to all ranks
      if (petsc_np < mpi_np) {
	MPI_Bcast((void*)&(Efield.phi_rho(INDICIES(0-phix_guard_size[0],0,0))),
		  Efield.phi_rho.size(),
		  MPI_FTYPE,
		  0,
		  MPI_COMM_WORLD);
      }
#endif

      break;
    case 3:
      efield_inject(Efield,rho);
      break;
    case 4:
      efield_inject_parallel(Efield,rho);
      break;
    case 5:
      efield_inject_tridiag(Efield,rho);
      break;
    case 6:
      efield_inject_orig(Efield,rho,0);
      break;
    case 7:
      efield_p3dfft(Efield, rho);
      break;
    case 8:
      // No Efield
      break;
    case 9:
      efield_multigrid(Efield, rho);
      break;
    case 10:
      efield_glenntest(Efield, rho);
      break;
    }

#ifdef DEBUG
  eppic_system sys;
  sys.dx=dx;
  sys.dy=dy;
  sys.dz=dz;
  sys.nx=nx;
  sys.ny=ny;
  sys.nz=nz;
  sys.eps=eps;

  /*  if ((fsteep == 0) and (fwidth ==0)) {
    efield_test_phi_rho(Efield,rho,sys,static_cast<int>(itime/dt));
  } 
  */
#endif

  /* Timer */
  efield_time += times(&times_buf)-start_time;
}
