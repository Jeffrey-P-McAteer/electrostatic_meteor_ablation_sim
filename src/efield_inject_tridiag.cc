

//----------------------------------------------------------------------------
//   Solves for the electric potential (phi) on a 2D (or 3D) grid.
//   It uses a spectral method in y (and z) and a tridiagonal solver in x for
//   each ky (and kz). This allows for non-periodic solutions to phi in the x 
//   direction. The possible conditions in x are open (natural), Dirichlet, 
//   Neumann and periodic. 
//
//   INPUTS:
//   Efield -- field    -- See eppic-efield.h for definition.
//   rho    -- FArrayND -- passes in the charge density; this is a working 
//                         array and will be over written. 
//
//   Notes:
//   Efield.boundary[0] and Efield.boundary[1] define the left and right x
//   boundary conditions, respectively. Efield.phi_rho holds the solution. 
//
//   Since this algorithm is long, each step is broken into several helper 
//   routines. For more information on the algorithm, please see Birdsall and
//   Langdon , Plasma Physics via Computer Simulation, page 318 (First edition,
//   published 1991).
//----------------------------------------------------------------------------

#include "eppic.h"
#include "eppic-times.h" // to allow clocking 
#include "eppic-mpi.h" // to know about domain decomposition 
#include "eppic-io.h" // for terminate
#include "efield_inject_tridiag.h" // for helper functions 
#include "timerFile.h"

void efield_inject_tridiag(field &Efield , FArrayND &rho)
{
  // This is stuff for another timer, one specific to efield_inject_tridiag
  static timerFile timer;
  clock_t before_time;
  static bool firstEntry=true;
  if (firstEntry) {
    firstEntry = false;
    init_timerFile(timer,"efield_times.dat",nout,0);
  }

  
  index_type array_size[]={Efield.ny,Efield.nz}; 
  //  cout <<mpi_rank << ": rfftw_many("  <<array_size[0] <<","<<array_size[1]<<","<<Efield.nx<<")"<<endl; 

  MPI_Barrier(MPI_COMM_WORLD);
  static rfftw_many working_array = 
    rfftw_many(subdomain.internal_comm,array_size,Efield.nx);
  static rfftw_many boundary_working_array=
    rfftw_many(subdomain.internal_comm,array_size,3);
  MPI_Barrier(MPI_COMM_WORLD);

  before_time=times(&timer.refTimer);
  efield_inject_tridiag_init(rho,Efield,working_array);
  timer.times["init"] += times(&timer.refTimer)-before_time;

  before_time=times(&timer.refTimer);
  working_array.transform();
  timer.times["1st_transform"] += times(&timer.refTimer)-before_time;
  
  before_time=times(&timer.refTimer);
  if ((Efield.boundary[0]==open)&&(Efield.boundary[1]==open))
    efield_inject_tridiag_solve_open(working_array,
				     boundary_working_array,
                                     Efield);
  else if ((Efield.boundary[0]==dirichlet)&&(Efield.boundary[1]==dirichlet))
    efield_inject_tridiag_solve_dirichlet(working_array,
					  boundary_working_array,
					  Efield);
  else if ((Efield.boundary[0]==neumann)&&(Efield.boundary[1]==dirichlet))
    efield_inject_tridiag_solve_mixed(working_array,
					  boundary_working_array,
					  Efield);
  else if ((Efield.boundary[0]==dirichlet)&&(Efield.boundary[1]==neumann))
    terminate(-1,"Dirichlet,Neumann not implemented, try Neumann,Dirichlet.");
  else if ((Efield.boundary[0]==neumann)&&(Efield.boundary[1]==neumann))
    terminate(-1,"Neumann,Neumann not implemented, try Neumann,Dirichlet.");
  else 
    efield_inject_tridiag_solve_periodic(working_array,
					 boundary_working_array,
					 Efield);
  timer.times["solve"] += times(&timer.refTimer)-before_time;
#if NDIM==2
  before_time=times(&timer.refTimer);
  efield_inject_tridiag_sync_subdomain_batch(working_array,
				       boundary_working_array,
				       Efield);
  timer.times["2dsync"] += times(&timer.refTimer)-before_time;
#endif

  before_time=times(&timer.refTimer);
  working_array.invtransform();

  boundary_working_array.invtransform();
  timer.times["2nd_transform"] += times(&timer.refTimer)-before_time;

  before_time=times(&timer.refTimer);
  efield_inject_tridiag_tophi(working_array,boundary_working_array,Efield);
  timer.times["tophi"] += times(&timer.refTimer)-before_time;

  flush_timerFile(timer,-1);

}




