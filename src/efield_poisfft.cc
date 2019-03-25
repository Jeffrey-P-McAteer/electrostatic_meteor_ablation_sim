 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Solves for phi, the electric potential, using Poisson's equation.
    
    We use the PETSc library, which allows for distributed solutions to
    linear (used in this case, ie Ax=b), non-linear, and matrix free 
    problems. Currently this will use the multigrid petsc method.
    
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#include "eppic.h"

#ifdef USE_MPI
#include "eppic-mpi.h"
#endif


#if HAVE_POISFFT
#include "poisfft.h"

void efield_poisfft(field &Efield, FArrayND &rho)
{

  int *BCs;
  PoisFFT::Solver<3,FTYPE> pois_solver(ns, Ls, BCs)
                                       {PoisFFT::DIRICHLET, 
                                           PoisFFT::PERIODIC, 
                                           PoisFFT::PERIODIC});
                                       
}
#else
void efield_poisfft(field &Efield, FArrayND &rho){
  terminate(-1, "Error: efield_multigrid can only handle NDIM=3 and HAVE_PETSC=1");
}
#endif
