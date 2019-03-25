/* Terminate the program elegantly */

#include <csignal>
#include <stdio.h>
#include "eppic.h"
#include "eppic-mpi.h"
#include "eppic-types.h"

void terminate(int n, const char *message)
{
  if (n != 0 || mpi_rank == 0) {
    printf("\n%s\n\nEPPIC process %d terminating all %d processes\n",
	   message, mpi_rank, mpi_np);
    fprintf(stderr,"\n%s\n\nEPPIC process %d terminating all %d processes\n",
	    message, mpi_rank, mpi_np);
  }
  
#ifdef CHECK
  // Raise a trap signal to allow a debugger to catch this problem.
  if (n!=0) raise(SIGTRAP);
#endif

#ifdef USE_MPI
  int mpi_err;
  if (n != 0) MPI_Abort(MPI_COMM_WORLD,n);
  mpi_err=MPI_Finalize();
  if (mpi_rank == 0) printf("terminate.cc\tMPI finalized with mpi_err=%d\n",
			    mpi_err);
#endif

  exit(n);
}
