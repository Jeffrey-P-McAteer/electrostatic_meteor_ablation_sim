#include "eppic.h"
#ifdef USE_MPI
#include "eppic-mpi.h"
#include <stdio.h> 

void mpi_error(int mpi_rank, int mpi_err, char *mpi_msg) {
  
  int stringlen=MPI_MAX_ERROR_STRING;
  char string[MPI_MAX_ERROR_STRING];
  printf("ERROR in processor %d: MPI error %d from %s\n",
	 mpi_rank,mpi_err,mpi_msg);
  fprintf(stderr,"ERROR in processor %d: MPI error %d from %s\n",
	  mpi_rank,mpi_err,mpi_msg);
  MPI_Error_string(mpi_err, string, &stringlen);
  printf("                             : %s\n", string);
  terminate(-1,"output termination");
  
}
#endif
