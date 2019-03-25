#include <stdio.h>
// #include <mpi.h>
#include "eppic-mpi.h"
#include "eppic.h"

void mpi_check(int error_code)
{

  if (error_code != MPI_SUCCESS) {
    char code_string[BUFSIZ],class_string[BUFSIZ];
    int code_strlen,class_strlen,error_class;

    MPI_Error_class(error_code,&error_class);
    MPI_Error_string(error_class,class_string,&class_strlen);
    MPI_Error_string(error_code,code_string,&code_strlen);
    printf("[%d] MPI_ERROR %s: %s\n",mpi_rank,class_string,code_string);
    fflush(stdout);

    MPI_Abort(MPI_COMM_WORLD,error_code);
  }

}
