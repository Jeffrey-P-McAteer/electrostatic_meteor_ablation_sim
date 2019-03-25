/*
  Simple routine to echo the PETSc options file so it will show up in the output log.
  This is designed to make it easier to keep track of which PETSc options were
  used to create a particular EPPIC run. Initially, this routine will assume that
  the options file was written correctly (e.g. -petsc_opt_name value). We may want 
  to add more rigorous error-checking later.

  Started 31Jul14 (may)
*/

#include <cstdio>
#include <iostream>
#include <string>
#include "eppic.h"

#if HAVE_PETSC
#include "petsc.h"

PetscErrorCode echo_petsc(char *petscOpts) {
  FILE *fp;
  int nLines,strLen=130;
  char opt[strLen],val[strLen],line[strLen],*fstatus,*space;
  char *fname=petscOpts;


  /* Check for file and open */
  if ((fp=fopen(petscOpts,"r"))==NULL) {
    if (mpi_rank==0 && subdomain.rank==0) printf("Cannot open %s for reading",fname);
    std::terminate();
  }
    
  /* Print leading message */
  if (mpi_rank==0 && subdomain.rank==0) {
    printf(";======================= PETSc =======================\n");
    printf(";Echoing options database input file...\n");
    printf(";Filename: %s\n",fname);
  }
  /* Read line */
  fstatus=fgets(line,strLen,fp);
  while (fstatus != (char *)EOF && fstatus != (char *)NULL) {
    if (mpi_rank==0 && subdomain.rank==0) printf(";%s",line);
    fstatus=fgets(line,strLen,fp);
  }
  
  //  if (mpi_rank==0 && subdomain.rank==0) printf(";\nEND PETSc\n\n");
  if (mpi_rank==0 && subdomain.rank==0)
    printf("\n;=====================================================\n\n");

  /* Close the file */ 
  fclose(fp);
  fflush(NULL);
  return 0;
}
#else
typedef int PetscErrorCode;
PetscErrorCode echo_petsc() { return 0; }

#endif
