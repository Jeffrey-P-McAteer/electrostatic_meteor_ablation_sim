/* 
   Print a string corresponding to the PETSc convergence/divergence
   reason returned by KSPGetConvergedReason 
   
   PETSc may have a routine for this. I just copied the info from
   the manual entry for KSPConvergedReason
*/

#include <iostream>
#include <string.h>

using namespace std;

#if HAVE_PETSC
#include "petsc.h"

PetscErrorCode interpret_reason(int reason, string& strreason)
{

  //  string strreason;
  PetscErrorCode perr;

  if (reason==0) strreason="KSP_CONVERGED_ITERATING";

  /* Convergence */
  else if (reason==1) strreason="KSP_CONVERGED_RTOL_NORMAL";
  else if (reason==9) strreason="KSP_CONVERGED_ATOL_NORMAL";
  else if (reason==2) strreason="KSP_CONVERGED_RTOL";
  else if (reason==3) strreason="KSP_CONVERGED_ATOL";
  else if (reason==4) strreason="KSP_CONVERGED_ITS";
  else if (reason==5) strreason="KSP_CONVERGED_CG_NEG_CURVE";
  else if (reason==6) strreason="KSP_CONVERGED_CG_CONSTRAINED";
  else if (reason==7) strreason="KSP_CONVERGED_STEP_LENGTH";
  else if (reason==8) strreason="KSP_CONVERGED_HAPPY_BREAKDOWN";

  /* Divergence */
  else if (reason==-2) strreason="KSP_DIVERGED_NULL";
  else if (reason==-3) strreason="KSP_DIVERGED_ITS";
  else if (reason==-4) strreason="KSP_DIVERGED_DTOL";
  else if (reason==-5) strreason="KSP_DIVERGED_BREAKDOWN";
  else if (reason==-6) strreason="KSP_DIVERGED_BREAKDOWN_BICG";
  else if (reason==-7) strreason="KSP_DIVERGED_NONSYMMETRIC";
  else if (reason==-8) strreason="KSP_DIVERGED_INDEFINITE_PC";
  else if (reason==-9) strreason="KSP_DIVERGED_NANORINF";
  else if (reason==-10) strreason="KSP_DIVERGED_INDEFINITE_MAT";

  else strreason="No reason found";

  /* Print Reason */
  //  perr=PetscPrintf(PETSC_COMM_WORLD,"KSPGetConvergedReason: %d (%s)\n",reason,strreason.c_str());CHKERRQ(perr);
  /*
  perr=PetscPrintf(PETSC_COMM_SELF,
		   "proc[%d] KSPGetConvergedReason: %d (%s)\n",
		   mpi_rank,reason,strreason.c_str());CHKERRQ(perr);
  */
  return(0);
}

#else
/* If not using PETSc, here is a dummy routine. */
typedef int PetscErrorCode;

PetscErrorCode interpret_reason(int reason, string& strreason)
{ return 0; }

#endif
