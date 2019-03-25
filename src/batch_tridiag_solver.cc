

#include "eppic-types.h"
#include "tridiag.h"

/* 
   Algorithm taken from Press et al., Numerical Recipies.
   Solves for x in Ax=b, where A is a tridiagonal matrix.
   INPUTS:
   a -- ArrayNd<FTYPE,1> -- lower stripe of tridiag matirx
   b -- ArrayNd<FTYPE,1> -- diagonal, middle stripe of tridiag matrix
   c -- ArrayNd<FTYPE,1> -- upper stripe of tridiag matrix
   rhs -- ArrayNd<FTYPE,1> -- corresponds to 'b' in Ax=b 
   lhs -- ArrayNd<FTYPE,1> -- solution and corresponds to 'x' in Ax=b

   Note: all inputs must be predefined and have the same dimensions. There 
   is no test done to verify this.
*/

void batch_tridiag_solver(ArrayNd<FTYPE,2> &a,ArrayNd<FTYPE,2> &b, 
			  ArrayNd<FTYPE,2> &c,
			  ArrayNd<FTYPE,2>&rhs, ArrayNd<FTYPE,2> &lhs) 
{
  int n_solves = rhs.size(0);  
  int n_unknown = rhs.size(1);  
  FTYPEAVec gam=FTYPEAVec(n_unknown)=0;
  int j;
  FTYPE bet;

  for (int isolve=0;isolve<n_solves; isolve++) {
    // We need to do a forward sweep
    bet=b(isolve,0); 
    gam=0;
    if (b[0] == 0.0) 
      throw("Error 1 in tridiag_solver");
    
    lhs(isolve,0) = rhs(isolve,0)/bet;
    
    /* decomposition and forward substitution */
    for (j=1;j<n_unknown;j++) {
      gam[j] = c(isolve,j-1)/bet;
      bet=b(isolve,j)-a(isolve,j)*gam[j];
      //#if DEBUG
      if (bet==0.0) 
        throw("Error 2 in tridiag_solver");
      //#endif    
      lhs(isolve,j)=(rhs(isolve,j)-a(isolve,j)*lhs(isolve,j-1))/bet;
    }
    /* backsubstitution */
    for (j=(n_unknown-2);j>=0;j--){
      lhs(isolve,j) -= gam[j+1]*lhs(isolve,j+1);
    }
  }
}
