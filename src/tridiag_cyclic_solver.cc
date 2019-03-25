

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
   alpha -- FTYPE -- value of lower, left corner of tridiag matrix
   beta -- FTYPE -- value of upper, right corner of tridiag matrix

   Note: all input arryas must be predefined and have the same dimensions. 
   There is no test done to verify this.

   only notation differences from Press code:
   (this code) -> (Press code)
   n_unknown -> n
   rhs -> r
   lhs -> x
*/

void tridiag_cyclic_solver(ArrayNd<FTYPE,1> &a,ArrayNd<FTYPE,1> &b, 
			   ArrayNd<FTYPE,1> &c,
			   ArrayNd<FTYPE,1>&rhs, ArrayNd<FTYPE,1> &lhs,
			   FTYPE alpha, FTYPE beta) 
{

  if ((alpha==0)&&(beta==0)) {
    // not cyclic!!
    tridiag_solver(a,b,c,rhs,lhs);
    return;
  }
  int n_unknown = rhs.size(0);
  int i;
  FTYPE fact,gamma;
  if (n_unknown <= 2 ) 
    throw("n_unknown too small in cyclic");
  FTYPEAVec bb=FTYPEAVec(n_unknown)=0;  
  FTYPEAVec u=FTYPEAVec(n_unknown)=0;  
  FTYPEAVec z=FTYPEAVec(n_unknown)=0;  
  gamma = -b[0];
  bb[0] = b[0]-gamma;
  bb[n_unknown-1] = b[n_unknown-1]-alpha*beta/gamma;
  for (i=1;i<n_unknown-1;i++) bb[i]=b[i];
  tridiag_solver(a,bb,c,rhs,lhs);
  u[0]=gamma;
  u[n_unknown-1]=alpha;
  tridiag_solver(a,bb,c,u,z);
  fact=(lhs[0]+beta*lhs[n_unknown-1]/gamma)/
    (1.0+z[0]+beta*z[n_unknown-1]/gamma);
  for (i=0;i<n_unknown;i++) lhs[i] -= fact*z[i];
  
}



// ASSUMES a=c=1, b is a constant
void tridiag_cyclic_solver_complex(FTYPE b,
                                   std::complex<FTYPE> *rhs,
                                   std::complex<FTYPE> *lhs,
                                   const FTYPE alpha, 
                                   const FTYPE beta,
                                   int n_unknown) 
{
  if ((alpha==0)&&(beta==0)) {
    // not cyclic!!
    //tridiag_solver(a,b,c,rhs,lhs);
    tridiag_solver_complex(b,rhs,lhs,n_unknown);
    return;
  }
  int i;
  std::complex<FTYPE> fact,gamma;
  if (n_unknown <= 2 )
    throw("n_unknown too small in cyclic");
  //std::complex<FTYPE> bb[n_unknown]=0;
  std::complex<FTYPE> bb[n_unknown];
  //std::complex<FTYPE> u[n_unknown]=0;
  std::complex<FTYPE> u[n_unknown];
  //std::complex<FTYPE> z[n_unknown]=0;
  std::complex<FTYPE> z[n_unknown];
  //std::complex<FTYPE> a[n_unknown]=1;
  std::complex<FTYPE> a[n_unknown];
  //std::complex<FTYPE> c[n_unknown]=1;
  std::complex<FTYPE> c[n_unknown];
  //gamma = -b[0];
  gamma = -b;
  //bb[0] = b[0]-gamma;
  bb[0] = b-gamma;
  //bb[n_unknown-1] = b[n_unknown-1]-alpha*beta/gamma;
  bb[n_unknown-1] = b-alpha*beta/gamma;
  //for (i=1;i<n_unknown-1;i++) bb[i]=b[i];
  for (i=1;i<n_unknown-1;i++) bb[i]=b;
  for (i=0;i<n_unknown;i++){
    u[i] = 0;
    z[i] = 0;
    a[i] = 1;
    c[i] = 1;
  }
  tridiag_solver_complex(a,bb,c,rhs,lhs,n_unknown);
  u[0]=gamma;
  u[n_unknown-1]=alpha;
  tridiag_solver_complex(a,bb,c,u,z,n_unknown);
  fact=(lhs[0]+beta*lhs[n_unknown-1]/gamma)/
    (1.0+z[0]+beta*z[n_unknown-1]/gamma);
  for (i=0;i<n_unknown;i++) lhs[i] -= fact*z[i];
}
