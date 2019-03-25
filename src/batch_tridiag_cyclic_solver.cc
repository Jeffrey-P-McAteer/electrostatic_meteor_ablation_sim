

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

void batch_tridiag_cyclic_solver(ArrayNd<FTYPE,2> &a,ArrayNd<FTYPE,2> &b, 
				 ArrayNd<FTYPE,2> &c,
				 ArrayNd<FTYPE,2> &rhs, ArrayNd<FTYPE,2> &lhs,
				 ArrayNd<FTYPE,1> &alpha, 
				 ArrayNd<FTYPE,1> &beta) 
{
  
  int n_solve = rhs.size(0);
  int n_unknown = rhs.size(1);
  int i;
  FTYPE fact,gamma;
  if (n_unknown <= 2 ) 
    throw("n_unknown too small in cyclic");
  ArrayNd<FTYPE,2> bb=ArrayNd<FTYPE,2>(n_solve,n_unknown);
  ArrayNd<FTYPE,2> u=ArrayNd<FTYPE,2>(n_solve,n_unknown);
  ArrayNd<FTYPE,2> z=ArrayNd<FTYPE,2>(n_solve,n_unknown);
  
  for (int isolve = 0; isolve<n_solve;isolve++) {
    gamma = -b(isolve,0);
    bb(isolve,0) = b(isolve,0)-gamma;
    bb(isolve,n_unknown-1) = b(isolve,n_unknown-1)
      -alpha(isolve)*beta(isolve)/gamma;
    u(isolve,0)=gamma;
    u(isolve,n_unknown-1)=alpha(isolve);
    for (i=1;i<n_unknown-1;i++) bb(isolve,i)=b(isolve,i);
  }

  batch_tridiag_solver(a,bb,c,rhs,lhs);

  batch_tridiag_solver(a,bb,c,u,z);

  for (int isolve = 0; isolve<n_solve;isolve++) {
    gamma = -b(isolve,0); 
    fact=(lhs(isolve,0)+beta(isolve)*lhs(isolve,n_unknown-1)/
	  gamma)/(1.0+z(isolve,0)+beta(isolve)
		  *z(isolve,n_unknown-1)/gamma);
    for (i=0;i<n_unknown;i++) lhs(isolve,i) -= fact*z(isolve,i);
  }
  
}
