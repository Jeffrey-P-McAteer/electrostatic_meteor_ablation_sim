
#include "tridiag.h" /* for type defs */

void batch_reorder_tridiag(ArrayNd<FTYPE,2> &a, ArrayNd<FTYPE,2> &b,
		     ArrayNd<FTYPE,2> &c, ArrayNd<FTYPE,2> &rhs,
		     ArrayNd<FTYPE,2> &v_lk, ArrayNd<FTYPE,2> &v_uk, 
		     ArrayNd<FTYPE,1> &rhs_lk, ArrayNd<FTYPE,1> &rhs_uk)
{
  /*
    Per Moulton Paper:

    v_lk(0) -> v_lk,jk
    v_lk(1) -> v_lk,i
    v_lk(2) -> v_lk,i+1

  */

  int n_solve = a.size(0);
  int n_local_rows = a.size(1);
  
  for (int isolve=0;isolve<n_solve;isolve++) {
    /* build lower row, starting from 2nd local row */
    v_lk(isolve,0) = a(isolve,1); 
    v_lk(isolve,1) = b(isolve,1); 
    v_lk(isolve,2) = c(isolve,1); 
    rhs_lk(isolve)=rhs(isolve,1);
    
    FTYPE alpha = 0;
    for (int irow=2;irow<n_local_rows; irow++) {
      alpha = a(isolve,irow)/v_lk(isolve,1); /* = a(isolve,irow)/b(isolve,irow-1); */
      v_lk(isolve,0) = -alpha*v_lk(isolve,0); 
      v_lk(isolve,1) = b(isolve,irow)-alpha*v_lk(isolve,2);
      v_lk(isolve,2) = c(isolve,irow); 
      rhs_lk(isolve) = rhs(isolve,irow)-alpha*rhs_lk(isolve);
    }
    
    /* build upper row, starting from 2nd to last local row */
    v_uk(isolve,0) = a(isolve,n_local_rows-2); 
    v_uk(isolve,1) = b(isolve,n_local_rows-2);
    v_uk(isolve,2) = c(isolve,n_local_rows-2); 
    rhs_uk(isolve)=rhs(isolve,n_local_rows-2);
    
    for (int irow=n_local_rows-3;irow>=0;irow--) {
      alpha = c(isolve,irow)/v_uk(isolve,1);
      v_uk(isolve,2) = -alpha*v_uk(isolve,2);
      v_uk(isolve,1) = b(isolve,irow) - alpha*v_uk(isolve,0);
      v_uk(isolve,0) = a(isolve,irow);
      rhs_uk(isolve) = rhs(isolve,irow)-alpha*rhs_uk(isolve);
    }
  }
  
}
