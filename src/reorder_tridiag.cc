
#include "tridiag.h" /* for type defs */

void reorder_tridiag(ArrayNd<FTYPE,1> &a, ArrayNd<FTYPE,1> &b,
		     ArrayNd<FTYPE,1> &c, ArrayNd<FTYPE,1> &rhs,
		     ArrayNd<FTYPE,1> &v_lk, ArrayNd<FTYPE,1> &v_uk, 
		     FTYPE &rhs_lk, FTYPE &rhs_uk)
{
  /*
    Per Moulton Paper:

    v_lk(0) -> v_lk,jk
    v_lk(1) -> v_lk,i
    v_lk(2) -> v_lk,i+1

  */

  int n_local_rows = a.size(0);
  
  /* build lower row, starting from 2nd local row */
  v_lk(0) = a(1); v_lk(1) = b(1); v_lk(2) = c(1); rhs_lk=rhs(1);

  FTYPE alpha = 0;
  for (int irow=2;irow<n_local_rows; irow++) {
    alpha = a(irow)/v_lk(1); /* = a(irow)/b(irow-1); */
    v_lk(0) = -alpha*v_lk(0); 
    v_lk(1) = b(irow)-alpha*v_lk(2);
    v_lk(2) = c(irow); 
    rhs_lk = rhs(irow)-alpha*rhs_lk;
  }

  /* build upper row, starting from 2nd to last local row */
  v_uk(0) = a(n_local_rows-2); v_uk(1) = b(n_local_rows-2);
  v_uk(2) = c(n_local_rows-2); rhs_uk=rhs(n_local_rows-2);

  for (int irow=n_local_rows-3;irow>=0;irow--) {
    alpha = c(irow)/v_uk(1);
    v_uk(2) = -alpha*v_uk(2);
    v_uk(1) = b(irow) - alpha*v_uk(0);
    v_uk(0) = a(irow);
    rhs_uk = rhs(irow)-alpha*rhs_uk;
  }
  
}
