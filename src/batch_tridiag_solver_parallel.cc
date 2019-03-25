
#include "tridiag.h" // for serial solver, and reorder helper function 
#include "eppic-math.h" // for transpose2d routine
/* 
   This algorithm (and it's notation) is taken from Austin, Berndt and 
   Moulton, "A Memory Efficient Parallel Tridiagonal Solver", published online
   through a LANL website.
*/

void batch_tridiag_solver_parallel(ArrayNd<FTYPE,2> &a,ArrayNd<FTYPE,2> &b, 
				   ArrayNd<FTYPE,2> &c,
				   ArrayNd<FTYPE,2> &rhs, ArrayNd<FTYPE,2> &lhs,
				   MPI_Comm comm) 
{
  
  int n_solve = rhs.size(0);
  int n_unknown = rhs.size(1);

  int mpi_rank,mpi_np;
  MPI_Comm_rank(comm,&mpi_rank);
  MPI_Comm_size(comm,&mpi_np);
  if (mpi_np==1) {
    batch_tridiag_solver(a,b,c,rhs,lhs);
    return;
  }

  /* 
     These symbols (variables) are the same as those in the Moulton paper
     The rows of the tridiagonal are represented by a vector v_i in the 
     paper, but here each of the three terms that each row has are placed 
     in the three vectors a, b and c at their ith element. 

     The one exception is the vector b, which Moulton uses to represent the
     rhs of the Ax=b equation. So b_uk,b_lk in his notation correspond to 
     rhs_uk,rhs_lk in my notation. 
     
     The index i is a local index. 

     The index j_k is the lowest global row number on processor k.

     The subscripts lk and uk refere to lowest row (ie i=Nk-1) 
     and highest row on processor k (ie i=0). 
  */

  ArrayNd<FTYPE,2> v_lk,v_uk;
  ArrayNd<FTYPE,1> rhs_lk,rhs_uk;
  v_lk=v_uk=ArrayNd<FTYPE,2>(n_solve,3)=0;
  rhs_lk=rhs_uk=ArrayNd<FTYPE,1>(n_solve)=0;
  
  /* "phase 1" */
  /* reorder local tridiagonal */
  
  batch_reorder_tridiag(a,b,c,rhs,v_lk,v_uk,rhs_lk,rhs_uk);

  /* "phase 2" */
  /* communicate border rows */
  /* tridiag info for communication rows */
  ArrayNd<FTYPE,2> a_comm_rows,b_comm_rows,c_comm_rows,
    rhs_comm_rows,lhs_comm_rows, comm_rows_buffer;

  // initially setup in transposed form, for efficient MPI communication
  a_comm_rows = ArrayNd<FTYPE,2>(2*mpi_np,n_solve);
  comm_rows_buffer = a_comm_rows;
  b_comm_rows = ArrayNd<FTYPE,2>(2*mpi_np,n_solve);
  c_comm_rows = ArrayNd<FTYPE,2>(2*mpi_np,n_solve);
  rhs_comm_rows = ArrayNd<FTYPE,2>(2*mpi_np,n_solve);
  // no need to put in transposed form:
  lhs_comm_rows = ArrayNd<FTYPE,2>(n_solve,2*mpi_np);
  
  for (int isolve=0;isolve<n_solve;isolve++) {
    a_comm_rows(mpi_rank*2,isolve) = v_uk(isolve,0);
    a_comm_rows(mpi_rank*2+1,isolve) = v_lk(isolve,0);
    b_comm_rows(mpi_rank*2,isolve) = v_uk(isolve,1);
    b_comm_rows(mpi_rank*2+1,isolve) = v_lk(isolve,1);
    c_comm_rows(mpi_rank*2,isolve) = v_uk(isolve,2);
    c_comm_rows(mpi_rank*2+1,isolve) = v_lk(isolve,2);
    rhs_comm_rows(mpi_rank*2,isolve) = rhs_uk(isolve);
    rhs_comm_rows(mpi_rank*2+1,isolve) = rhs_lk(isolve);
  }


  int start_idx[2]={mpi_rank*2,0};
  MPI_Allgather(a_comm_rows.address(start_idx),2*n_solve,MPI_FTYPE,
		comm_rows_buffer.address(0),2*n_solve,MPI_FTYPE,comm);
  a_comm_rows = comm_rows_buffer;
  transpose2d(a_comm_rows);
  MPI_Allgather(b_comm_rows.address(start_idx),2*n_solve,MPI_FTYPE,
		comm_rows_buffer.address(0),2*n_solve,MPI_FTYPE,comm);
  b_comm_rows = comm_rows_buffer;
  transpose2d(b_comm_rows);
  MPI_Allgather(c_comm_rows.address(start_idx),2*n_solve,MPI_FTYPE,
		comm_rows_buffer.address(0),2*n_solve,MPI_FTYPE,comm);
  c_comm_rows = comm_rows_buffer;
  transpose2d(c_comm_rows);
  MPI_Allgather(rhs_comm_rows.address(start_idx),2*n_solve,MPI_FTYPE,
		comm_rows_buffer.address(0),2*n_solve,MPI_FTYPE,comm);
  rhs_comm_rows = comm_rows_buffer;
  transpose2d(rhs_comm_rows);

  /* solve border problem */
  batch_tridiag_solver(a_comm_rows,b_comm_rows,c_comm_rows,
		       rhs_comm_rows,lhs_comm_rows);


  // store boundary equations, which are overwritten, to replace later

  ArrayNd<FTYPE,1> a0_old = ArrayNd<FTYPE,1>(n_solve);
  ArrayNd<FTYPE,1> an_unknown_old = ArrayNd<FTYPE,1>(n_solve);
  ArrayNd<FTYPE,1> c0_old = ArrayNd<FTYPE,1>(n_solve);
  ArrayNd<FTYPE,1> cn_unknown_old = ArrayNd<FTYPE,1>(n_solve);
  ArrayNd<FTYPE,1> b0_old = ArrayNd<FTYPE,1>(n_solve);
  ArrayNd<FTYPE,1> bn_unknown_old = ArrayNd<FTYPE,1>(n_solve);
  ArrayNd<FTYPE,1> rhs0_old = ArrayNd<FTYPE,1>(n_solve);
  ArrayNd<FTYPE,1> rhsn_unknown_old = ArrayNd<FTYPE,1>(n_solve);

  for (int isolve=0;isolve<n_solve;isolve++) {
    a0_old(isolve) = a(isolve,0);
    an_unknown_old(isolve) = a(isolve,n_unknown-1);
    c0_old(isolve) = c(isolve,0);
    cn_unknown_old(isolve) = c(isolve,n_unknown-1);
    b0_old(isolve) = b(isolve,0);
    bn_unknown_old(isolve) = b(isolve,n_unknown-1);
    rhs0_old(isolve) = rhs(isolve,0);
    rhsn_unknown_old(isolve) = rhs(isolve,n_unknown-1);
    // setup first,last eqns to just copy solution found above
    b(isolve,0)=1;
    b(isolve,n_unknown-1) = 1;
    a(isolve,0)=0;
    a(isolve,n_unknown-1) = 0;
    c(isolve,0)=0;
    c(isolve,n_unknown-1) = 0;
    rhs(isolve,0)=lhs_comm_rows(isolve,mpi_rank*2);
    rhs(isolve,n_unknown-1) = lhs_comm_rows(isolve,mpi_rank*2+1);
  }

  /* "phase 3" */
  /* solve local tridiagonal */
  batch_tridiag_solver(a,b,c,rhs,lhs);
  
  /* reset variables that should not have been changed */
  for (int isolve=0;isolve<n_solve;isolve++) {
    a(isolve,0) = a0_old(isolve);
    a(isolve,n_unknown-1) = an_unknown_old(isolve);
    b(isolve,0) = b0_old(isolve);
    b(isolve,n_unknown-1) = bn_unknown_old(isolve);
    c(isolve,0) = c0_old(isolve);
    c(isolve,n_unknown-1) = cn_unknown_old(isolve);
    rhs(isolve,0) = rhs0_old(isolve);
    rhs(isolve,n_unknown-1) = rhsn_unknown_old(isolve);
  }
}

