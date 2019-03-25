
#include "tridiag.h" /* for serial solver, and reorder helper function */

/* 
   This algorithm (and it's notation) is taken from Austin, Berndt and 
   Moulton, "A Memory Efficient Parallel Tridiagonal Solver", published online
   through a LANL website.

   It uses a cyclic tridiagonal solver for the one global solver. Otherwise
   it is analogous to the tridiag_solver_parallel routine. 
*/

void tridiag_cyclic_solver_parallel(ArrayNd<FTYPE,1> &a,ArrayNd<FTYPE,1> &b, 
				    ArrayNd<FTYPE,1> &c,
				    ArrayNd<FTYPE,1> &rhs, 
				    ArrayNd<FTYPE,1> &lhs,
				    FTYPE alpha, FTYPE beta,
				    MPI_Comm comm) 
{
  
  int n_unknown = rhs.size(0);

  int mpi_rank,mpi_np;
  MPI_Comm_rank(comm,&mpi_rank);
  MPI_Comm_size(comm,&mpi_np);

  if (mpi_np==1) {
    tridiag_cyclic_solver(a,b,c,rhs,lhs,alpha,beta);
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

  ArrayNd<FTYPE,1> v_lk,v_uk;
  FTYPE rhs_lk=0,rhs_uk=0;
  v_lk=v_uk=ArrayNd<FTYPE,1>(3)=0;
  
  /* "phase 1" */
  /* reorder local tridiagonal */
  
  reorder_tridiag(a,b,c,rhs,v_lk,v_uk,rhs_lk,rhs_uk);

  /* "phase 2" */
  /* communicate border rows */
  /* tridiag info for communication rows */
  ArrayNd<FTYPE,1> local_a_comm_rows,local_b_comm_rows,local_c_comm_rows,
    local_rhs_comm_rows;
  ArrayNd<FTYPE,1> a_comm_rows,b_comm_rows,c_comm_rows,
    rhs_comm_rows,lhs_comm_rows;

  //cout << "rank " << mpi_rank << " before gather\n";
  
  if (mpi_rank == 0) {
    a_comm_rows = ArrayNd<FTYPE,1>(2*mpi_np);
    b_comm_rows = ArrayNd<FTYPE,1>(2*mpi_np);
    c_comm_rows = ArrayNd<FTYPE,1>(2*mpi_np);
    rhs_comm_rows = ArrayNd<FTYPE,1>(2*mpi_np);
    lhs_comm_rows = ArrayNd<FTYPE,1>(2*mpi_np);
  }
  local_a_comm_rows = ArrayNd<FTYPE,1>(2); 
  local_b_comm_rows = ArrayNd<FTYPE,1>(2); 
  local_c_comm_rows = ArrayNd<FTYPE,1>(2); 
  local_rhs_comm_rows = ArrayNd<FTYPE,1>(2); 

  local_a_comm_rows[0] = v_uk[0];
  local_a_comm_rows[1] = v_lk[0];
  MPI_Gather(local_a_comm_rows.address(0),2,MPI_FTYPE,a_comm_rows.address(0)
	     ,2,MPI_FTYPE,0,comm);

  local_b_comm_rows[0] = v_uk[1];
  local_b_comm_rows[1] = v_lk[1];
  MPI_Gather(local_b_comm_rows.address(0),2,MPI_FTYPE,b_comm_rows.address(0)
	     ,2,MPI_FTYPE,0,comm);

  local_c_comm_rows[0] = v_uk[2];
  local_c_comm_rows[1] = v_lk[2];
  MPI_Gather(local_c_comm_rows.address(0),2,MPI_FTYPE,c_comm_rows.address(0)
	     ,2,MPI_FTYPE,0,comm);

  local_rhs_comm_rows[0] = rhs_uk;
  local_rhs_comm_rows[1] = rhs_lk;
  MPI_Gather(local_rhs_comm_rows.address(0),2,MPI_FTYPE,
	     rhs_comm_rows.address(0),2,MPI_FTYPE,0,comm);
  
  /* solve border problem */
  //cout << "rank " << mpi_rank << " before common solver\n";
  if (mpi_rank==0) {
    tridiag_cyclic_solver(a_comm_rows,b_comm_rows,c_comm_rows,
			  rhs_comm_rows,lhs_comm_rows,alpha,beta);
    //    cout << "Solution to common tridiag\n";
    int ncommon = 2*mpi_np;
    for (int i = 0; i<ncommon;i++) {
      /*
      cout << i 
	   << ", a= " << a_comm_rows[i] 
	   << ", b= " << b_comm_rows[i] 
	   << ", c= " << c_comm_rows[i] 
	   << ", lhs= " << lhs_comm_rows[i] << endl;
      */
    }
  }

  /* communicate border solution */
      //  cout << "rank " << mpi_rank << " before scatter\n";
  MPI_Scatter(lhs_comm_rows.address(0),2,MPI_FTYPE,
	      local_rhs_comm_rows.address(0),2,MPI_FTYPE,0,comm);

  /*
    cout << "rank " << mpi_rank << " local_rhs_comm_rows " 
       << local_rhs_comm_rows[0] << ", "
       << local_rhs_comm_rows[1] << endl;
  */
  FTYPE a0_old = a[0];
  FTYPE an_unknown_old = a[n_unknown-1];
  a[0]=0;
  a[n_unknown-1] = 0;
  FTYPE c0_old = c[0];
  FTYPE cn_unknown_old = c[n_unknown-1];
  c[0]=0;
  c[n_unknown-1] = 0;
  FTYPE b0_old = b[0];
  FTYPE bn_unknown_old = b[n_unknown-1];
  b[0]=1;
  b[n_unknown-1] = 1;
  FTYPE rhs0_old = rhs[0];
  FTYPE rhsn_unknown_old = rhs[n_unknown-1];
  /* should be equation to local_lhs_comm_rows
     but only local_rhs_comm_rows exists and is reused, see above */
  rhs[0]=local_rhs_comm_rows[0];
  rhs[n_unknown-1] = local_rhs_comm_rows[1];

  /* "phase 3" */
  /* solve local tridiagonal */
  //  cout << "rank " << mpi_rank << " before local solver\n";
  tridiag_solver(a,b,c,rhs,lhs);
  
  /* reset variables that should not have been changed */
  a[0] = a0_old;
  a[n_unknown-1] = an_unknown_old;
  b[0] = b0_old;
  b[n_unknown-1] = bn_unknown_old;
  c[0] = c0_old;
  c[n_unknown-1] = cn_unknown_old;
  rhs[0] = rhs0_old;
  rhs[n_unknown-1] = rhsn_unknown_old;
  
}
