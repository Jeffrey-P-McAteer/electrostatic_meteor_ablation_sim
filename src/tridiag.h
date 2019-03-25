#ifndef TRIDIAG_H
#define TRIDIAG_H
#include "eppic-types.h"
#include <complex>

void tridiag_solver(ArrayNd<FTYPE,1> &a,ArrayNd<FTYPE,1> &b, 
		    ArrayNd<FTYPE,1> &c,
		    ArrayNd<FTYPE,1>&rhs, ArrayNd<FTYPE,1> &lhs);

// Assumes a=c=1, b=constant
void tridiag_solver_complex(FTYPE b,
                            std::complex<FTYPE> *rhs,
                            std::complex<FTYPE> *lhs,
                            int n_unknown);

void tridiag_solver_complex(std::complex<FTYPE> *a,
                            std::complex<FTYPE> *b,
                            std::complex<FTYPE> *c,
                            std::complex<FTYPE> *rhs,
                            std::complex<FTYPE> *lhs,
                            int n_unknown);

void tridiag_solver_complex_neumann_dirichlet(FTYPE boundary_condition_low,
                                              FTYPE boundary_condition_high,
                                              FTYPE b,
                                              std::complex<FTYPE> *rhs,
                                              std::complex<FTYPE> *lhs,
                                              int n_unknown);

void tridiag_solver_complex_dirichlet_neumann(FTYPE boundary_condition_low,
                                              FTYPE boundary_condition_high,
                                              FTYPE b,
                                              std::complex<FTYPE> *rhs,
                                              std::complex<FTYPE> *lhs,
                                              int n_unknown);

void tridiag_solver_complex_dirichlet_dirichlet(FTYPE boundary_condition_low,
                                                FTYPE boundary_condition_high,
                                                FTYPE b,
                                                std::complex<FTYPE> *rhs,
                                                std::complex<FTYPE> *lhs,
                                                int n_unknown);

void tridiag_solver_complex_neumann_neumann(FTYPE boundary_condition_low,
                                            FTYPE boundary_condition_high,
                                            FTYPE b,
                                            std::complex<FTYPE> *rhs,
                                            std::complex<FTYPE> *lhs,
                                            int n_unknown);

void tridiag_solver_complex_open_open(FTYPE b, FTYPE rinv,
                                      std::complex<FTYPE> *rhs,
                                      std::complex<FTYPE> *lhs,
                                      int n_unknown);

void batch_tridiag_solver(ArrayNd<FTYPE,2> &a,ArrayNd<FTYPE,2> &b, 
			  ArrayNd<FTYPE,2> &c,
			  ArrayNd<FTYPE,2>&rhs, ArrayNd<FTYPE,2> &lhs);


void tridiag_cyclic_solver(ArrayNd<FTYPE,1> &a,ArrayNd<FTYPE,1> &b, 
			   ArrayNd<FTYPE,1> &c,
			   ArrayNd<FTYPE,1>&rhs, ArrayNd<FTYPE,1> &lhs,
			   FTYPE alpha, FTYPE beta);

void tridiag_cyclic_solver_complex(FTYPE b,
                                   std::complex<FTYPE> *rhs,
                                   std::complex<FTYPE> *lhs,
                                   const FTYPE alpha, 
                                   const FTYPE beta,
                                   int n_unknown);

void tridiag_cyclic_solver_complex(std::complex<FTYPE> *a,
                                   std::complex<FTYPE> *b,
                                   std::complex<FTYPE> *c,
                                   std::complex<FTYPE> *rhs,
                                   std::complex<FTYPE> *lhs,
                                   const FTYPE alpha, 
                                   const FTYPE beta,
                                   int n_unknown);

void batch_tridiag_cyclic_solver(ArrayNd<FTYPE,2> &a,ArrayNd<FTYPE,2> &b, 
				 ArrayNd<FTYPE,2> &c,
				 ArrayNd<FTYPE,2> &rhs, ArrayNd<FTYPE,2> &lhs,
				 ArrayNd<FTYPE,1> &alpha, 
				 ArrayNd<FTYPE,1> &beta); 

void reorder_tridiag(ArrayNd<FTYPE,1> &a, ArrayNd<FTYPE,1> &b,
		     ArrayNd<FTYPE,1> &c, ArrayNd<FTYPE,1> &rhs,
		     ArrayNd<FTYPE,1> &v_lk, ArrayNd<FTYPE,1> &v_uk, 
		     FTYPE &rhs_lk, FTYPE &rhs_uk);

void batch_reorder_tridiag(ArrayNd<FTYPE,2> &a, ArrayNd<FTYPE,2> &b,
		     ArrayNd<FTYPE,2> &c, ArrayNd<FTYPE,2> &rhs,
		     ArrayNd<FTYPE,2> &v_lk, ArrayNd<FTYPE,2> &v_uk, 
		     ArrayNd<FTYPE,1> &rhs_lk, ArrayNd<FTYPE,1> &rhs_uk);


#ifdef USE_MPI
#include "eppic-mpi.h"

void tridiag_solver_parallel(ArrayNd<FTYPE,1> &a,ArrayNd<FTYPE,1> &b, 
			     ArrayNd<FTYPE,1> &c,
			     ArrayNd<FTYPE,1> &rhs, ArrayNd<FTYPE,1> &lhs,
			     MPI_Comm comm);

void batch_tridiag_solver_parallel(ArrayNd<FTYPE,2> &a,ArrayNd<FTYPE,2> &b, 
				   ArrayNd<FTYPE,2> &c,
				   ArrayNd<FTYPE,2> &rhs, ArrayNd<FTYPE,2> &lhs,
				   MPI_Comm comm);


void tridiag_cyclic_solver_parallel(ArrayNd<FTYPE,1> &a,ArrayNd<FTYPE,1> &b, 
				    ArrayNd<FTYPE,1> &c,
				    ArrayNd<FTYPE,1> &rhs, 
				    ArrayNd<FTYPE,1> &lhs,
				    FTYPE alpha, FTYPE beta,
				    MPI_Comm comm);

void batch_tridiag_cyclic_solver_parallel(ArrayNd<FTYPE,2> &a,ArrayNd<FTYPE,2> &b, 
					  ArrayNd<FTYPE,2> &c,
					  ArrayNd<FTYPE,2> &rhs, ArrayNd<FTYPE,2> &lhs,
					  ArrayNd<FTYPE,1> alpha, ArrayNd<FTYPE,1> beta,
					  MPI_Comm comm);


void tridiag_solver_parallel_recursive(ArrayNd<FTYPE,1> &a,
				       ArrayNd<FTYPE,1> &b, 
				       ArrayNd<FTYPE,1> &c,
				       ArrayNd<FTYPE,1> &rhs, 
				       ArrayNd<FTYPE,1> &lhs,
				       MPI_Comm comm,
				       int max_np);

#endif
#endif
