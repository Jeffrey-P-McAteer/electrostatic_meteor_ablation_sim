#!/bin/bash
# A test of efield subroutine.
# While the test is all taken care of in the .cc executable, this script is
# setup to handle mpi running. Further tests can be added at a later time, for 
# example tests involving different charge densities. 

echo "....... TESTING tridiag_solver_parallel_test.sh"

source run_mpi.tests

echo "....... Single processor"
nproc=1
run_mpi ./tridiag_solver_parallel_test || echo exit 1
echo "....... Multi processor"
nproc=2
run_mpi ./tridiag_solver_parallel_test || exit 1
echo "....... Multi processor (4)"
nproc=4
run_mpi ./tridiag_solver_parallel_test || exit 1
exit 0
