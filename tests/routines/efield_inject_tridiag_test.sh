#!/bin/bash
# A test of efield subroutine.
# While the test is all taken care of in the .cc executable, this script is
# setup to handle mpi running. Further tests can be added at a later time, for 
# example tests involving different charge densities. 

echo "....... TESTING efield_inject_tridiag_test.sh"
echo " This sould be printed!!!"
source run_mpi.tests

nproc=1
echo "....... Single processor open job"
echo run_mpi ./efield_inject_tridiag_test -1 -1 
run_mpi ./efield_inject_tridiag_test -1 -1 || exit 1
echo "....... Single processor dirichlet job"
run_mpi ./efield_inject_tridiag_test 0 0|| exit 1
#echo "....... Single processor neumann"
#./efield_inject_tridiag_test 1 1|| exit 1
echo "....... Single processor periodic job"
run_mpi ./efield_inject_tridiag_test 2 2|| exit 1


nproc=2
echo "....... Multi processor open job, using" $nproc "processors"
run_mpi "./efield_inject_tridiag_test -1 -1 " || exit 1
echo "....... Multi processor dirichlet job, using" $nproc "processors"
run_mpi "./efield_inject_tridiag_test 0 0" || exit 1
#echo "....... Multi processor neumann"
#run_mpi "./efield_inject_tridiag_test 1 1" || exit 1
echo "....... Multi processor periodic job, using" $nproc "processors"
run_mpi "./efield_inject_tridiag_test 2 2" || exit 1

nproc=4

echo "....... Multi processor open job, using" $nproc "processors"
run_mpi "./efield_inject_tridiag_test -1 -1 " || exit 1
echo "....... Multi processor dirichlet job, using" $nproc "processors"
run_mpi "./efield_inject_tridiag_test 0 0" || exit 1
#echo "....... Multi processor neumann"
#run_mpi "./efield_inject_tridiag_test 1 1" || exit 1
echo "....... Multi processor periodic job, using" $nproc "processors"
run_mpi "./efield_inject_tridiag_test 2 2" || exit 1

nproc=8

echo "....... Multi processor open job, using" $nproc "processors"
run_mpi "./efield_inject_tridiag_test -1 -1 " || exit 1
echo "....... Multi processor dirichlet job, using" $nproc "processors"
run_mpi "./efield_inject_tridiag_test 0 0" || exit 1
#echo "....... Multi processor neumann"
#run_mpi "./efield_inject_tridiag_test 1 1" || exit 1
echo "....... Multi processor periodic job, using" $nproc "processors"
run_mpi "./efield_inject_tridiag_test 2 2" || exit 1

exit 0
