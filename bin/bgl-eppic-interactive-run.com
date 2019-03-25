#!/bin/csh
mpirun -noallocate -partition R001 -np 1024 -mode VN -cwd `pwd` eppic.x >& eppic.log &
