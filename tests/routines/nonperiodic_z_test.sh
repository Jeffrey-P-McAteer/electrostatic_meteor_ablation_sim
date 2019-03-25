#!/bin/bash

echo "Compiling the test..."
g++ -std=c++11 -O1 -g -D NDIM=3 -I../../src -I../../src/classes -I../../include -c ../../src/tridiag_solver.cc -o tridiag_solver.o
g++ -std=c++11 -O1 -g -I../../src -I../../src/classes -I../../include -D NDIM=3 tridiag_solver.o nonperiodic_z_test.cc -o test.x
echo "Running the test..."
if ./test.x
then
    echo "Test passed."
else
    echo "Test failed."
fi
