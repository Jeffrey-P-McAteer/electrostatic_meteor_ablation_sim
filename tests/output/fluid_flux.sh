#!/usr/bin/bash

source ../test_routines.sh

echo "Test for fluid flux output"
printf "Running EPPIC ......... "
../../src/eppic -procs 4 fluid_flux.i > /dev/null 2>&1
test_result
printf "Copying input file .... "
cp fluid_flux.i eppic.i
test_result
printf "Running IDL ........... "
idl fluid_flux.pro > /dev/null 2>&1
test_result
exit 0

