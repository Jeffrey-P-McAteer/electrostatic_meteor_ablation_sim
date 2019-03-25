#!/bin/csh
#
# Execute this script in the fftw-2.1 directory

./configure --enable-type-prefix --prefix=$HOME/ebeam/eppic
make
make install
make clean
./configure --enable-float --enable-type-prefix --prefix=$HOME/ebeam/eppic
make
make install
make clean

