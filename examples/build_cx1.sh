#!/bin/sh 
module load gsl
module load intel-suite
export LIBRARY_PATH=$LD_LIBRARY_PATH
cd src
make -f Makefile.cx1 2> /dev/null
cd ..
