#!/bin/bash

# Tested in Ubuntu 16.04

: ${WRF_ROOT:=${HOME}/thirdparty}

mkdir -p $WRF_ROOT/wrf_interp
cd $WRF_ROOT/wrf_interp

test -f WRF_INTERP.TAR.gz || wget http://www2.mmm.ucar.edu/wrf/src/WRF_INTERP.TAR.gz
test -f /wrf_interp || tar xvf WRF_INTERP.TAR.gz

sed -i 's/\.eq\. \.T/.eqv. .T/' wrf_interp.F90
sed -i 's/\(call handle_err(rcode\)[^)]*/\1/'  wrf_interp.F90
sed -i '2550s/REAL/INTEGER/' wrf_interp.F90

gfortran -o wrf_interp.exe wrf_interp.F90 -I/usr/include -free -L/usr/lib -lnetcdff

