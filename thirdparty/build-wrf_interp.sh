#!/bin/bash

# Tested in Ubuntu 16.04

: ${WRF_ROOT:=${HOME}/thirdparty}

mkdir -p $WRF_ROOT/wrf_interp
cd $WRF_ROOT/wrf_interp

test -f WRF_INTERP.TAR.gz || wget http://www2.mmm.ucar.edu/wrf/src/WRF_INTERP.TAR.gz
tar xvf WRF_INTERP.TAR.gz wrf_interp.F90

# Little fixes + change so that the namelist is read from stdin instead of a
# file namelist.input...

sed -i wrf_interp.F90 \
    -e 's/\.eq\. \.T/.eqv. .T/' \
    -e 's/\(call handle_err(rcode\)[^)]*/\1/' \
    -e '2550s/REAL/INTEGER/' \
    -e '858,866d;867i  READ(*,io)\n  READ(*,interp_in)' \
    -e '2439s/time_ivar/times_id/' \
    -e '2339a  integer :: times_id\n  rcode=nf_inq_varid(ncid,"Times",times_id)'

gfortran -o wrf_interp.exe wrf_interp.F90 -I/usr/include -free -L/usr/lib -lnetcdff

