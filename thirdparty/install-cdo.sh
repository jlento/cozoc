#!/bin/bash

set -e

# Tested in Ubuntu 16.04

: ${CDO_ROOT:=${HOME}/thirdparty}

test -f cdo-1.9.4.tar.gz || wget https://code.mpimet.mpg.de/attachments/download/17374/cdo-1.9.4.tar.gz
test -f cdo-1.9.4 || tar xvf cdo-1.9.4.tar.gz

cd cdo-1.9.4
./configure --prefix=${CDO_ROOT} --with-netcdf
make
make install
