#!/bin/bash

# Tested in Ubuntu 16.04

: ${WRF_ROOT:=${HOME}/thirdparty}

cd $WRF_ROOT

test -f WRFV3.9.1.1.TAR.gz || wget http://www2.mmm.ucar.edu/wrf/src/WRFV3.9.1.1.TAR.gz
test -d WRFV3 || tar xvf WRFV3.9.1.1.TAR.gz

cd WRFV3
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
./configure <<<'
/usr/include
/usr/lib/x86_64-linux-gnu
32

'

./compile em_b_wave


