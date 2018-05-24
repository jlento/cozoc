#!/bin/bash

# Tested in Ubuntu 16.04

: ${WRF_ROOT:=${HOME}/thirdparty}

# OS X workaround
readlink-f () {
    local absname
    readlink "$1" || ( test -d "$1" && cd "$1" && pwd -P ) || (
        cd $(dirname "$1") 2> /dev/null && echo "$(pwd -P)/$(basename "$1")"
    )
}

srcdir="$(dirname $(readlink-f ${BASH_SOURCE[0]}))"

cd $WRF_ROOT/WRFV3/test/em_b_wave

test -f iofield_list.txt || cp ${srcdir}/iofield_list.txt .
test -f namelist.input.orig || cp namelist.input namelist.input.orig
cp ${srcdir}/namelist.input .

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
./run_me_first.csh
./ideal.exe
./wrf.exe
