#!/bin/bash

set -ex

# Tested in Ubuntu 16.04


# Tools
NCKS=$(which ncks)
: ${NCKS:?No ncks in PATH}

which ${WRF_INTERP} || \
    : ${WRF_INTERP:=$HOME/thirdparty/wrf_interp/wrf_interp.exe}

nc-select-steps () {
    first=$1
    last=$2
    infile=$3
    outfile=$4
    ncks -O -a -h -d Time,${first},${last} $infile $outfile
}

wrf-interpolate () {
    dp=$1
    date=$2
    infile=$3
    outfile=$4
    ${WRF-INTERP} <<"EOF"
&io
 path_to_input = '$(dirname $2)'
 path_to_output = '$(dirname $3)'
 root_name = '${2%%_d01_*}'
 grid_id = 1
 start_date = '${date}'
 leap_year  = .FALSE.
 debug = .TRUE.
/
&interp_in
  interp_levels = 1000,-100,${dp}
  extrapolate = 1
  unstagger_grid = .TRUE.
  vert_coordinate = 'pres'
/
EOF
}
