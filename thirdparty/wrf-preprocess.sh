#!/bin/bash

set -ex

source "$(dirname $(readlink -f ${BASH_SOURCE[0]}))/wrf-tools.sh"

usage="
Usage: wrf-preprocess.sh <pressure level spacing (hPa)> <first step> <last step> <infile> <outfile>
"

first=${1:?$usage}
last=${2:?$usage}
dp=${3:?$usage}
infile=${4:?$usage}
outfile=${5:?$usage}


timeslice=$(nc-select-steps $first $last $infile)

#wrf-interpolate $dp $timeslice

wrf-fix-dimensions-etc ${timeslice}_INTRP ${outfile}


