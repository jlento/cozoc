#!/bin/bash

set -ex

# Tested in Ubuntu 16.04

: ${WRF_ROOT:=${HOME}/thirdparty}
: ${CDO:=${HOME}/thirdparty/bin/cdo}
: ${NCCOPY:=${HOME}/thirdparty/bin/nccopy}

# OS X workaround
readlink-f () {
    local absname
    readlink "$1" || ( test -d "$1" && cd "$1" && pwd -P ) || (
        cd $(dirname "$1") 2> /dev/null && echo "$(pwd -P)/$(basename "$1")"
    )
}

srcdir="$(dirname $(readlink-f ${BASH_SOURCE[0]}))"

cd $WRF_ROOT/WRFV3/test/em_b_wave/interpolated

infile=$(ls -1 *_INTRP | head -1)
fname=wrf.nc

formulas="U=UU;
V=VV;
T=TT;
Z=GHT;
SP=PSFC;
Q=((RTHRATEN+RTHBLTEN)/(MU+MUB)+H_DIABATIC) * (clev(RTHRATEN)/100000)^(287/1004);
FU=RUBLTEN/(MU+MUB);
FV=RVBLTEN/(MU+MUB);
F=F;
LEV=LEV;"

# Start processing

# Fix time coordinate
ncks -h -C -O -x -v time ${infile} ${infile}_1
ncatted -O -a coordinates,,d,, ${infile}_1
ncrename -O -v XTIME,time ${infile}_1

# Compute fields
IGNORE_ATT_COORDINATES=1 ${CDO} expr,"$formulas" ${infile}_1 ${fname}

# Remove unnecessary attributes from the variables
ncatted -O -a units,,d,, -a long_name,,d,, $fname

# Add units and long names to the variables
arr=(
LEV  '"Pa"'       '"Pressure levels"'
U    '"m s**-1"'  '"U velocity"'
V    '"m s**-1"'  '"V velocity"'
T    '"K"'        '"Temperature"'
Z    '"m"'        '"Geopotential height"'
SP   '"hPa"'      '"Surface pressure"'
Q    '"K s**-1"'  '"Diabatic heating"'
FU   '"m s**-2"'  '"Frictional U-tendency"'
FV   '"m s**-2"'  '"Frictional U-tendency"'
)

cmd="ncatted -O"
i=0
while (( i < ${#arr[@]} )); do
    cmd="$cmd -a units,${arr[$i]},c,c,${arr[$((i+1))]} -a long_name,${arr[$i]},c,c,${arr[$((i+2))]}"
    (( i += 3 ))
done

eval $cmd $fname

# Remove the long history-attribute
ncatted -O -a history,global,d,, $fname
ncatted -O -a history_of_appended_files,global,d,, $fname

# Rename dimensions
ncrename -d vlevs,lev -d west_east,x -d south_north,y ${fname}

# Convert the output file to netcdf-4 format
$NCCOPY -k nc4 $fname ${fname}4
rm -rf ${infile}_1
