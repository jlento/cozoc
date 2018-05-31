#!/bin/bash

# Tested in Ubuntu 16.04

source "$(dirname $(readlink -f ${BASH_SOURCE[0]}))/bash-utils.sh"

require ncks wrf_interp.exe ncatted ncrename

nc-select-steps () {
    exec 3>&1
    exec 1> /dev/null
    local first="$1"
    local last="$2"
    local infile="$3"
    local tmpfile=$(mktemp --suffix=.nc)
    #trap "{ rm -f $tmpfile; }" EXIT
    ${NCKS} -O -h -d Time,${first},${last} $infile $tmpfile
    ${NCATTED} -O -h -a history,global,d,, $tmpfile
    local startdate=$(ncdump -v Times $tmpfile \
                             | sed -nr '/^ *Times = *$/{n;s/^.*([0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}:[0-9]{2}:[0-9]{2}).*/\1/p;q}')
    local outfile=${infile%%_d01_*}_d01_${startdate}
    mv $tmpfile $outfile
    exec 1>&3
    echo $outfile
}

wrf-interpolate () {
    local dp="$1"
    local infile="$2"
    local rootname id startdate time
    IFS=_ read rootname id startdate time <<<"$infile"
    rm -f "${infile}_INTRP"
    ${WRF_INTERP_EXE} <<<"
        &io
          path_to_input = '$(dirname $infile)'
          path_to_output = '$(dirname $infile)'
          root_name = '${rootname}'
          grid_id = 1
          start_date = '${startdate}_${time%%:*}'
          leap_year  = .FALSE.
          debug = .TRUE.
        /
        &interp_in
          interp_levels = 1000,-100,${dp}
          extrapolate = 1
          unstagger_grid = .TRUE.
          vert_coordinate = 'pres'
        /"
}

wrf-fix-dimensions-etc () {
    local infile="$1"
    local outfile="$2"
    ${NCKS} -O -h -x -C -v time ${infile} ${outfile}
    ${NCRENAME} -O -h \
                -v XTIME,time \
                -d vlevs,lev \
                -d west_east,lon \
                -d south_north,lat ${outfile}
    ${NCATTED} -O -h \
            -a coordinates,,d,, \
            -a history,global,d,, ${outfile}
}

return

formulas="
U=UU;
V=VV;
T=TT;
Z=GHT;
PS=PSFC;
Q=((RTHRATEN+RTHBLTEN)/(MU+MUB)+H_DIABATIC) * (clev(RTHRATEN)/100000)^(287/1004);
FU=RUBLTEN/(MU+MUB);
FV=RVBLTEN/(MU+MUB);
F=F;
lev=LEV;
"

# Start processing

# Get grid spacings
dxdydt=$(ncdump -h wrfout_d01_0001-01-05_01\:00\:00_INTRP | sed -nr 's/:D. = ([0-9.]+).*/\1/p')


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
PS   '"hPa"'      '"Surface pressure"'
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

# Add attributes
set -- $dxdydt
ncatted -O -a DX,global,o,d,$1 -a DY,global,o,d,$2 -a DT,global,o,d,$3 ${fname}

# Rename dimensions


# Convert the output file to netcdf-4 format
$NCCOPY -k nc4 $fname ${fname}4
rm -rf ${infile}_1
