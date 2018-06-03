#!/bin/bash

# Tested in Ubuntu 16.04

source "$(dirname $(readlink -f ${BASH_SOURCE[0]}))/bash-utils.sh"

require ncks wrf_interp.exe ncatted ncrename cdo nccopy

nc-select-steps () {

    local first="$1"
    local last="$2"
    local infile="$3"
    local result="$4"

    local tmpfile=$(mktemp --suffix=.nc)
    trap "{ rm -f $tmpfile; }" EXIT

    ${NCKS} -O -h -d Time,${first},${last} $infile $tmpfile

    ${NCATTED} -O -h -a history,global,d,, $tmpfile

    local startdate=$(ncdump -v Times $tmpfile \
                             | sed -nr '/^ *Times = *$/{n;s/^.*([0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}:[0-9]{2}:[0-9]{2}).*/\1/p;q}')

    local outfile=${infile%%_d01_*}_d01_${startdate}

    mv $tmpfile $outfile

    printf -v "$result" "%s" "$outfile"
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
    local formulas="
        U = UU;
        V = VV;
        T = TT;
        Z = GHT;
        PS = PSFC;
        Q = ((RTHRATEN + RTHBLTEN) / (MU + MUB) + H_DIABATIC) * (clev(RTHRATEN) / 100000) ^ (287 / 1004);
        FU = RUBLTEN / (MU + MUB);
        FV = RVBLTEN / (MU + MUB);
        F = F;
        lev = LEV;"
    local dxdydt=$(ncdump -h $infile | sed -nr 's/:D. = ([0-9.]+).*/\1/p')

    local tmpfile=$(mktemp --suffix=.nc)
    trap "{ rm -f $tmpfile; }" EXIT

    IGNORE_ATT_COORDINATES=1 ${CDO} expr,"$formulas" ${infile} ${tmpfile}

    ${NCATTED} -O -a units,,d,, -a long_name,,d,, $tmpfile

    attr='lev  "Pa"       "Pressure levels"
          U    "m s**-1"  "U velocity"
          V    "m s**-1"  "V velocity"
          T    "K"        "Temperature"
          Z    "m"        "Geopotential height"
          PS   "hPa"      "Surface pressure"
          Q    "K s**-1"  "Diabatic heating"
          FU   "m s**-2"  "Frictional U-tendency"
          FV   "m s**-2"  "Frictional U-tendency"'

    ${NCATTED} -O $(awk '{printf "-a units,%s,c,c,%s -a long_name,%s,c,c,%s\n", $1, $2, $1, $3}' <<<"$attr") $tmpfile

    #${NCATTED} -O -h -a history,global,d,, -a history_of_appended_files,global,d,, $tmpfile
    ${NCATTED} -O -h -a ,global,d,, $tmpfile

    ${NCRENAME} -h -d vlevs,lev -d west_east,x -d south_north,y ${tmpfile}

    set -- $dxdydt
    ncatted -O -h -a DX,global,o,d,$1 -a DY,global,o,d,$2 -a DT,global,o,d,$3 ${tmpfile}

    ${NCCOPY} -k nc4 $tmpfile $outfile
}
