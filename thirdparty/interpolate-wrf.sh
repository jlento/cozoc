#!/bin/bash

set -e

# Tested in Ubuntu 16.04

: ${WRF_ROOT:=${HOME}/thirdparty}

mkdir -p $WRF_ROOT/WRFV3/test/em_b_wave/interpolated
cd $WRF_ROOT/WRFV3/test/em_b_wave/interpolated

simulation_start_date='0001-01-01_00'
dt=1 # hour
day=5
hour=1
nsteps=4
(( first = 24 * (day - 1) + hour ))
(( last = first + nsteps - 1 ))
interpolation_start_date=0001-01-$(printf "%02u" $day)_$(printf "%02u" $hour)

echo "
&io
 path_to_input = '$PWD'
 path_to_output = '$PWD'
 root_name = 'wrfout'
 grid_id = 1
 start_date = '$interpolation_start_date'
 leap_year  = .FALSE.
 debug = .TRUE.
/
&interp_in
  interp_levels = 1000,-100,20
  extrapolate = 0
  unstagger_grid = .TRUE.
  vert_coordinate = 'pres'
/" > namelist.vinterp

ncks -O -a -h -d Time,${first},${last} ../wrfout_d01_${simulation_start_date}:00:00 wrfout_d01_${simulation_start_date}:00:00

rm -f wrfout_d01_${interpolation_start_date}:00:00_INTRP*
$WRF_ROOT/wrf_interp/wrf_interp.exe
