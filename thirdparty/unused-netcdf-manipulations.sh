#!/bin/bash

# Modify global attributes without changing history attribute
#
#ncatted -h -a START_DATE,global,m,c,${startdate}:00:00 \
#        -a SIMULATION_START_DATE,global,m,c,${startdate}:00:00 \
#        wrfout_d01_${startdate}:00:00

# Remove history and NCO global attributes
#
#ncatted -hO -a history,global,d,, -a NCO,global,d,, wrfout_d01_${startdate}:00:00

