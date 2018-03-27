#!/bin/bash

curdir=$(readlink -f $(dirname ${BASH_SOURCE[0]}))

ncgen=$(which ncgen 2> /dev/null)
: ${ncgen:=${curdir}/../netcdf/bin/ncgen}

# Defaults
size=16
template='wrf'

usage="
Usage: $0 [-s SIZE] [-o OUTPUT] [-t TEMPLATE]

Optional arguments:
  SIZE     -- Each field will have SIZE^4/2 points (Default: $size)
  TEMPLATE -- The basename for the ncgen (CDL) input file
              'TEMPLATE.template' (Default: 'wrf')
  OUTPUT   -- Name of the output file (Default: 'TEMPLATE.nc4')
"

# Command line options
OPTS=$(getopt -o hs:o:t: --long help,output:,size:,template: -n 'parse-options' -- "$@")
if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; echo "$usage" >&2; exit 1 ; fi
eval set -- "$OPTS"
while true; do
    case "$1" in
        -h | --help ) echo "$usage"; exit 0;;
        -o | --output ) output="$2"; shift; shift ;;
        -s | --size ) size="$2"; shift; shift ;;
        -t | --template ) template="$2"; shift; shift ;;
        -- ) shift; break ;;
        * ) break ;;
  esac
done

: ${output:=${template}.nc4}

nx=$size
ny=$size
((nz = size / 2))
nt=$size
dt=60

# Functions and variables (fields) that can be used in the template
#
# Loosely from the steady state vortex solution
#   https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations

field () {
    bc -ql <<EOF | tr -d '[\\\n]' | head -c -1
scale = 2
nt=$nt
nz=$nz
ny=$ny
nx=$nx

define pow (a,b) { return (e(b * l(a))); }
define div (a,b) { auto s,r; s = scale; scale = 0; r = a / b; scale = s; return r; }
define rem (a,b) { auto s,r; s = scale; scale = 0; r = a % b; scale = s; return r; }
for (i = 0;i < nx; i++) periodic[i] = i - div(nx,2)

for (t = 0; t < nt; t++)
  for (z = - div(nz,2); z < nz - div(nz,2); z++)
    for (y = - div(ny,2); y < ny - div(ny,2); y++)
      for (i = 0; i < nx; i++) {
        x = periodic[rem(i+div(nt,2)-t+nt*nx,nx)]
        print $1, ","
      }
EOF
}

A=1
B=1
(( r = size / 8 ))

rho=$(field "3 * $B / ($r^2 + x^2 + y^2 + z^2)")
p=$(field "-$A^2 * $B / ($r^2 + x^2 + y^2 + z^2)^3")
u=$(field "2 * $A / ($r^2 + x^2 + y^2 + z^2)^2 * (-$r * y + x * z)")
v=$(field "2 * $A / ($r^2 + x^2 + y^2 + z^2)^2 * ($r * x + y * z)")
w=$(field "2 * $A / ($r^2 + x^2 + y^2 +z^2)^2 * ($r^2 - x^2 - y^2 - x^2)")
T=$(field "-$A^2 * $B / ($r^2 + x^2 + y^2 +z^2)^3 - y + 273")

function csv_repeat () {
    for ((i=1;i<$1;i++)); do printf "$2,"; done; printf "$2"
}

eval 'echo "'"$(sed 's/\\/\\\\/g;s/\"/\\\"/g' ${curdir}/${template}.template)"'"' | ${ncgen} -k nc4 -o "$output"
