#!/bin/bash

curdir=$(readlink -f $(dirname ${BASH_SOURCE[0]}))
output_filename=${1:-$(basename -s .bash "$0").nc4}
: ${ncgen:=${curdir}/../netcdf/bin/ncgen}

nx=8
ny=8
nz=4
nt=8
dt=60

# Loosely from the steady state vortex solution
#   https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations

field () {
    bc -ql <<EOF
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
        $1
      }
EOF
}

A=1
B=1
r=1

rho=( $(field "3 * $B / ($r^2 + x^2 + y^2 + z^2)") )

p=(   $(field "-$A^2 * $B / ($r^2 + x^2 + y^2 + z^2)^3") )

u=(   $(field "2 * $A / ($r^2 + x^2 + y^2 + z^2)^2 * (-$r * y + x * z)") )

v=(   $(field "2 * $A / ($r^2 + x^2 + y^2 + z^2)^2 * ($r * x + y * z)") )

w=(   $(field "2 * $A / ($r^2 + x^2 + y^2 +z^2)^2 * ($r^2 - x^2 - y^2 - x^2)") )

T=(   $(field "-$A^2 * $B / ($r^2 + x^2 + y^2 +z^2)^3 - y + 273") )

IFS=,

function csv_repeat () {
    for ((i=1;i<$1;i++)); do printf "$2,"; done; printf "$2"
}

cat <<EOF | ${ncgen} -k nc4 -o "$output_filename"
netcdf wrf {
dimensions:
	west_east = ${nx} ;
	south_north = ${ny} ;
	vlevs = ${nz} ;
	time = UNLIMITED ; // (${nt} currently)
variables:
	float LEV(vlevs) ;
		LEV:units = "PaH" ;
		LEV:description = "Pressure Levels" ;
	float UU(time, vlevs, south_north, west_east) ;
		UU:units = "m s-1" ;
		UU:description = "x-wind component" ;
	float VV(time, vlevs, south_north, west_east) ;
		VV:units = "m s-1" ;
		VV:description = "y-wind component" ;
	float MU(time, south_north, west_east) ;
		MU:units = "Pa" ;
		MU:description = "perturbation dry air mass in column" ;
	float MUB(time, south_north, west_east) ;
		MUB:units = "Pa" ;
		MUB:description = "base state dry air mass in column" ;
	float PSFC(time, south_north, west_east) ;
		PSFC:units = "Pa" ;
		PSFC:description = "SFC PRESSURE" ;
	float H_DIABATIC(time, vlevs, south_north, west_east) ;
		H_DIABATIC:units = "K s-1" ;
		H_DIABATIC:description = "MICROPHYSICS LATENT HEATING" ;
	float F(time, south_north, west_east) ;
		F:units = "s-1" ;
		F:description = "Coriolis sine latitude term" ;
	float RUCUTEN(time, vlevs, south_north, west_east) ;
		RUCUTEN:units = "Pa m s-2" ;
		RUCUTEN:description = "COUPLED X WIND TENDENCY DUE TO CUMULUS PARAMETERIZATION" ;
	float RVCUTEN(time, vlevs, south_north, west_east) ;
		RVCUTEN:units = "Pa m s-2" ;
		RVCUTEN:description = "COUPLED Y WIND TENDENCY DUE TO CUMULUS PARAMETERIZATION" ;
		RVCUTEN:stagger = "" ;
	float RTHCUTEN(time, vlevs, south_north, west_east) ;
		RTHCUTEN:units = "Pa K s-1" ;
		RTHCUTEN:description = "COUPLED THETA TENDENCY DUE TO CUMULUS SCHEME" ;
	float RTHRATEN(time, vlevs, south_north, west_east) ;
		RTHRATEN:units = "Pa K s-1" ;
		RTHRATEN:description = "COUPLED THETA TENDENCY DUE TO RADIATION" ;
	float RUBLTEN(time, vlevs, south_north, west_east) ;
		RUBLTEN:units = "Pa m s-2" ;
		RUBLTEN:description = "COUPLED X WIND TENDENCY DUE TO PBL PARAMETERIZATION" ;
	float RVBLTEN(time, vlevs, south_north, west_east) ;
		RVBLTEN:units = "Pa m s-2" ;
		RVBLTEN:description = "COUPLED Y WIND TENDENCY DUE TO PBL PARAMETERIZATION" ;
	float RTHBLTEN(time, vlevs, south_north, west_east) ;
		RTHBLTEN:units = "Pa K s-1" ;
		RTHBLTEN:description = "COUPLED THETA TENDENCY DUE TO PBL PARAMETERIZATION" ;
	float TT(time, vlevs, south_north, west_east) ;
		TT:units = "K" ;
		TT:description = "Temperature" ;
	float GHT(time, vlevs, south_north, west_east) ;
		GHT:units = "M" ;
		GHT:description = "Geopotential Height" ;
	double XTIME(time) ;

// global attributes:
		:history = "Some almost arbitrary numbers for testing purposes" ;
		:TITLE = " SIMULATED OUTPUT FROM WRF V3.8.1 MODEL - ON PRES LEVELS" ;
		:START_DATE = "0001-01-01_00:00:00" ;
		:SIMULATION_START_DATE = "0001-01-01_00:00:00" ;
		:WEST-EAST_GRID_DIMENSION = 41 ;
		:SOUTH-NORTH_GRID_DIMENSION = 81 ;
		:BOTTOM-TOP_GRID_DIMENSION = 19 ;
		:DX = 100000.f ;
		:DY = 100000.f ;
		:USE_Q_DIABATIC = 0 ;
		:GRIDTYPE = "C" ;
		:CU_PHYSICS = 1 ;
		:DT = 600.f ;
		:MAP_PROJ_CHAR = "Cartesian" ;
data:
	LEV = $(seq -s , 100000 -$((90000/(nz-1))) 10000);
	UU = ${u[*]};
	VV = ${v[*]};
	MU = $(csv_repeat $((nt*ny*nx)) 1);
	MUB = $(csv_repeat $((nt*ny*nx)) 0);
	PSFC = $(csv_repeat $((nt*ny*nx)) 10000);
	H_DIABATIC = ${w[*]};
	F = $(csv_repeat $((nt*ny*nx)) 1);
	RUCUTEN = ${w[*]};
	RVCUTEN = $(csv_repeat $((nt*nz*ny*nx)) 1);
	RTHCUTEN = $(csv_repeat $((nt*nz*ny*nx)) 1);
	RTHRATEN = $(csv_repeat $((nt*nz*ny*nx)) 1);
	RUBLTEN = $(csv_repeat $((nt*nz*ny*nx)) 1);
	RVBLTEN = $(csv_repeat $((nt*nz*ny*nx)) 1);
	RTHBLTEN = $(csv_repeat $((nt*nz*ny*nx)) 1);
	TT = ${T[*]};
	GHT = ${p[*]};
	XTIME = $(seq -s , 0 $dt $((dt*(nt-1))));
}
EOF
