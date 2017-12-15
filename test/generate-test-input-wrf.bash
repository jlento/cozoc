#!/bin/bash

nx=8
ny=8
nz=4
nt=4
dt=60

function csv_repeat () {
    for ((i=1;i<$1;i++))
    do
        printf "$2,"
    done
    echo "$2"
}

cat <<EOF
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
		:history = "Fri Sep  9 10:43:37 2016: ncrename -O -d x,DateStrLen test_wrf.nc\nFri Sep  9 10:27:52 2016: ncrename -O -d x_2,west_east test_wrf.nc\nFri Sep  9 10:27:26 2016: ncrename -O -d y,south_north test_wrf.nc\nFri Sep  9 10:26:53 2016: ncrename -O -d lev,vlevs test_wrf.nc\nFri Sep 09 10:22:24 2016: cdo seltimestep,117,118,119,120 wrfout_d01_0001-01-01_INTRP.nc test_wrf.nc" ;
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
	UU = $(csv_repeat $((nt*nz*ny*nx)) 1);
	VV = $(csv_repeat $((nt*nz*ny*nx)) 0);
	MU = $(csv_repeat $((nt*ny*nx)) 0);
	MUB = $(csv_repeat $((nt*ny*nx)) 0);
	PSFC = $(csv_repeat $((nt*ny*nx)) 10000);
	H_DIABATIC = $(csv_repeat $((nt*nz*ny*nx)) 0);
	F = $(csv_repeat $((nt*ny*nx)) 0);
	RUCUTEN = $(csv_repeat $((nt*nz*ny*nx)) 0);
	RVCUTEN = $(csv_repeat $((nt*nz*ny*nx)) 0);
	RTHCUTEN = $(csv_repeat $((nt*nz*ny*nx)) 0);
	RTHRATEN = $(csv_repeat $((nt*nz*ny*nx)) 0);
	RUBLTEN = $(csv_repeat $((nt*nz*ny*nx)) 0);
	RVBLTEN = $(csv_repeat $((nt*nz*ny*nx)) 0);
	RTHBLTEN = $(csv_repeat $((nt*nz*ny*nx)) 0);
	TT = $(csv_repeat $((nt*nz*ny*nx)) 293);
	GHT = $(csv_repeat $((nt*nz*ny*nx)) 1);
	XTIME = $(seq -s , 0 $dt $((dt*(nt-1))));
}
EOF
