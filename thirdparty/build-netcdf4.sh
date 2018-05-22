#!/bin/bash

# Build and install parallel I/O enabled NetCDF4
#
# Requires:
# - HDF5 build with '--enable-parallel'
#
# Usage:
#   [URL=<URL>] [INSTALLDIR=<INSTALLDIR>] bash <this file>
# where
#   <URL>        -- the download URL of the netcdf source
#   <INSTALLDIR> -- the install directory

: ${URL:=ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.tar.gz}
: ${INSTALLDIR:=$PWD/netcdf}

PKG=${URL##*/}
DIR=${PKG%.tar.gz}

test -f $PKG || wget $URL
test -d $DIR || tar xvf $PKG
cd $DIR
./configure --prefix=${INSTALLDIR} \
            CFLAGS="$(pkg-config --cflags hdf5-openmpi)" \
            LDFLAGS="$(pkg-config --libs hdf5-openmpi)"
make install
