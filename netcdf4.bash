#!/bin/bash

# Build and install parallel I/O enabled NetCDF4
#
# Requires:
# - HDF5 build with '--enable-parallel'
#
# Usage:
#   bash <this file> [URL=<URL>] [DESTDIR=<DESTDIR>]
# where
#   <URL>     -- the download URL of the netcdf source
#   <DESTDIR> -- the install directory

: ${URL:=ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.tar.gz}
: ${DESTDIR:=$PWD/netcdf}

PKG=${URL##*/}
DIR=${PKG%.tar.gz}

test -f $PKG || wget $URL
test -d $DIR || tar xvf $PKG
cd $DIR
./configure --prefix=${DESTDIR} \
            CFLAGS="$(pkg-config --cflags hdf5)" \
            LDFLAGS="$(pkg-config --libs hdf5)"
make install
