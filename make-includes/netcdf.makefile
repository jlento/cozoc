# Netcdf4 library with parallel IO support using HDF5

NETCDF_PREFIX ?= $(realpath $(CURDIR))

NETCDF_FTP    = ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.tar.gz
NETCDF_TARGETS       = $(NETCDF_PREFIX)/include/netcdf_par.h $(NETCDF_PREFIX)/lib/libnetcdf.a

export CPPFLAGS != pkg-config --cflags hdf5
export LDFLAGS  != pkg-config --libs hdf5

.PHONY : all

all : $(NETCDF_TARGETS)

$(NETCDF_TARGETS) : $(notdir $(NETCDF_FTP))
	tar xf $<
	cd $(<:.tar.gz=); ./configure --disable-shared --enable-parallel-tests --prefix=$(NETCDF_PREFIX); $(MAKE) install

$(notdir $(NETCDF_FTP)) :
	wget $(NETCDF_FTP)

$(NETCDF_PREFIX)/lib/libnetcdf.a : $(NETCDF_PREFIX)/include/netcdf_par.h
