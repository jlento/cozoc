ROOTDIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
SHELL   := /bin/bash
PROG    := cozoc

define USAGE

Usage:

  make -f <makefile> [help|$(PROG)|test|tags] \$(eval)
      [PETSC_DIR=...] [NETCDF_DIR=...] [MPIEXEC=...]

endef


.PHONY : all help clean test tags edit netcdf


# Regular trgets

all :
	$(MAKE) -f $(ROOTDIR)/src/makefile PROG=$(PROG)

help :
	@:$(info $(USAGE))

clean :
	rm -rf $(PROG) *.o *.d.*


# Test case

INPUT     := wrf.nc4
INPUT_URL := https://bitbucket.org/mikarant/ozo/downloads/test_wrf.nc
NCCOPY    := $(shell compgen -G $(CURDIR)/bin/nccopy || compgen -c nccopy)

test : all $(INPUT)
	$(MPIEXEC) ./$(PROG) -Q -G

$(INPUT) : $(notdir $(INPUT_URL))
ifeq ($(NCCOPY),)
	$(error Command nccopy not found. Install netcdf and/or set variable NCCOPY.)
else
	$(NCCOPY) -k nc4 $< $@
endif

$(notdir $(INPUT_URL)) :
	wget -O $@ $(INPUT_URL)


# GTAGS

tags :
	ln -sfT $(PETSC_DIR)/include $(ROOTDIR)/src/petsc-includes
	ln -sfT $(NETCDF_DIR)/include $(ROOTDIR)/src/netcdf-includes
	cd $(ROOTDIR)/src; gtags


# Emacs IDE

edit :
	cd $(ROOTDIR)/src; emacs -q -l $(ROOTDIR)/ide/init.el &


# Build Netcdf4/HDF5

netcdf :
	$(MAKE) -f $(ROOTDIR)/make-includes/netcdf.makefile
