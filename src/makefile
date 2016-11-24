SRCDIR     := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
SHELL      := /bin/bash
PROG       ?= cozoc
PKG_CONFIG  = PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --silence-errors

include $(SRCDIR)/../make-includes/constants.makefile
include $(SRCDIR)/../make-includes/findpetsc.makefile
include $(SRCDIR)/../make-includes/findnetcdf.makefile

NETCDF_INCS != $(PKG_CONFIG) --cflags netcdf
NETCDF_LIBS != $(PKG_CONFIG) --libs netcdf

PETSC_INCS  != $(PKG_CONFIG) --cflags PETSc
PETSC_LIBS  != $(PKG_CONFIG) --libs PETSc

###########################################################################

CC     != $(PKG_CONFIG) --variable=ccompiler PETSc
COPTS   = -g -pg -Wall -fbounds-check
CFLAGS  = $(COPTS) -I$(SRCDIR) $(NETCDF_INCS) $(PETSC_INCS)
LDFLAGS = $(NETCDF_LIBS) $(PETSC_LIBS)

SRC     = cozoc.c context.c idtest.c solve.c omegaQG.c omega.c io.c vecops.c \
          daslice.c derivatives.c utils.c constants.c coriolis_parameter.c \
          sigma_parameter.c vorticity.c wrfnc.c
OBJ     = $(SRC:.c=.o)

vpath %.c $(SRCDIR)
vpath %.h $(SRCDIR)
vpath %.py $(SRCDIR)
vpath %.inc $(SRCDIR)

.PHONY : all

all : $(PROG)

$(PROG) : $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

%.d : %.c
	set -e; rm -f $@; \
         $(CC) -MM $(CFLAGS) $< > $@.$$$$; \
         sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
         rm -f $@.$$$$

include $(SRC:.c=.d)

omega.d : omega_stencil.inc

$(SRCDIR)/omega_stencil.inc : omega_stencil.py
	python $< > $@



# Configure (generate config.h) with make

config.h :
	$(file >config.h,$(config-current))
	$(info -- Wrote config.h)

define config-current
#ifndef CONFIG_H
#define CONFIG_H

#$(if $(USE_PARALLEL_NETCDF),define,undef) USE_PARALLEL_NETCDF

#endif /* CONFIG_H */
endef

ifneq ($(strip $(config-current)),$(strip $(shell cat config.h 2>&-)))
    $(shell rm -f config.h)
endif


