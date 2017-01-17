COZOC
=====

Cozoc is parallel re-implementation of [OZO](https://github.com/mikarant/ozo).


Status
------

- Computes the quasi-geostrophic omega equation and the first term
  (vorticity advection) of the generalized omega equation, only
- very much under construction...


Obtaining COZOC source
----------------------

Clone the COZOC GitHub repository:

    git clone https://jlento@bitbucket.org/jlento/cozoc.git


COZOC requirements
------------------

Building COZOC requires GNU C compiler (gcc), CMake, and the following
libraries:

- MPI
- PETSc
- (Parallel) Netcdf/HDF5

On regular Linux workstations these dependencies can be installed
using the system's package manager. The names of the required packages
for Ubuntu 16.04 LTS, for example, can be found from file
[playbook.yml](playbook.yml), along with some of the development tools
that I use.

In practice, these libraries are readily available in supercomputer
environments through environment module system. For example, in Cray
XC40, the systems commands

    module swap PrgEnv-cray PrgEnv-gnu
    module load cmake cray-petsc
    module load cray-hdf5-parallel cray-netcdf-hdf5parallel

are enough to prepare the environment for building COZOC.


Testing and developing COZOC in a virtual machine
-------------------------------------------------

The easiest way to get started with COZOC is to clone the whole
development environment. Install Git, VirtualBox and Vagrant on your
laptop (the host machine). Then:

    vagrant up

The first boot of the virtual machine takes a while. After the dust
has settled, you should have Ubuntu 16.04 LTS virtual machine running
in the background, with all packages necessary for testing and
developing COZOC already installed.


Configure COZOC with CMake
-------------------------

COZOC is configured as usual with CMake:

    cd build
    cmake ..

The configuration should work out of the box in the provided virtual machine.
Tips for configuring in other environments are provided below.


Building COZOC
---------------

After CMake configuration, just type

    make

to build COZOC. In case of troubles, first try `make VERBOSE=1`.


Spacemacs C IDE
----------------

The project includes `.spacemacs` configuration file that can be used
with [Spacemacs](http://spacemacs.org).

The above `cmake`-command also writes files `.clang_complete` and
`.dir-locals.el` files to the project root, that are used by the Spacemacs
packages. Thus, to get the full features of the Spacemacs C IDE working, run the
CMake configuration first.


Building parallel NetCDF4/HDF5
-------------------------------

*NOTE: This chapter is outdated and under construction!!!*

Basic ubuntu repositories do not(?) have the parallel Netcdf4/HDF5, so
we need to build that library first from the sources:

    mkdir -p build
    cd $_
    make -f ../cozoc/makefile netcdf

With netcdf in place

    make -f ../cozoc/makefile test

should build cozoc executable, load the test input file from a WRF
simulation, and run cozoc generating the omega and height-tendency
fields.



Building and testing COZOC in Cray XC40
---------------------------------------

*NOTE: This chapter is under construction!!!*

*NOTE: Should definitely write Toolchain files for cmake in Cray*

INTEL:

    cmake .. -DCMAKE_C_FLAGS_RELEASE="-std=gnu99" \
        -DPETSC_INCLUDE_DIR= -DNETCDF_INCLUDE_DIR= \
	    -DNETCDF_LIBRARY= -DPETSC_LIBRARY= -DUSE_PARALLEL_NETCDF=1 \
		-DCMAKE_SYSTEM_NAME=Cray -DCMAKE_C_COMPILER=cc

GNU & CRAY:

    cmake .. -DPETSC_INCLUDE_DIR= -DNETCDF_INCLUDE_DIR= \
        -DNETCDF_LIBRARY= -DPETSC_LIBRARY= -DUSE_PARALLEL_NETCDF=1


TODO:
-----

- omegaQG.c was written before omega.c. Now omegaQC should be updated to
  same style and using the same subroutines as omega.
- drop vecops.c?
- put wrf netcdf input file field names into wrfnc.h and wrfnc.c files,
  similar to dimension names

