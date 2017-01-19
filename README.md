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

Building COZOC requires C compiler, CMake, and the following
libraries:

- MPI
- PETSc
- HDF5 (parallel)
- NetCDF4/HDF5 (parallel)


Testing and developing COZOC in Ubuntu 16.04 LTS
---------------------------------------------------------

The easiest way to get started with COZOC is to clone the whole
development environment. Install Git, VirtualBox and Vagrant on your
laptop (the host machine). Then:

    vagrant up

The first boot of the virtual machine takes a while. After the dust
has settled, you should have Ubuntu 16.04 LTS virtual machine running
in the background, with most of the packages necessary for testing and
developing COZOC already installed. See files [Vagrantfile](Vagrantfile)
and[playbook.yml](playbook.yml).

Login to the virtual machine with

    vagrant ssh -- -Y

At the time of writing, Ubuntu/Debian does have a package for the parallel
version of HDF5, but not for NetCDF4, which needs to be compiled separately
from the sources. If parallel version of NetCDF4/HDF5 library needs to be
build from the sources, build and install it into the cozoc build directory
using the provided [netcdf4.bash](netcdf4.bash) script:

    cd /vagrant/build
    bash ../netcdf4.bash


Configure COZOC with CMake
-------------------------------

COZOC build is configured using CMake. In the provided virtual machine,
with the NetCDF4 library build from sources as described above,

    cd build
    cmake ..
    
should work.

If cmake cannot auto-detect the location of the PETSc or NetCDF4 libraries,
you can set environment variables PETSC_DIR or NETCDF_DIR, to give cmake
a hint where to search for the libraries. For example,

    NETCDF_DIR=$HOME/my-netcdf cmake ..

In supercomputer environments the libraries are often made available through
environment module system. For example, in Cray XC40, the commands

    module swap PrgEnv-cray PrgEnv-gnu
    module load cmake cray-petsc
    module load cray-hdf5-parallel cray-netcdf-hdf5parallel

set up the environment for building COZOC. As the compile wrapper `cc` now
takes care of the compile and link flags, we need to 




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

