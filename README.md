COZOC
=====

Cozoc is parallel re-implementation of [OZO](https://github.com/mikarant/ozo).

Status
------

- Computes the quasi-geostrophic omega equation and the first term
  (vorticity advection) of the generalized omega equation, only
- very much under construction...


Developing or testing COZOC in a virtual machine
------------------------------------------------

The easiest way to get started with COZOC is to clone the whole
development environment. Install Git, VirtualBox and Vagrant on your
laptop (the host machine). Then:

    git clone https://jlento@bitbucket.org/jlento/cozoc.git
    vagrant up

The first boot of the virtual machine takes a while. After the dust
has settled, you should have Ubuntu 16.04 LTS desktop running in a
virtual machine, with all packages necessary for COZOC already
installed. The username and password for the virtual machine are both
"vagrant", with sudo priviledges, as usual. You may need to fix the
keyboard mapping "System Settings.../Text Entry" and time zone "System
Settings.../Time & Date" in the virtual machine.

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

The parallel version of Netcdf4/HDF5 is already provided, so:

    module swap PrgEnv-gnu/5.2.82 PrgEnv-gnu
    module load cray-petsc cray-hdf5-parallel cray-netcdf-hdf5parallel
    make -f ../cozoc/makefile test MPIEXEC=aprun


Build with cmake
----------------------------

...on the way...


TODO:
-----

- fix gtags support and emacs rc-file init.el
- omegaQG.c was written before omega.c. Now omegaQC should be updated to
  same style and using the same subroutines as omega.
- drop vecops.c?
- put wrf netcdf input file field names into wrfnc.h and wrfnc.c files,
  similar to dimension names

