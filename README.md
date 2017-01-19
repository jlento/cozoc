# COZOC

Cozoc is parallel re-implementation of [OZO](https://github.com/mikarant/ozo).


## Status

- Computes the quasi-geostrophic omega equation and the first term
  (vorticity advection) of the generalized omega equation, only
- very much under construction...


## Obtaining COZOC source

Clone the COZOC GitHub repository:

    git clone https://jlento@bitbucket.org/jlento/cozoc.git


## COZOC requirements

Building COZOC requires C compiler, CMake, and the following
libraries:

- MPI
- PETSc
- HDF5 (parallel)
- NetCDF4/HDF5 (parallel)


## Testing and developing COZOC in Ubuntu 16.04 LTS

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

The `cozoc` source directory in the host machine is shared with the
virtual machine. In the virtual machine it is visible as a directory
`/vagrant`. 

It is a good practice to have separate source and build directories
(out-of-source builds). Let's create a build directory and continue in it,
with for example,

    mkdir $HOME/build
    cd $_

At the time of writing, Ubuntu/Debian does have a package for the parallel
version of HDF5, but not for NetCDF4, which needs to be compiled separately
from the sources. If parallel version of NetCDF4/HDF5 library needs to be
build from the sources, you can use the provided [netcdf4.bash](netcdf4.bash)
script:

    bash /vagrant/netcdf4.bash


## Configure COZOC build with CMake

### In Ubuntu 16.04 LTS

COZOC build is configured using CMake. In the provided virtual machine,
with the NetCDF4 library build from sources as described above,

    cmake /vagrant
    
command in the build directory should configure COZOC automatically.

If cmake cannot auto-detect the location of the PETSc or NetCDF4 libraries,
you can set environment variables PETSC_DIR or NETCDF_DIR, to give cmake
a hint where to search for the libraries. For example,

    NETCDF_DIR=$HOME/my-netcdf cmake ..

### In Cray XC40

In supercomputer environments the libraries are often made available through
environment module system. For example, in Cray XC40, the commands

    module load cmake cray-petsc
    module load cray-hdf5-parallel cray-netcdf-hdf5parallel

set up the environment for building COZOC. As the compile wrapper `cc` now
takes care of the compile and link flags, we need to instruct cmake not to
try to add all the flags a second time. Depending on the loaded compile
environment, the cmake command could be

INTEL:

    cmake .. -DCMAKE_C_FLAGS_RELEASE="-std=gnu99" \
        -DPETSC_INCLUDE_DIR= -DNETCDF_INCLUDE_DIR= \
	    -DNETCDF_LIBRARY= -DPETSC_LIBRARY= \
		-DCMAKE_SYSTEM_NAME=Cray -DCMAKE_C_COMPILER=cc

GNU & CRAY:

    cmake .. -DPETSC_INCLUDE_DIR= -DNETCDF_INCLUDE_DIR= \
        -DNETCDF_LIBRARY= -DPETSC_LIBRARY=



## Building COZOC

After CMake configuration, just type

    make

to build COZOC. In case of troubles, first try `make VERBOSE=1`.


## Spacemacs C IDE

The project includes `.spacemacs` configuration file that can be used
with [Spacemacs](http://spacemacs.org).

The above `cmake`-command also writes files `.clang_complete` and
`.dir-locals.el` files to the project root, that are used by the Spacemacs
packages. Thus, to get the full features of the Spacemacs C IDE working, run the
CMake configuration first.


## TODO:

- omegaQG.c was written before omega.c. Now omegaQC should be updated to
  same style and using the same subroutines as omega.
