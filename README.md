# COZOC

Cozoc is parallel re-implementation of [OZO](https://github.com/mikarant/ozo),
described in the paper [OZO v.1.0](https://doi.org/10.5194/gmd-10-827-2017).

## Status

- Computes the diabatic heating component of the generalize omega equation
- very much under construction...


## Obtaining COZOC source

Clone the COZOC GitHub repository, with for example:

    mkdir -p ~/github
    git clone https://jlento@bitbucket.org/jlento/cozoc.git ~/github/cozoc


## COZOC requirements

Building COZOC requires C compiler, CMake, and the following
libraries:

- MPI
- PETSc
- HDF5 (parallel)
- NetCDF4/HDF5 (parallel)


### (WORK IN PROGRESS, DOES NOT WORK PROPERLY, YET!!!") Setting up the environment for COZOC using conda package manager

After cloning the cozoc repository, install the requirements and development
tools by running

    bash setup.sh

This takes a while as it is installing the complete software stack for the
dependencies. Tested on Ubuntu 16.04 LTS, and should work on many other Linux
OSes, too. Likely will work on OS X with minor modifications.

Then, setup the environment with

    source env.sh

If everything went fine, one should be able to build and test cozoc with

    mkdir build
    cd build
    cmake ..
    make check
    

### Setting up the environment for COZOC in Ubuntu 16.04 LTS

The file [playbook.yml](playbook.yml) should give a good idea which additional
packages need to be present.

At the time of writing, Ubuntu/Debian does have a package for the parallel
version of HDF5, but not for NetCDF4, which needs to be compiled separately from
the sources. If parallel version of NetCDF4/HDF5 library needs to be build from
the sources, you can use the provided [netcdf4.bash](netcdf4.bash) script, for
example:

    mkdir -p ~/thirdparty
    INSTALLDIR=~/thirdparty bash thirdparty/build-netcdf4.sh


## Configure COZOC using CMake

COZOC build is configured using CMake.

It is a good practice to have separate source and build directories
(out-of-source builds). Let's create a build directory and continue in it, with
for example,

    mkdir $HOME/build
    cd $_


### In Ubuntu 16.04 LTS

Run cmake in the build directory, with for example

    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD \
             -DCMAKE_PREFIX_PATH=~/thirdparty

The CMake variable `CMAKE_INSTALL_PREFIX` sets the install root, and
`CMAKE_PREFIX_PATH` helps cmake to find libraries that are not in the standard
system directories.


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


## Tests

The provided tests serve as examples on how to use COZOC in general. To run just
a single test, for example the quickest with simulated input file, run

    make check-wrf-simulated

The available tests can listed with

    ctest -N
    
The whole test suite can be run with

    ctest


### WRF baroclinic wave

An example script to build WRF model is TODO




## Preparing the input file

Cozoc reads the input fields for the computation of the generalized omega
equations (6) - (11) in [OZO v.1.0](https://doi.org/10.5194/gmd-10-827-2017)
from a netcdf file. The definitions of the expected input fields in the netcdf
file can be read from preprocessing scripts for different models in directory
[preprocess](preprocess). We have chosen to keep the pre-processing of the model
input data to separate.







Currently COZOC can be used to analyse WRF simulations run on rectangular
grid, such as the baroclinic wave test case. Before the WRF output can be
processed by COZOC,

1. the fields need to be interpolated to the pressure levels

    TODO: `<wrfinterp command here>`

2. the file needs to be converted to NetCDF4/HDF5 format

    `nccopy -k nc4 <input> <output>`


## Running COZOC

In Ubuntu 16.04 LTS:

    # If running on machine without Infiniband network
    mkdir ~/.openmpi
    echo "btl_base_exclude = openib" >> ~/.openmpi/mca-params.config
    
    mpiexec -n 3 -mca btl ^openib ./src/cozoc -Q -G

In Cray XC40:

    aprun -n 3 ./src/cozoc -Q -G


## Profiling in Cray XC40 (...rough overview)

Load performance profiling tools before building cozoc

    module load perftools-base perftools

and instrument (for sampling profile) code after build

    pat_build src/cozoc

Run

    aprun -n 4 ./build/cozoc+pat -G

Read the profile

    pat_report cozoc+path+*

To instrument the code for MPI message timeline profiling

    path_build -g mpi,netcdf

and run the tracing experiment

    aprun -e PAT_RT_SUMMARY=0 -n 4 ./cozoc+pat -G

and view results

    pat_report *.xf
    app2 *.ap2


## Spacemacs C IDE

The project includes `.spacemacs` configuration file that can be used
with [Spacemacs](http://spacemacs.org).

The above `cmake`-command also writes files `.clang_complete` and
`.dir-locals.el` files to the project root, that are used by the Spacemacs
packages. Thus, to get the full features of the Spacemacs C IDE working, run the
CMake configuration first.

### Requirements

- clang
- clang-format
- astyle
-spacemacs

    git clone https://github.com/syl20bnr/spacemacs ./.emacs.d

### Use

    HOME=/directory/where/.spacemacs/and/.emacs.s/are emacs &


## Coding style

PETSc library functions are in `CamelCase`, COZOC application functions are in
`snake_case`. Source code formatting is influenced more by Python and Lisp,
where the curly braces have less visual emphasis, than the usual C styles.
