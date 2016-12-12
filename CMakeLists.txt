cmake_minimum_required (VERSION 3.2)
project (cozoc C)

message ("
COZOC's cmake configuration uses pkg-config to determine build flags for
PETSc and Netcdf libraries. Set environment variable PKG_CONFIG_PATH to
control where the libraries are searched first, with
'export PKG_CONFIG_PATH=\"\${NETCDF_DIR}/lib/pkgconfig\"', for example.
")

find_package (PkgConfig)
pkg_check_modules (PKG_CONFIG_PETSC PETSc netcdf)


find_file (NETCDF_PAR_H netcdf_par.h
  DOC "Parallel Netcdf4/HDF5 has netcdf_par.h include file"
  HINTS ${PKG_CONFIG_NETCDF_INCLUDE_DIRS}
  )
if (NETCDF_PAR_H OR USE_PARALLEL_NETCDF)
  message (STATUS "Using parallel netcdf")
  message (STATUS "  Override with '-DUSE_PARALLEL_NETCDF=OFF'")
else ()
  message (STATUS "Using sequential netcdf.")
  message (STATUS "  Override with '-DUSE_PARALLEL_NETCDF=ON'")
endif ()

#if (NOT NETCDF_PAR_H_FOUND)
#  include (ExternalProject)
#  pkg_check_modules (PKG_CONFIG_HDF5 hdf5)
#  ExternalProject_Add (netcdf
#    EXCLUDE_FROM_ALL 1
#    URL ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.tar.gz
#    CONFIGURE_COMMAND eval ../netcdf/configure --prefix=${PROJECT_BINARY_DIR} --disable-shared -disable-dap -disable-utilities CC=${CMAKE_C_COMPILER} CFLAGS="${PKG_CONFIG_HDF5_CFLAGS}" LDFLAGS="${PKG_CONFIG_HDF5_LDFLAGS}"
#    INSTALL_COMMAND make install)
#endif ()

#    GIT_REPOSITORY https://github.com/Unidata/netcdf-c.git
#    PREFIX "${PROJECT_BINARY_DIR}"
#    CMAKE_CACHE_ARGS -DENABLE_DAP:STRING=OFF -DBUILD_UTILITIES:STRING=OFF -DBUILD_SHARED:STRING=OFF -DENABLE_TESTS:STRING=OFF -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}

# CMAKE_C_COMPILER cache default needs to be set before project is declared

#pkg_get_variable(PKG_CONFIG_PETSC_C_COMPILER PETSc ccompiler)
#set (CMAKE_C_COMPILER ${PKG_CONFIG_PETSC_C_COMPILER} CACHE FILEPATH "C compiler")

# Set cache defaults for dependencies

set (PETSC_INCLUDE_DIRS ${PKG_CONFIG_PETSC_INCLUDE_DIRS} CACHE STRING "PETSc include directories (w/o the '-I')")
set (PETSC_LIBRARY_DIRS ${PKG_CONFIG_PETSC_LIBRARY_DIRS} CACHE STRING "PETSc library paths (w/o the '-L')")
set (PETSC_LIBRARIES ${PKG_CONFIG_PETSC_LIBRARIES} CACHE STRING "PETSc libraries (w/o the '-l')")

#set (HDF5_CFLAGS ${PKG_CONFIG_HDF5_CFLAGS_S} CACHE STRING "HDF5 include flags (used only if netcdf4 build)")
#set (HDF5_LDFLAGS ${PKG_CONFIG_HDF5_LDFLAGS_S} CACHE STRING "HDF5 linker flags (used only if netcdf4 build)")

set (NETCDF_INCLUDE_DIRS ${PKG_CONFIG_NETCDF_INCLUDE_DIRS} CACHE STRING "PETSc include directories (w/o the '-I')")
set (NETCDF_LIBRARY_DIRS ${PKG_CONFIG_NETCDF_LIBRARY_DIRS} CACHE STRING "PETSc library paths (w/o the '-L')")
set (NETCDF_LIBRARIES ${PKG_CONFIG_NETCDF_LIBRARIES} CACHE STRING "PETSc libraries (w/o the '-l')")

if (NETCDF_PAR_H)
  set (USE_PARALLEL_NETCDF ON CACHE BOOL "Use parallel netcdf4/HDF5")
endif ()

add_subdirectory (src)
