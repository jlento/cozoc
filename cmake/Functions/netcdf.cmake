function (find_netcdf_par_h)

  find_path(NETCDF_PAR_INCLUDE_DIR
    NAMES netcdf_par.h
    PATHS ${NETCDF_INCLUDE_DIRS})

  if (NETCDF_PAR_INCLUDE_DIR)
    message (STATUS "  Found netcdf_par.h")
  else ()
    message (STATUS "  Did not find netcdf_par.h")
    set (USE_PARALLEL_NETCDF OFF CACHE BOOL "Use parallel NetCDF4/HDF5")
  endif ()
  unset (NETCDF_PAR_INCLUDE_DIR CACHE)

  if (USE_PARALLEL_NETCDF)
    message (STATUS "  Using parallel netcdf")
  else ()
    message (STATUS "  Using sequential netcdf.")
  endif ()

endfunction ()


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
