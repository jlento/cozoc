include (LibFindMacros)

libfind_pkg_check_modules (NETCDF_PKG_CONFIG netcdf)

find_path(NETCDF_INCLUDE_DIR
  NAMES netcdf_par.h
  HINTS
  $ENV{NETCDF_DIR}/include
  ${CMAKE_CURRENT_BINARY_DIR}/../netcdf/include
  ${NETCDF_PKG_CONFIG_INCLUDE_DIRS})

find_library(NETCDF_LIBRARY
  NAMES netcdf
  HINTS
  $ENV{NETCDF_DIR}/lib
  ${CMAKE_CURRENT_BINARY_DIR}/../netcdf/lib
  ${NETCDF_PKG_CONFIG_LIBRARY_DIRS})
set (NETCDF_PROCESS_INCLUDES
  NETCDF_INCLUDE_DIR)
set (NETCDF_PROCESS_LIBS
  NETCDF_LIBRARY)
libfind_process (NETCDF)
