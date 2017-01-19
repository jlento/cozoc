include (LibFindMacros)

libfind_pkg_check_modules (PETSC_PKG_CONFIG PETSc)

find_path(PETSC_INCLUDE_DIR
  NAMES petscsys.h
  HINTS $ENV{PETSC_DIR}/include ${PETSC_PKG_CONFIG_INCLUDE_DIRS})

if (PETSC_PKG_CONFIG_FOUND)
  set (petsc_libs ${PETSC_PKG_CONFIG_LIBRARIES})
else ()
  set (petsc_libs petsc)
endif ()

find_library(PETSC_LIBRARY
  NAMES ${petsc_libs}
  HINTS $ENV{PETSC_DIR}/lib ${PETSC_PKG_CONFIG_LIBRARY_DIRS})

set (PETSC_PROCESS_INCLUDES
  PETSC_INCLUDE_DIR)
set (PETSC_PROCESS_LIBS
  PETSC_LIBRARY)
libfind_process (PETSC)
