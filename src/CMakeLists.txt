include (ide_support)
include (netcdf)

find_package (MPI)
find_package (PETSC)
find_package (NETCDF)

find_netcdf_par_h ()

set (CMAKE_C_COMPILER ${MPI_C_COMPILER})

include_directories (
  ${CMAKE_CURRENT_BINARY_DIR}
  ${PETSC_INCLUDE_DIR}
  ${NETCDF_INCLUDE_DIR} )

link_libraries (
  ${PETSC_LIBRARY}
  ${NETCDF_LIBRARY} )

configure_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

write_dir_locals ()
write_clang_complete_cflags ()

add_executable (cozoc constants.c context.c coriolis_parameter.c
  cozoc.c daslice.c derivatives.c io.c omega.c omegaQG.c
  sigma_parameter.c utils.c vorticity.c wrfnc.c field.c)

target_link_libraries (cozoc m)

add_custom_command( OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/omega_stencil.inc
  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/omega_stencil.py > ${CMAKE_CURRENT_SOURCE_DIR}/omega_stencil.inc
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/omega_stencil.py
  COMMENT "Generating source code" )
add_custom_target(omega_stencil
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/omega_stencil.inc)
add_dependencies(cozoc omega_stencil)
