include_directories (
  ${PROJECT_BINARY_DIR}/src
  ${PETSC_INCLUDE_DIRS}
  ${NETCDF_INCLUDE_DIRS})

link_directories (${PETSC_LIBRARY_DIRS} ${NETCDF_LIBRARY_DIRS})

link_libraries (${PETSC_LIBRARIES} ${NETCDF_LIBRARIES})

add_executable (cozoc arrays.c constants.c context.c coriolis_parameter.c
  cozoc.c daslice.c derivatives.c idtest.c io.c omega.c omegaQG.c
  sigma_parameter.c solve.c utils.c vecops.c vorticity.c wrfnc.c)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

function (write_clang_complete_cflags)
  get_directory_property (includes INCLUDE_DIRECTORIES)
  get_directory_property (definitions COMPILE_DEFINITIONS)
  if (includes)
    string (REPLACE ";" "\n-I" includes "-I${includes}")
  endif ()
  if (definitions)
    string (REPLACE ";" "\n-D" definitions "-I${definitions}")
  endif ()
  file (WRITE "${CMAKE_SOURCE_DIR}/.clang_complete" "${includes}\n${definitions}")
endfunction ()

write_clang_complete_cflags ()
