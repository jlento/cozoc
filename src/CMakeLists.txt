include_directories (${PETSC_INCLUDE_DIRS})
link_directories (${PETSC_LIBRARY_DIRS})
link_libraries (${PETSC_LIBRARIES})

add_executable (cozoc arrays.c constants.c context.c coriolis_parameter.c cozoc.c daslice.c derivatives.c idtest.c io.c omega.c omegaQG.c sigma_parameter.c solve.c utils.c vecops.c vorticity.c wrfnc.c)
