include (ide_support)
include (AlwaysRPATH)
include (BuildType)

find_package (MPI REQUIRED)
find_package (PETSC REQUIRED)
find_package (NETCDF REQUIRED)

set (CMAKE_C_COMPILER ${MPI_C_COMPILER})

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra -pedantic -Wformat=2 -Wno-unused-parameter -Wshadow -Wwrite-strings -Wstrict-prototypes -Wold-style-definition -Wredundant-decls -Wnested-externs -Wmissing-include-dirs -Wjump-misses-init -Wlogical-op")
  endif ()
endif ()

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${PETSC_INCLUDE_DIR}
  ${NETCDF_INCLUDE_DIR} )

link_libraries (
  ${PETSC_LIBRARY}
  ${NETCDF_LIBRARY} )

write_dir_locals ()
write_clang_complete_cflags ()

add_executable (cozoc constants.c options.c context.c
  equation.c cozoc.c daslice.c ops.c io.c omega.c omegaQG.c)

target_link_libraries (cozoc m)

install (TARGETS cozoc DESTINATION bin)
