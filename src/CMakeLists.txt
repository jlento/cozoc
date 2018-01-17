include (ide_support)
include (AlwaysRPATH)

find_package (MPI REQUIRED)
find_package (PETSC REQUIRED)
find_package (NETCDF REQUIRED)

set (CMAKE_C_COMPILER ${MPI_C_COMPILER})

# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

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
  cozoc.c daslice.c ops.c io.c omega.c omegaQG.c)

target_link_libraries (cozoc m)

install (TARGETS cozoc DESTINATION bin)
