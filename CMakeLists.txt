cmake_minimum_required (VERSION 3.2)

project (cozoc C)

list (APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/Modules
  ${CMAKE_SOURCE_DIR}/cmake/Functions)

enable_testing()

include(ExternalProject)

add_subdirectory (src)
add_subdirectory (test)
