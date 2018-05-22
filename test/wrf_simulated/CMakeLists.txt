configure_file(prepare-tests.make.in prepare-tests.make)
add_test (NAME cozoc
  COMMAND make -f prepare-tests.make)
add_custom_target (check-wrf-simulated
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS cozoc)
