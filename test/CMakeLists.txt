configure_file(CTestCustom.cmake CTestCustom.cmake)
configure_file(prepare-tests.make.in prepare-tests.make)
add_test (NAME cozoc
  COMMAND cozoc -G -f test-input-wrf.nc4)
add_custom_target (check
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS cozoc)

