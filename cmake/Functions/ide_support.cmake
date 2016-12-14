function (write_clang_complete_cflags)
  message (STATUS "Writing ${CMAKE_SOURCE_DIR}/.clang_complete")
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


function (write_dir_locals)
  message (STATUS "Writing ${CMAKE_SOURCE_DIR}/.dir-locals.el")
  file (WRITE "${CMAKE_SOURCE_DIR}/.dir-locals.el"
    "((c-mode (helm-make-build-dir . \"${PROJECT_BINARY_DIR}\")))\n")
endfunction ()
