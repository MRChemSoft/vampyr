# we glob the Python test files and let CMake add a rule such that
# the glob is repeated every time we rebuild.
# This is NOT RECOMMENDED by CMake
# (https://cmake.org/cmake/help/v3.16/command/file.html#filesystem) but you only
# live once!
file(
   GLOB
    _vmp_pytests
  LIST_DIRECTORIES
    FALSE
  CONFIGURE_DEPENDS
  RELATIVE
    ${CMAKE_CURRENT_LIST_DIR}
  test_*.py
  )

foreach(_f IN LISTS _vmp_pytests)
  # link the Python files under the build folder
  file(
     CREATE_LINK
      ${CMAKE_CURRENT_LIST_DIR}/${_f}
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/tests/${_f}
    COPY_ON_ERROR
    SYMBOLIC
    )

  # extract the name and label for the test from the file name
  # by removing the "test_" prefix and the ".py" suffix
  string(REGEX REPLACE "^test_(.*).py$" "\\1" _name ${_f})
  add_test(
    NAME
      ${_name}
    COMMAND
      ${Python_EXECUTABLE} -m pytest -vv -s --pyargs vampyr -k test_${_name}
    WORKING_DIRECTORY
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/..
    )

  set_tests_properties(${_name}
    PROPERTIES
      LABELS
        "integration;${_name}"
    )
endforeach()

# install the current folder
install(
  DIRECTORY
    ${CMAKE_CURRENT_LIST_DIR}/
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}/tests
  PATTERN "CMakeLists.txt" EXCLUDE
  PATTERN "__pycache__" EXCLUDE
  )

