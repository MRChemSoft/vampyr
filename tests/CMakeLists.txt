set(_tests_in_build_tree ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/..)

list(APPEND _tests
  basis
  mra_1d
  mra_3d
  gauss_1d
  gauss_3d
  trees_1d
  trees_3d
  grids_and_project_1d
  grids_and_project_3d
  arithmetics_1d
  arithmetics_3d
  derivatives_1d
  derivatives_3d
  )

foreach(_test IN LISTS _tests)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_${_test}.py
    ${_tests_in_build_tree}/test_${_test}.py
    COPYONLY
    )
  add_test(
    NAME
      ${_test}
    COMMAND
      COMMAND ${PYTHON_EXECUTABLE} -m pytest -x -rws -vv ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/../test_${_test}.py
    )
  install(
    FILES
      test_${_test}.py
    DESTINATION
      ${PYMOD_INSTALL_FULLDIR}/..
    )
endforeach()
