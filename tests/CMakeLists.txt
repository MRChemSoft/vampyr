set(_tests_in_build_tree ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/..)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/test_vampyr3d.py
  ${_tests_in_build_tree}/test_vampyr3d.py
  COPYONLY
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/pot_gen.py
  ${_tests_in_build_tree}/pot_gen.py
  COPYONLY
  )

# Define test
add_test(
  NAME
    vampyr3d
  COMMAND
    COMMAND ${PYTHON_EXECUTABLE} -m pytest -x -rws -vv ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/../test_vampyr3d.py
  )

install(
  FILES
    test_vampyr3d.py
    pot_gen.py
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}/..
  )
