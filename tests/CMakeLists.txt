set(_tests_in_build_tree ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/..)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/test_gauss.py
  ${_tests_in_build_tree}/test_gauss.py
  COPYONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/test_mra.py
  ${_tests_in_build_tree}/test_mra.py
  COPYONLY
  )

# Define tests
add_test(
  NAME
    mra
  COMMAND
    COMMAND ${PYTHON_EXECUTABLE} -m pytest -x -rws -vv ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/../test_mra.py
  )
add_test(
  NAME
    gauss
  COMMAND
    COMMAND ${PYTHON_EXECUTABLE} -m pytest -x -rws -vv ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/../test_gauss.py
  )

install(
  FILES
    test_mra.py
    test_gauss.py
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}/..
  )
