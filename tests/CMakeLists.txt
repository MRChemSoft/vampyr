configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_vampyr3d.py ${CMAKE_CURRENT_BINARY_DIR})

# Define test
add_test(
  NAME
    vampyr3d
  COMMAND
    ${CMAKE_COMMAND} -E env VAMPYR_MODULE_PATH=$<TARGET_FILE_DIR:vampyr3d>
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/test_vampyr3d.py
  )
