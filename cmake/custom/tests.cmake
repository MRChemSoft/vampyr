# Turn on testing
enable_testing()

# Define test
add_test(
  NAME
    python_test
  COMMAND
    ${CMAKE_COMMAND} -E env VAMPYR_MODULE_PATH=$<TARGET_FILE_DIR:vampyr>
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/vampyr/test_vampyr3d.py
  )

# This must come last!!
#add_subdirectory(tests)
