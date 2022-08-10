# read PROGRAM_VERSION from file
if(EXISTS "${PROJECT_SOURCE_DIR}/VERSION")
  file(READ "${PROJECT_SOURCE_DIR}/VERSION" PROGRAM_VERSION)
  string(STRIP "${PROGRAM_VERSION}" PROGRAM_VERSION)
else()
  message(FATAL_ERROR "File ${PROJECT_SOURCE_DIR}/VERSION not found")
endif()

# define where to install the Python module
if(NOT DEFINED PYMOD_INSTALL_FULLDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_FULLDIR: lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/vampyr")
  file(TO_NATIVE_PATH "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/vampyr" PYMOD_INSTALL_FULLDIR)
else()
  message(STATUS "Setting option PYMOD_INSTALL_FULLDIR: ${PYMOD_INSTALL_FULLDIR}")
endif()

add_subdirectory(src)
