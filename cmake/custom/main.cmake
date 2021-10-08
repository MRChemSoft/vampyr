# read PROGRAM_VERSION from file
if(EXISTS "${PROJECT_SOURCE_DIR}/VERSION")
  file(READ "${PROJECT_SOURCE_DIR}/VERSION" PROGRAM_VERSION)
  string(STRIP "${PROGRAM_VERSION}" PROGRAM_VERSION)
else()
  message(FATAL_ERROR "File ${PROJECT_SOURCE_DIR}/VERSION not found")
endif()

# define where to install the Python module
if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_LIBDIR: python")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location within lib to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location within lib to which Python modules are installed" FORCE)
endif()
file(TO_NATIVE_PATH "lib/${PYMOD_INSTALL_LIBDIR}/vampyr" PYMOD_INSTALL_FULLDIR)

add_subdirectory(src)
