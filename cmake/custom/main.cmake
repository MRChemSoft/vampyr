set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
message(STATUS "${PROJECT_NAME} staged install: ${STAGED_INSTALL_PREFIX}")

if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_LIBDIR: python")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
endif()
file(TO_NATIVE_PATH "${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}/vampyr" PYMOD_INSTALL_FULLDIR)

add_subdirectory(external)
add_subdirectory(vampyr)
