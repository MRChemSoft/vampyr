include(FetchContent)

set(PYBIND11_PYTHON_VERSION 3.5)
if(MSVC)
  set(PYBIND11_CPP_STANDARD "/std:c++${CMAKE_CXX_STANDARD}")
else()
  set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")
endif()

find_package(pybind11 2.3 CONFIG QUIET)
if(pybind11_FOUND)
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
else()
  message(STATUS "Suitable pybind11 could not be located. Fetching and building!")
  FetchContent_Declare(pybind11_sources
    QUIET
    URL
      https://github.com/pybind/pybind11/archive/v2.3.0.tar.gz
    )

  FetchContent_GetProperties(pybind11_sources)

  set(PYBIND11_PYTHON_VERSION ${PYBIND11_PYTHON_VERSION})
  set(PYBIND11_TEST OFF CACHE BOOL "")
  set(PYMOD_INSTALL_FULLDIR ${PYMOD_INSTALL_FULLDIR})

  if(NOT pybind11_sources_POPULATED)
    FetchContent_Populate(pybind11_sources)

    add_subdirectory(
      ${pybind11_sources_SOURCE_DIR}
      ${pybind11_sources_BINARY_DIR}
      )
  endif()
endif()

set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
message(STATUS "${PROJECT_NAME} staged install: ${STAGED_INSTALL_PREFIX}")

find_package(MRCPP 1.0 CONFIG QUIET)
if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
  message(STATUS "Found MRCPP: ${_loc} (found version ${MRCPP_VERSION})")
else()
  message(STATUS "Suitable MRCPP could not be located. Fetching and building!")
  FetchContent_Declare(mrcpp_sources
    QUIET
    URL
      https://github.com/MRChemSoft/mrcpp/archive/master.tar.gz
    )

  FetchContent_GetProperties(mrcpp_sources)

  set(CMAKE_INSTALL_PREFIX ${STAGED_INSTALL_PREFIX})
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
  set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(ENABLE_OPENMP ${ENABLE_OPENMP} CACHE BOOL "")
  set(ENABLE_MPI ${ENABLE_MPI} CACHE BOOL "")
  set(ENABLE_TESTS FALSE CACHE BOOL "")
  set(ENABLE_EXAMPLES FALSE CACHE BOOL "")
  set(SHARED_LIBRARY_ONLY TRUE CACHE BOOL "")

  if(NOT mrcpp_sources_POPULATED)
    FetchContent_Populate(mrcpp_sources)

    add_subdirectory(
      ${mrcpp_sources_SOURCE_DIR}
      ${mrcpp_sources_BINARY_DIR}
      )
  endif()
endif()
