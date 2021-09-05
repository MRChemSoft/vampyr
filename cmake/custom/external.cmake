# We save CMAKE_BUILD_TYPE, as we will set it to Release for externals
set(_build_type ${CMAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE Release)

include(FetchContent)

set(PYBIND11_PYTHON_VERSION 3.6)
set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")

find_package(pybind11 2.5 CONFIG QUIET)
if(pybind11_FOUND)
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
else()
  message(STATUS "Suitable pybind11 could not be located. Fetching and building!")
  FetchContent_Declare(pybind11_sources
    QUIET
    URL
      https://github.com/pybind/pybind11/archive/v2.5.0.tar.gz
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

find_package(MRCPP 1.3.5 CONFIG QUIET)
if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
  message(STATUS "Found MRCPP: ${_loc} (found version ${MRCPP_VERSION})")
else()
  message(STATUS "Suitable MRCPP could not be located. Fetching and building!")
  FetchContent_Declare(mrcpp_sources
    QUIET
    GIT_REPOSITORY
      https://github.com/stigrj/mrcpp.git
    GIT_TAG
      06aaa3f81ffbc5fa3eb2f250838ee5f6d0ed89d7
    )

  FetchContent_GetProperties(mrcpp_sources)

  set(CMAKE_INSTALL_PREFIX ${STAGED_INSTALL_PREFIX})
  set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  # Always build with OpenMP and without MPI
  set(ENABLE_OPENMP TRUE CACHE BOOL "")
  set(ENABLE_MPI FALSE CACHE BOOL "")
  set(ENABLE_TESTS FALSE CACHE BOOL "")
  set(ENABLE_EXAMPLES FALSE CACHE BOOL "")

  if(NOT mrcpp_sources_POPULATED)
    FetchContent_Populate(mrcpp_sources)

    add_subdirectory(
      ${mrcpp_sources_SOURCE_DIR}
      ${mrcpp_sources_BINARY_DIR}
      )
  endif()
endif()

# reset CMAKE_BUILD_TYPE to whatever it was for VAMPyR
set(CMAKE_BUILD_TYPE ${_build_type})
