# define minimum cmake version
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

# project name and supported language
project(VAMPYR LANGUAGES CXX)

# require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# for testing we will need the python interpreter
find_package(PythonInterp REQUIRED)

# we need ExternalProject functionality
include(ExternalProject)

set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX})

# fetch pybind11 sources
ExternalProject_Add(pybind11_external
  GIT_REPOSITORY
    https://github.com/pybind/pybind11.git
  GIT_TAG
    v2.2
  CMAKE_ARGS
    -DPYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}  # unused, but needs working compiler
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DPYBIND11_TEST=OFF
  LOG_DOWNLOAD 1
  LOG_UPDATE 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  )
set(PYBIND11_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/pybind11 CACHE PATH "Path to internally built pybind11Config.cmake" FORCE)

# fetch mrcpp sources
ExternalProject_Add(mrcpp_external
  GIT_REPOSITORY
    https://github.com/MRChemSoft/mrcpp.git
  GIT_TAG
    master
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DENABLE_TESTS=FALSE
    -DENABLE_EXAMPLES=FALSE
  CMAKE_CACHE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  )
set(MRCPP_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/MRCPP CACHE PATH "Path to internally built MRCPPConfig.cmake" FORCE)

# process our source code as if it was an external project
ExternalProject_Add(vampyr
  DEPENDS pybind11_external
  DEPENDS mrcpp_external
  SOURCE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/vampyr
  BINARY_DIR
    ${CMAKE_CURRENT_BINARY_DIR}/vampyr
  CMAKE_ARGS
    -Dpybind11_DIR=${PYBIND11_DIR}
    -DMRCPP_DIR=${MRCPP_DIR}
    -DPYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  )

# turn on testing
enable_testing()

# define test
add_test(
  NAME
    python_test
  COMMAND
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vampyr/test_vampyr1d.py
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vampyr/test_vampyr2d.py
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vampyr/test_vampyr3d.py
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )

# test needs to know the location of the compiled module
set_tests_properties(python_test
  PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PYTHONPATH}
  )
