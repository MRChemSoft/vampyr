# minimum cmake version
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

# project name and supported language
project(account-project LANGUAGES CXX)

# require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# the pybind11 package will provide the pybind11::module imported target
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 3.3 CONFIG REQUIRED)
find_package(MRCPP CONFIG REQUIRED)

# create python module
add_library(vampyr
  MODULE
    vampyr.cpp
    PyBoundingBox.cpp
    bases.cpp
    fundamentalTypes.cpp
    methods.cpp
    operators.cpp
    project.cpp
  )

if(WIN32 AND (NOT MSVC))
  target_compile_definitions(vampyr
    PUBLIC
      "_hypot=hypot"
    )
endif()

target_link_libraries(vampyr
  PUBLIC
    pybind11::module
    Eigen3::Eigen
    MRCPP::mrcpp
  )

set_target_properties(vampyr
  PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
  )

file(
      COPY
        __init__.py
        version.py
      DESTINATION
        ${CMAKE_CURRENT_BINARY_DIR}
      )

install(
  TARGETS
    vampyr
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
