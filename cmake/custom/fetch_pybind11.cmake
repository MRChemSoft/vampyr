find_package(pybind11 2.9 CONFIG QUIET)

set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")

if(pybind11_FOUND)
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
else()
  message(STATUS "Suitable pybind11 could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(pybind11
    QUIET
    URL
      https://github.com/pybind/pybind11/archive/v2.10.0.tar.gz
    )
  set(PYBIND11_TEST OFF CACHE BOOL "")
  FetchContent_MakeAvailable(pybind11)
endif()
