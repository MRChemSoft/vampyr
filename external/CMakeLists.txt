include(FetchContent)

find_package(pybind11 2.2 CONFIG QUIET)
if(pybind11_FOUND)
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
else()
  message(STATUS "Suitable pybind11 could not be located. Fetching and building!")
  FetchContent_Declare(
    pybind11_sources
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.2.4
  )

  FetchContent_GetProperties(pybind11_sources)

  if(NOT pybind11_sources_POPULATED)
    FetchContent_Populate(pybind11_sources)

    add_subdirectory(
      ${pybind11_sources_SOURCE_DIR}
      ${pybind11_sources_BINARY_DIR}
      )
  endif()
endif()

find_package(MRCPP CONFIG QUIET)
if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
  message(STATUS "Found MRCPP: ${_loc} (found version ${MRCPP_VERSION})")
else()
  message(STATUS "Suitable MRCPP could not be located. Fetching and building!")

    FetchContent_Populate(mrcpp_sources
      QUIET
      GIT_REPOSITORY
        https://github.com/MRChemSoft/mrcpp.git
      GIT_TAG
        master # Preferable to have a tag for a release
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DENABLE_OPENMP=${ENABLE_OPENMP}
        -DENABLE_MPI=${ENABLE_MPI}
        -DEigen3_DIR=${Eigen3_DIR}
        -DENABLE_TESTS=TRUE
        -DENABLE_EXAMPLES=FALSE
        -DSHARED_LIBRARY_ONLY=TRUE # Seems to be ignored...
      CMAKE_CACHE_ARGS
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      TEST_BEFORE_INSTALL # Also this one seems to be ignored...
        1
      )

    add_subdirectory(
      ${mrcpp_sources_SOURCE_DIR}
      ${mrcpp_sources_BINARY_DIR}
      )

    add_library(MRCPP::mrcpp ALIAS mrcpp-shared)
endif()
