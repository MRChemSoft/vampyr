add_custom_target(git-submodule-init
    COMMAND git submodule init
    COMMAND git submodule update
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

find_package(MRCPP CONFIG QUIET)
if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
  message(STATUS "Found MRCPP: ${_loc} (found version ${MRCPP_VERSION})")
  add_library(mrcpp_external INTERFACE)  # dummy
else()
  message(STATUS "Suitable MRCPP could not be located. Building the local git submodule!")
  include(ExternalProject)
  ExternalProject_Add(mrcpp_external
    DOWNLOAD_COMMAND
      ""
    SOURCE_DIR
      ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DENABLE_OPENMP=${ENABLE_OPENMP}
      -DENABLE_MPI=${ENABLE_MPI}
      -DEigen3_DIR=${Eigen3_DIR}
      -DENABLE_TESTS=TRUE
      -DENABLE_EXAMPLES=TRUE
    CMAKE_CACHE_ARGS
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    TEST_BEFORE_INSTALL
      1
    )
  set(MRCPP_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/MRCPP CACHE PATH "Path to internally built MRCPPConfig.cmake" FORCE)
endif()

install(
  DIRECTORY
    ${STAGED_INSTALL_PREFIX}/
  DESTINATION
    .
  USE_SOURCE_PERMISSIONS
  PATTERN
    "*.pyc" EXCLUDE
  )
