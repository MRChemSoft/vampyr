# We save CMAKE_BUILD_TYPE, as we will set it to Release for externals
set(_build_type ${CMAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE Release)

find_package(MRCPP 1.4 CONFIG QUIET)

# whether MRCPP was fetched and built locally
set(MRCPP_FETCHED FALSE)

if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)

  get_target_property(MRCPP_HAS_OMP MRCPP::mrcpp MRCPP_HAS_OMP)
  get_target_property(MRCPP_HAS_MPI MRCPP::mrcpp MRCPP_HAS_MPI)

  message(STATUS "Found MRCPP (OpenMP: ${MRCPP_HAS_OMP}; MPI: ${MRCPP_HAS_MPI}): ${_loc} (found version ${MRCPP_VERSION})")

  if(NOT MRCPP_HAS_OMP)
    message(FATAL_ERROR "VAMPyR needs MRCPP with OpenMP features enabled!")
  endif()

  if(NOT MRCPP_HAS_MPI AND VAMPYR_WITH_MPI)
    message(FATAL_ERROR "You enabled MPI for VAMPyR, but MRCPP does not have MPI features enabled!")
  endif()
else()
  message(STATUS "Suitable MRCPP could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(mrcpp
    QUIET
    GIT_REPOSITORY
      https://github.com/MRChemSoft/mrcpp.git
    GIT_TAG
      1a6d450ac356feffba5c5cb2d02b92a875dbb9aa
    )

  set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  # Always build with OpenMP and without MPI
  set(ENABLE_OPENMP TRUE CACHE BOOL "")
  set(ENABLE_MPI ${VAMPYR_WITH_MPI} CACHE BOOL "")
  set(ENABLE_TESTS FALSE CACHE BOOL "")
  set(ENABLE_EXAMPLES FALSE CACHE BOOL "")

  FetchContent_MakeAvailable(mrcpp)
  set(MRCPP_FETCHED TRUE)
endif()

# reset CMAKE_BUILD_TYPE to whatever it was for VAMPyR
set(CMAKE_BUILD_TYPE ${_build_type})
