# We save CMAKE_BUILD_TYPE, as we will set it to Release for externals
set(_build_type ${CMAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE Release)

find_package(MRCPP 1.4 CONFIG QUIET)

# whether MRCPP was fetched and built locally
set(MRCPP_FETCHED FALSE)

if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
  message(STATUS "Found MRCPP: ${_loc} (found version ${MRCPP_VERSION})")
else()
  message(STATUS "Suitable MRCPP could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(mrcpp
    QUIET
    GIT_REPOSITORY
      https://github.com/MRChemSoft/mrcpp.git
    GIT_TAG
      e2f0e7e897261dbf
    )

  set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  # Always build with OpenMP and without MPI
  set(ENABLE_OPENMP TRUE CACHE BOOL "")
  set(ENABLE_MPI FALSE CACHE BOOL "")
  set(ENABLE_TESTS FALSE CACHE BOOL "")
  set(ENABLE_EXAMPLES FALSE CACHE BOOL "")

  FetchContent_MakeAvailable(mrcpp)
  set(MRCPP_FETCHED TRUE)
endif()

# reset CMAKE_BUILD_TYPE to whatever it was for VAMPyR
set(CMAKE_BUILD_TYPE ${_build_type})
