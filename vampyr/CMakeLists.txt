include(FetchContent)

FetchContent_Declare(
  pybind11_sources
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.2
)

FetchContent_GetProperties(pybind11_sources)

if(NOT pybind11_sources_POPULATED)
  FetchContent_Populate(pybind11_sources)

  add_subdirectory(
    ${pybind11_sources_SOURCE_DIR}
    ${pybind11_sources_BINARY_DIR}
    )
endif()

# Much disgusting
link_directories(${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

# create python module
add_library(vampyr3d
  MODULE
    vampyr3d.cpp
    bases.cpp
    project.cpp
  )

target_include_directories(vampyr3d
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

# This command will become unnecessary, eventually
target_include_directories(vampyr3d
  PUBLIC
    $<BUILD_INTERFACE:${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP>
  )

add_dependencies(vampyr3d
  mrcpp_external
  )

target_link_libraries(vampyr3d
  PUBLIC
    pybind11::module
    Eigen3::Eigen
    mrcpp
  )

set_target_properties(vampyr3d
  PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
  )

install(
  TARGETS
    vampyr3d
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  )
