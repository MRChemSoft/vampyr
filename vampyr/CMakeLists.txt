# create python module
pybind11_add_module(vampyr3d MODULE THIN_LTO "")

target_sources(vampyr3d
  PRIVATE
    vampyr3d.cpp
    bases.cpp
    project.cpp
    math.cpp
    world.cpp
    trees.cpp
    derivative.cpp
    applys.cpp
    representable_functions.cpp
  )

target_link_libraries(vampyr3d
  PUBLIC
    Eigen3::Eigen
    MRCPP::mrcpp
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
