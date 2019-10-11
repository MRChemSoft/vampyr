# create python module
pybind11_add_module(vampyr3d MODULE THIN_LTO "")

target_sources(vampyr3d
  PRIVATE
    vampyr3d.cpp
  )
target_include_directories(vampyr3d
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )


target_link_libraries(vampyr3d
  PUBLIC
    Eigen3::Eigen
    MRCPP::mrcpp
  )


add_subdirectory(applys)
add_subdirectory(bases)
add_subdirectory(functions)
add_subdirectory(math)
add_subdirectory(operators)
add_subdirectory(project)
add_subdirectory(world)
add_subdirectory(trees)


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
