# create python module
pybind11_add_module(vampyr3d MODULE THIN_LTO
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

file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${PYMOD_INSTALL_FULLDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(VAMPyR_RPATH "@loader_path/${_rel}")
else()
  set(VAMPyR_RPATH "\$ORIGIN/${_rel}")
endif()
set_target_properties(vampyr3d
  PROPERTIES
  MACOSX_RPATH ON
  SKIP_BUILD_RPATH OFF
  BUILD_WITH_INSTALL_RPATH OFF
  INSTALL_RPATH "${VAMPyR_RPATH}${CMAKE_INSTALL_LIBDIR}"
  INSTALL_RPATH_USE_LINK_PATH ON
  )
set_target_properties(vampyr3d
  PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
  RESOURCE
    "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
  )

foreach(_f __init__.py)
  configure_file(
    ${_f}
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/${_f}
    COPYONLY
    )
endforeach()

install(
  TARGETS
    vampyr3d
  LIBRARY
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  RESOURCE
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  )
