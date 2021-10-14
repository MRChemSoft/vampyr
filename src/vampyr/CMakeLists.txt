# create python module
pybind11_add_module(_vampyr MODULE THIN_LTO
    export_vampyr.cpp
  )

target_compile_definitions(_vampyr
  PRIVATE
    VERSION_INFO="${PROGRAM_VERSION}"
  )

target_include_directories(_vampyr
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )

target_link_libraries(_vampyr
  PUBLIC
    MRCPP::mrcpp
  )

file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${PYMOD_INSTALL_FULLDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(_RPATH "@loader_path/${_rel}/lib;@loader_path/${CMAKE_INSTALL_LIBDIR}")
else()
  set(_RPATH "\$ORIGIN/${_rel}/lib;\$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
endif()

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})

# we glob any Python file present and let CMake add a rule such that
# the glob is repeated every time we rebuild.
# This is NOT RECOMMENDED by CMake
# (https://cmake.org/cmake/help/v3.16/command/file.html#filesystem) but you only
# live once!
file(
   GLOB
    _vmp_pys
  LIST_DIRECTORIES
    FALSE
  CONFIGURE_DEPENDS
  RELATIVE
    ${CMAKE_CURRENT_LIST_DIR}
  *.py
  )

set_target_properties(_vampyr
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${_RPATH}"
    INSTALL_RPATH_USE_LINK_PATH ON
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
    RESOURCE "${_vmp_pys}"
  )

foreach(_f IN LISTS _vmp_pys)
  # link the Python files under the build folder
  file(
     CREATE_LINK
      ${CMAKE_CURRENT_LIST_DIR}/${_f}
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/${_f}
    COPY_ON_ERROR
    SYMBOLIC
    )
endforeach()

install(
  TARGETS
    _vampyr
  LIBRARY
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  RESOURCE
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  )
