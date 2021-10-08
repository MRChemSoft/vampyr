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
  set(_RPATH "@loader_path/${_rel}")
else()
  set(_RPATH "\$ORIGIN/${_rel}")
endif()

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})

set_target_properties(_vampyr
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${_RPATH}lib"
    INSTALL_RPATH_USE_LINK_PATH ON
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
    RESOURCE "${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py"
  )

add_custom_command(
  TARGET 
    _vampyr 
  PRE_BUILD
  COMMAND 
    ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_LIST_DIR}/__init__.py $<TARGET_FILE_DIR:_vampyr>
  )

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
