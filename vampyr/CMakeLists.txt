include(FetchContent) 

FetchContent_Declare(
    pybind11_sources
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.2
)

FetchContent_Declare(
    mrcpp_sources
    GIT_REPOSITORY https://github.com/bjorgve/mrcpp.git
    GIT_TAG        mrcpp_evalf 
)

FetchContent_GetProperties(pybind11_sources)

if(NOT pybind11_sources_POPULATED)
    FetchContent_Populate(pybind11_sources)

    add_subdirectory(
        ${pybind11_sources_SOURCE_DIR}
        ${pybind11_sources_BINARY_DIR}

        )
endif()

FetchContent_GetProperties(mrcpp_sources)

if(NOT mrcpp_sources_POPULATED)
    FetchContent_Populate(mrcpp_sources)

    add_subdirectory(
        ${mrcpp_sources_SOURCE_DIR}
        ${mrcpp_sources_BINARY_DIR}

        )
endif()

add_library(vampyr3d
  MODULE
    vampyr3d.cpp
    PyBoundingBox.cpp
    bases.cpp
    project.cpp
  )

target_link_libraries(vampyr3d
  PUBLIC
    pybind11::module
    MRCPP::mrcpp
  )

set_target_properties(vampyr3d
  PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
  )

