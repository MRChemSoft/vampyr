find_package (Python COMPONENTS Interpreter Development)

# Pybind11 uses upper case version of these variables
# while cmake 3.12 uses mixed cases
# Note: We must use cmake 3.12 because it can
# find python more robustly
set(PYTHONLIBS_FOUND ${Python_FOUND})
set(PYTHON_INCLUDE_DIRS ${Python_INCLUDE_DIRS})
set(PYTHON_LIBRARIES ${Python_LIBRARIES})
