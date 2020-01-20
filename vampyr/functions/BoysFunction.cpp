#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>

#include <MRCPP/functions/BoysFunction.h>

#include "PyRepresentableFunction.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {
void boys_function(py::module &m) {
    py::class_<BoysFunction, RepresentableFunction<1>>(m, "BoysFunction")
        .def(py::init<int, double>(), "nTerms"_a = 0, "Boys_func_prec"_a = 1.0e-10);
}
} // namespace vampyr
