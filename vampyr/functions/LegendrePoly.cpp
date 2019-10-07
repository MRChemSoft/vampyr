#include "pybind11/pybind11.h"

#include "PyRepresentableFunction.h"

#include "functions/LegendrePoly.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void legendre_poly(py::module &m, py::class_<RepresentableFunction<1>, PyRepresentableFunction<1>> &repfunc1d) {


py::class_<LegendrePoly>(m, "LegendrePoly", repfunc1d)
        .def(py::init<int, double, double>(), "k"_a, "n"_a = 1.0, "l"_a = 0.0)
        .def("firstDerivative", &LegendrePoly::firstDerivative)
        .def("secondDerivative", &LegendrePoly::secondDerivative);


}
} // namespace vampyr