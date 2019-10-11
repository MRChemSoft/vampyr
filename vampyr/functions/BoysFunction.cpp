#include "pybind11/pybind11.h"

#include "PyRepresentableFunction.h"

#include "functions/BoysFunction.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void boys_function(py::module &m, py::class_<RepresentableFunction<3>, PyRepresentableFunction<3>> &repfunc) {
    const auto D = 3;

py::class_<BoysFunction>(m, "BoysFunction", repfunc)
        .def(py::init<int, double>(), "nTerms"_a = 0, "Boys_func_prec"_a = 1.0e-10)
        .def("evalf", py::overload_cast<const Coord<1> &>(&BoysFunction::evalf, py::const_));


}
} // namespace vampyr
