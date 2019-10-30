#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "PyRepresentableFunction.h"

#include "functions/Polynomial.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void polynomial(py::module &m, py::class_<RepresentableFunction<1>, PyRepresentableFunction<1>> &repfunc1d) {

// The double * wont work properly in vampyr since raw pointers does not exist in python. 
// This should probably be changed to an std::vector in mrcpp. 
py::class_<Polynomial>(m, "Polynomial", repfunc1d)
        .def(py::init<int, const double *, const double *>(), "order"_a, "a"_a = nullptr, "b"_a = nullptr)
        .def("evalf", py::overload_cast<double>(&Polynomial::evalf, py::const_))
        .def("setCoefs", &Polynomial::setCoefs)
        .def("normalize", &Polynomial::normalize)
        .def("getCoefs", py::overload_cast<>(&Polynomial::getCoefs))
        .def("size", &Polynomial::size);


}
} // namespace vampyr
