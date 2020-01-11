#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "PyRepresentableFunction.h"

#include <MRCPP/functions/Polynomial.h>

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {
void polynomial(py::module &m) {
    py::class_<Polynomial, RepresentableFunction<1>>(m, "Polynomial")
        .def(py::init<int, const std::vector<double> &, const std::vector<double> &>(),
             "order"_a = 0,
             "a"_a = std::vector<double>(),
             "b"_a = std::vector<double>())
        .def(py::init<const Eigen::VectorXd &, const std::vector<double> &, const std::vector<double> &>(),
             "c"_a,
             "a"_a = std::vector<double>(),
             "b"_a = std::vector<double>())
        .def(py::init<double, int, const std::vector<double> &, const std::vector<double> &>(),
             "c"_a,
             "k"_a = 0,
             "a"_a = std::vector<double>(),
             "b"_a = std::vector<double>())
        .def("evalf", py::overload_cast<double>(&Polynomial::evalf, py::const_))
      //.def("evalf", py::overload_cast<const Coord<1> &>(&Polynomial::evalf, py::const_))
        .def("setCoefs", &Polynomial::setCoefs)
        .def("normalize", &Polynomial::normalize)
        .def("getCoefs", py::overload_cast<>(&Polynomial::getCoefs))
        .def("size", &Polynomial::size);
}
} // namespace vampyr
