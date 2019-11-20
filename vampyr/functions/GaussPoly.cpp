#include <array>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "functions/Gaussian.h"
#include "functions/GaussPoly.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void gauss_poly(py::module &m, py::class_<Gaussian<3>> &gaussian) {
    const auto D = 3;


py::class_<GaussPoly<D>>(m, "GaussPoly", gaussian)
        .def(py::init<double, double, const Coord<D> &, const std::array<int, D> &>())
        .def("setPoly", &GaussPoly<D>::setPoly)
        .def("differentiate", &GaussPoly<D>::differentiate)
        .def("calcOverlap", py::overload_cast<GaussPoly<D> &>(&GaussPoly<D>::calcOverlap))
        .def("evalf", py::overload_cast<const Coord<D> &>(&GaussPoly<D>::evalf, py::const_))
        .def("evalf", py::overload_cast<const double, int>(&GaussPoly<D>::evalf, py::const_));


}
} // namespace vampyr
